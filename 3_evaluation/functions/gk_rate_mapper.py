"""
gk_rate_mapper.py
=================
Maps Model RATE___xxx column names (from global_kin.tec) to their
corresponding reactions in the .dat chemistry file.

Key insight: The RATE index IS the sequential reaction number when counting
only non-comment, non-empty lines containing ">" after the "*" delimiter.

Usage:
    from gk_rate_mapper import build_reaction_lookup, get_rate_label

    rxn_map = build_reaction_lookup("path/to/global_kin.dat")
    print(rxn_map[859])   # -> ("CO2 + M > CO + O + M", "N01  Eremin1997,Nist2015")

    # Or get a short label for a plot legend
    label = get_rate_label(rxn_map, 859)
    # -> "CO2+M>CO+O+M"
"""

import re
from pathlib import Path


def build_reaction_lookup(dat_path: str) -> dict[int, tuple[str, str]]:
    """
    Parse a Model .dat file and return a dict mapping
        reaction_index (int) -> (reaction_string, comment_string)

    The reaction_index matches the RATE___xxx number in the .tec output.
    Counting starts at 1 and increments for every non-comment line
    containing ">" after the "*" section delimiter.
    """
    dat_path = Path(dat_path)
    reactions = {}
    rxn_counter = 0

    with open(dat_path, "r", errors="replace") as f:
        lines = f.readlines()

    # Find the "*" delimiter that separates species from reactions
    start_line = None
    for i, line in enumerate(lines):
        if line.strip() == "*":
            start_line = i + 1
            break

    if start_line is None:
        raise ValueError(f"No '*' delimiter found in {dat_path}")

    for line in lines[start_line:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("!") or stripped == "*":
            continue
        if ">" not in stripped:
            continue

        rxn_counter += 1

        # Split off the kinetics part after ":"
        rxn_part = stripped.split(":")[0].strip()

        # Extract comment after "%" (Model label field)
        comment = ""
        if "%" in stripped:
            comment = stripped.split("%", 1)[-1].strip()

        reactions[rxn_counter] = (rxn_part, comment)

    return reactions


def get_rate_label(
    rxn_map: dict,
    rate_index: int,
    max_len: int = 40,
    include_comment: bool = False,
) -> str:
    """
    Return a short string label for a given RATE index, suitable for plot legends.

    Parameters
    ----------
    rxn_map       : output of build_reaction_lookup()
    rate_index    : the integer from RATE___xxx
    max_len       : truncate reaction string to this length
    include_comment: if True, append the .dat comment tag (e.g. "N01")
    """
    if rate_index not in rxn_map:
        return f"RATE_{rate_index} [not found]"

    rxn, comment = rxn_map[rate_index]

    # Compact: remove spaces around "+" and ">"
    compact = re.sub(r"\s*\+\s*", "+", rxn)
    compact = re.sub(r"\s*>\s*", ">", compact)

    label = compact[:max_len]
    if include_comment and comment:
        # Take only the first token of the comment (e.g. "N01")
        tag = comment.split()[0]
        label = f"{label}  [{tag}]"
        #print(label)

    return label


def print_mapping(
    rxn_map: dict,
    rate_indices: list[int] | None = None,
) -> None:
    """
    Pretty-print the mapping for a list of RATE indices
    (or all entries if rate_indices is None).
    """
    indices = rate_indices if rate_indices is not None else sorted(rxn_map.keys())
    header = f"{'RATE col':<12} {'RXN #':<7} {'Reaction':<60} Comment"
    print(header)
    print("-" * len(header))
    for idx in indices:
        if idx in rxn_map:
            rxn, comment = rxn_map[idx]
            print(f"RATE_{idx:<7d} {idx:<7d} {rxn:<60} {comment}")
        else:
            print(f"RATE_{idx:<7d} {'':7} *** NOT FOUND ***")


# ---------------------------------------------------------------------------
# Example / quick-run
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    dat_file = sys.argv[1] if len(sys.argv) > 1 else "global_kin.dat"

    print(f"Building reaction lookup from: {dat_file}")
    rxn_map = build_reaction_lookup(dat_file)
    print(f"  Total reactions indexed: {len(rxn_map)}\n")

    # Example: the indices that appear in the .tec header for CRATE_PLOT='CO2'
    example_indices = [1, 28, 55, 82, 109, 136, 163, 190, 217, 243,
                       269, 270, 271, 272, 859, 886, 913]

    print("Example mapping (key CO2 reactions):")
    print_mapping(rxn_map, example_indices)

    print("\nShort labels for plot legends:")
    for idx in example_indices:
        print(f"  RATE_{idx}: {get_rate_label(rxn_map, idx, include_comment=True)}")


def count_species_occurrences(
    dat_path: str | Path,
    species_list: list[str],
) -> dict[str, dict[str, int]]:
    """
    Parameters
    ----------
    dat_path     : Pfad zum global_kin.dat
    species_list : Liste der Spezies-Namen (ohne M, TE, TGAS etc.)
 
    Returns
    -------
    dict mit { species: {'educt': int, 'product': int, 'total': int} }
    """
    dat_path = Path(dat_path)
 
    # --- Reaktionen einlesen ---
    with open(dat_path, "r", errors="replace") as f:
        lines = f.readlines()
 
    start_line = None
    for i, line in enumerate(lines):
        if line.strip() == "*":
            start_line = i + 1
            break
    if start_line is None:
        raise ValueError(f"Kein '*'-Delimiter in {dat_path} gefunden.")
 
    reaction_sides = []   # Liste von (left, right) Strings
    for line in lines[start_line:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("!") or stripped == "*":
            continue
        if ">" not in stripped:
            continue
        rxn_part = stripped.split(":")[0].strip()
        left, right = rxn_part.split(">", 1)
        reaction_sides.append((left, right))
 
    # --- Zählen ---
    def sp_in(sp: str, side: str) -> bool:
        # Lookahead/behind verhindert false positives (z.B. CO2 in CO2v01)
        pat = r"(?<![A-Za-z0-9])" + re.escape(sp) + r"(?![A-Za-z0-9^])"
        return bool(re.search(pat, side))
 
    result = {}
    for sp in species_list:
        n_educt = 0
        n_product = 0
        for left, right in reaction_sides:
            if sp_in(sp, left):
                n_educt += 1
            if sp_in(sp, right):
                n_product += 1
        result[sp] = {
            "educt":   n_educt,
            "product": n_product,
            "total":   n_educt + n_product,
        }
 
    return result

# def count_species_reactions(dat_path: str | Path, species_list: list[str]) -> dict[str, dict[str, int]]:
#     """
#     Count the number of reactions and the number of unique species occurrences in the dat-file.

#     Parameters
#     ----------
#     dat_path     : Path to the global_kin.dat file
#     species_list : List of species names (without M, TE, TGAS etc.)

#     Returns
#     -------
#     dict with { species: {'educt': int, 'product': int, 'total': int} }
#     """
#     dat_path = Path(dat_path)
#     reaction_count = 0

#     with open(dat_path, "r", errors="replace") as f:
#         lines = f.readlines()

#     for line in lines:
#         stripped = line.strip()
#         if stripped.startswith("!") or stripped == "*":
#             continue
#         reaction_count += 1

#     print(f"Total reactions: {reaction_count}")

#     return count_species_occurrences(dat_path, species_list)