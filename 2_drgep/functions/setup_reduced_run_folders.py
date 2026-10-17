import re
import shutil
from pathlib import Path

def _extract_numeric_key(path: Path) -> tuple:
    """
    Extracts all numbers (epsilon, delta, ...) from the filename for correct
    NUMERIC sorting.

    Note: sorted(list(src_dir.glob("*.dat"))) sorts Path objects
    lexicographically as strings. For filenames such as
    "model_0.5_0.0.dat" and "model_0.52_0.0.dat", this produces the wrong
    order because '_' (ASCII 95) is sorted lexicographically after digits such
    as '2', '4', '6', and '8' (ASCII 50-56):
        "...0.5_..."  ->  the character after '5' is '_'  (95)
        "...0.52_..." ->  the character after '5' is '2'  (50)
    => "0.5_..." is sorted after "0.52_...", "0.54_...", "0.56_...", "0.58_..."
    Numeric sorting by the extracted float values fixes this.
    """
    stem = path.stem.replace("model_", "")
    numbers = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", stem)
    return tuple(float(n) for n in numbers) if numbers else (float("inf"),)


def prepare_reduced_run_folders(source_dat_dir: str | Path, target_base_dir: str | Path):
    src_dir = Path(source_dat_dir)
    target_dir = Path(target_base_dir)

    if not src_dir.exists():
        raise FileNotFoundError(f"Quellordner nicht gefunden: {src_dir}")

    # Find all .dat files without sorted() (lexicographic); see
    # _extract_numeric_key() for the reasoning. Read them unsorted for now.
    dat_files = list(src_dir.glob("*.dat"))
    if not dat_files:
        print(f"Keine .dat-Dateien in {src_dir} gefunden.")
        return

    print(f"Gefundene Dateien ({len(dat_files)}): {[f.name for f in dat_files]}")

    run1_file = None
    other_files = []

    # Check specifically for the epsilon=0.0 / delta=0.0 file
    for f in dat_files:
        if "_0.0_0.0.dat" in f.name or f.name == "model_0.0.dat":
            run1_file = f
        else:
            other_files.append(f)

    # Sort NUMERICALLY (rather than lexicographically) using the numbers in the
    # filename (e.g. epsilon, delta)
    other_files.sort(key=_extract_numeric_key)

    # Combine: run 1 first, followed by the remaining files in numeric order
    ordered_files = []
    if run1_file:
        ordered_files.append(run1_file)
    else:
        print("WARNUNG: Keine explizite 0.0_0.0 Datei für Run 1 gefunden!")

    ordered_files.extend(other_files)

    print(f"Reihenfolge nach numerischer Sortierung: {[f.name for f in ordered_files]}")
    print(f"Starte Erstellung der Run-Ordner in: {target_dir.resolve()}\n")

    # Create folders and copy files
    for idx, dat_path in enumerate(ordered_files, start=1):
        # Extract the parameter suffix (e.g. "0.52_0.0" from "model_0.52_0.0.dat")
        param_suffix = dat_path.stem.replace("model_", "").replace("model", "").strip("_")
        if not param_suffix:
            param_suffix = "0.0_0.0"
        # NOTE (fix): the original
        #   param_suffix, _ = param_suffix.split("_", 1)
        # kept only the FIRST part (e.g. "0.52") for multiple parameters such
        # as "0.52_0.0" and silently discarded delta. This went unnoticed only
        # because delta had always been 0.0.
        # -> param_suffix is now kept in full (e.g. "0.52_0.0").

        param_suffix, _ = param_suffix.split("_", 1) if "_" in param_suffix else (param_suffix, "")
        folder_name = f"Run_{idx}_{param_suffix}"
        run_folder = target_dir / folder_name

        run_folder.mkdir(parents=True, exist_ok=True)
        target_file_path = run_folder / "model.dat"

        shutil.copy2(dat_path, target_file_path)
        print(f"  [+] Ordner '{folder_name}' erstellt -> 'model.dat' gespeichert.")

    print(f"\nFertig! {len(ordered_files)} Run-Ordner wurden erfolgreich vorbereitet.")