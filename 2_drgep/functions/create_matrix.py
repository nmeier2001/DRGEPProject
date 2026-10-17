import pickle as pkl
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import heapq


def create_matrix(labels_path: str) -> tuple[np.ndarray, dict, list, list]:
    """
    Loads reaction labels and builds stoichiometry matrix.

    Args:
        labels_path (str): Path to the pickle file containing reaction labels.

    Returns:
        nu:           np.ndarray [n_species × n_reactions]
        rxn_info:     dict with reaction information
        species_list: list of species in row order
        rxn_list:     list of reaction keys in column order
    """
    with open(labels_path, 'rb') as f:
        labels = pkl.load(f)

    rxn_info = {}
    for label, reaction in labels.items():
        reaction, group = reaction.split('  ')
        reactants_str, products_str = reaction.split('>')

        reactants_raw = [s.strip() for s in reactants_str.split('+')]
        products_raw  = [s.strip() for s in products_str.split('+')]

        # Build stoichiometry by counting duplicate entries
        reactants_dict = {}
        for r in reactants_raw:
            reactants_dict[r] = reactants_dict.get(r, 0) - 1

        products_dict = {}
        for p in products_raw:
            products_dict[p] = products_dict.get(p, 0) + 1

        rxn_info[label] = {
            'group':      group.strip(),
            'reaction':   reaction,
            'reactants':  list(reactants_dict.keys()),
            'stoch_reac': list(reactants_dict.values()), 
            'products':   list(products_dict.keys()),    
            'stoch_prod': list(products_dict.values()),  
        }

    # Ordered lists for matrix indices
    rxn_list    = list(rxn_info.keys())
    rxn_idx     = {label: i for i, label in enumerate(rxn_list)}

    species_set  = set()
    for info in rxn_info.values():
        species_set.update(info['reactants'])
        species_set.update(info['products'])
    species_list = sorted(species_set)
    species_idx  = {s: i for i, s in enumerate(species_list)}

    # Build the matrix as a NumPy array
    nu = np.zeros((len(species_list), len(rxn_list)), dtype=float)

    for label, info in rxn_info.items():
        j = rxn_idx[label]
        for reactant, stoch in zip(info['reactants'], info['stoch_reac']):
            nu[species_idx[reactant], j] = stoch
        for product, stoch in zip(info['products'], info['stoch_prod']):
            nu[species_idx[product], j] = stoch

    print(f"Matrix created with shape: "
          f"{len(species_list)} species x {len(rxn_list)} reactions")

    return nu, rxn_info, species_list, rxn_list

def load_rate_vector(study_path: str, time: float, run: int) -> np.ndarray:
    """
    Load the rate vector from a pickle file.

    Args:
        study_path (str): Path to the pickle file containing the rate vector.
        time (float): The time at which to load the rate vector.
        run (int): The run number.

    Returns:
        np.ndarray: Rate vector, ordered by reaction keys.
    """
    with open(study_path, 'rb') as f:
        study = pkl.load(f)

    tec = study['runs'][run]['tec']
    time_series = tec['T (s)']

    idx_time = (time_series - time).abs().argmin()

    rate_cols = [col for col in tec.columns if col.startswith('RATE_')]
    rate_vector = tec.iloc[idx_time][rate_cols].values

    return rate_vector, rate_cols


# Start of the functions for DRGEP
def compute_rAB(nu: np.ndarray, omega: np.ndarray, species_list: list, tp_name: str) -> np.ndarray:
    """
    Computes the r_AB coupling matrix for a single time point.

    Args:
        nu:    np.ndarray [n_species × n_reactions]  Stoichiometry matrix
        omega: np.ndarray [n_reactions]              Rate vector at time t
        species_list: list of species in row order
        tp_name: name of the time point

    Returns:
        r_AB:  np.ndarray [n_species × n_species]   Coupling coefficients
    """
    # Step 1: |ν| · |ω| → [n_species × n_reactions]
    weighted = np.abs(nu) * np.abs(omega)[np.newaxis, :]

    # Step 2: denominator → total turnover per species [n_species]
    denom = weighted.sum(axis=1)
    denom_safe = np.where(denom > 0, denom, 1.0)  # Avoid division by zero

    # Step 3: B_mask → 1 if species B occurs in reaction i [n_species × n_reactions]
    B_mask = (nu != 0).astype(float)

    # Step 4: r_AB → [n_species × n_species]
    # weighted @ B_mask.T → numerator: turnover of A in reactions that also contain B
    r_AB = (weighted @ B_mask.T) / denom_safe[:, np.newaxis]

    return r_AB

def compute_generalized_R(r_primary: np.ndarray, species_list: list, 
                           targets: list) -> np.ndarray:
    """
    True DRGEP step: generalized coefficient R_{T,S} via
    widest-path search (max-product) starting from each target species.
    """
    n = len(species_list)
    species_idx = {s: i for i, s in enumerate(species_list)}
    R = np.zeros(n)  # Final maximum over all targets

    for t in targets:
        ti = species_idx[t]
        best = np.zeros(n)
        best[ti] = 1.0
        visited = np.zeros(n, dtype=bool)
        heap = [(-1.0, ti)]  # Max-heap using negated values

        while heap:
            neg_val, i = heapq.heappop(heap)
            val = -neg_val
            if visited[i]:
                continue
            visited[i] = True
            for j in range(n):
                if visited[j]:
                    continue
                candidate = val * r_primary[i, j]
                if candidate > best[j]:
                    best[j] = candidate
                    heapq.heappush(heap, (-candidate, j))

        R = np.maximum(R, best)

    return R

def compute_RTS(nu: np.ndarray, dic_rates: dict, species_list: list, 
                 targets: list) -> np.ndarray:
    """
    Calculates R_T,S: first compute primary r_AB for each time point,
    then immediately run the path search (compute_generalized_R),
    and only then take the maximum over time points.
    """
    n_species = len(species_list)
    R_T_S = np.zeros(n_species)

    for tp_name, omega in dic_rates.items():
        r_AB_t = compute_rAB(nu, omega[0], species_list, tp_name)

        # Run the path search now for each time point, not only at the end
        R_T_S_t = compute_generalized_R(r_AB_t, species_list, targets)

        # Take the maximum over time points after the path search
        R_T_S = np.maximum(R_T_S, R_T_S_t)

        print(f"{tp_name:35s} | aktive Reaktionen: {(np.abs(omega[0]) > 1e-30).sum()}")

    return R_T_S


def select_species(R_T_S: np.ndarray, species_list: list, targets: list, epsilon: float) -> list:
    kept = [s for s, r in zip(species_list, R_T_S) if r >= epsilon]
    return kept

def filter_reactions(rxn_info: dict, kept_species: list) -> dict:
    """
    Removes reactions containing eliminated species.

    Args:
        rxn_info:      full reaction dictionary
        kept_species:  species to keep from select_species()

    Returns:
        skeletal_rxn_info: reduced reaction dictionary
    """
    kept_set = set(kept_species)
    skeletal_rxn_info = {}

    for label, info in rxn_info.items():
        all_species = set(info['reactants'] + info['products'])
        
        if all_species.issubset(kept_set):
            skeletal_rxn_info[label] = info

    print(f"\nReaktionen: {len(rxn_info)} → {len(skeletal_rxn_info)} "
          f"({len(rxn_info) - len(skeletal_rxn_info)} entfernt)")

    return skeletal_rxn_info

def filter_reactions_with_sigma(nu: np.ndarray,
                                 rxn_info: dict,
                                 kept_species: list,
                                 species_list: list,
                                 dic_rates: dict,
                                 sigma: float = 0.01) -> dict:

    kept_set    = set(kept_species)
    species_idx = {s: i for i, s in enumerate(species_list)}
    rxn_list    = list(rxn_info.keys())
    rxn_idx     = {label: j for j, label in enumerate(rxn_list)}

    # --- Stage 1: species check ---
    stage1_rxns = {
        label: info for label, info in rxn_info.items()
        if set(info['reactants'] + info['products']) - {'M'} <= kept_set
    }
    stage1_labels  = list(stage1_rxns.keys())
    stage1_indices = [rxn_idx[l] for l in stage1_labels]

    print(f"Stufe 1: {len(rxn_info)} → {len(stage1_rxns)} Reaktionen")

    if sigma == 0:
        return stage1_rxns

    kept_rows  = [species_idx[s] for s in kept_species if s in species_idx]
    n_kept     = len(kept_rows)
    n_stage1   = len(stage1_labels)
    kept_list  = [kept_species[i] for i in range(len(kept_species))
                  if kept_species[i] in species_idx]

    nu_net = np.zeros((n_kept, n_stage1))
    sp_to_row = {s: i for i, s in enumerate(kept_list)}

    for j, label in enumerate(stage1_labels):
        info = stage1_rxns[label]

        # Reactants: negative
        for sp, stoch in zip(info['reactants'], info['stoch_reac']):
            if sp in sp_to_row:
                nu_net[sp_to_row[sp], j] += stoch  # stoch_reac is already negative

        # Products: positive
        for sp, stoch in zip(info['products'], info['stoch_prod']):
            if sp in sp_to_row:
                nu_net[sp_to_row[sp], j] += stoch  # stoch_prod is positive

    print(f"nu_net nonzero: {np.count_nonzero(nu_net)}")
    print(f"nu_net min/max: {nu_net.min():.2f} / {nu_net.max():.2f}")

    # Sanity check: Katalytische Reaktionen (E+CO2>E+CO2) sollten nu_net=0 haben
    cat_count = np.sum(np.all(nu_net == 0, axis=0))
    print(f"Katalytische Reaktionen (nu_net=0): {cat_count} "
          f"→ werden automatisch entfernt wenn sigma > 0")

    max_importance = np.zeros(n_stage1)

    for tp_name, (rate_vector, rate_cols) in dic_rates.items():

        rate_map = {col: val for col, val in zip(rate_cols, rate_vector)}

        omega_sub = np.array([
            abs(rate_map.get(label, 0.0))
            for label in stage1_labels
        ])

        nonzero = np.count_nonzero(omega_sub)
        if nonzero == 0:
            print(f"  ⚠ Zeitpunkt {tp_name}: alle Raten = 0, überspringe")
            continue

        # P_ij: production contribution (nu_net > 0)
        # C_ij: consumption contribution (nu_net < 0)
        prod_contrib = np.maximum(0,  nu_net) * omega_sub[np.newaxis, :]
        cons_contrib = np.maximum(0, -nu_net) * omega_sub[np.newaxis, :]

        prod_total = prod_contrib.sum(axis=1)
        cons_total = cons_contrib.sum(axis=1)

        prod_total_s = np.where(prod_total > 0, prod_total, 1.0)
        cons_total_s = np.where(cons_total > 0, cons_total, 1.0)

        P_ij = prod_contrib / prod_total_s[:, np.newaxis]
        C_ij = cons_contrib / cons_total_s[:, np.newaxis]

        importance_t   = np.maximum(P_ij, C_ij).max(axis=0)
        max_importance = np.maximum(max_importance, importance_t)

        print(f"  Zeitpunkt {tp_name}: "
              f"omega nonzero={nonzero}, "
              f"max importance={importance_t.max():.4f}")

    keep_mask = max_importance >= sigma
    skeletal_rxn_info = {
        label: stage1_rxns[label]
        for j, label in enumerate(stage1_labels)
        if keep_mask[j]
    }

    print(f"Stufe 2 (σ={sigma:.3f}): "
          f"{len(stage1_rxns)} → {len(skeletal_rxn_info)} "
          f"({len(stage1_rxns) - len(skeletal_rxn_info)} entfernt)")

    return skeletal_rxn_info

def plot_RAB_fraction_vs_epsilon(
        study_path, labels_path, run_indices,
        timepoint_fn, targets,
        epsilon_list=None, save_dir=None):
    """
    R_TS = compute_RTS(...) is already the fully generalized
    R_{T,S} vector per species (maximum over all paths, targets,
    and time points). No additional BFS or indexing step is needed.
    """
    plt.rcParams.update({'font.size': 14})
    if epsilon_list is None:
        epsilon_list = [0.05, 0.1, 0.3]
    save_dir = Path(save_dir) if save_dir is not None else Path('.')

    nu, rxn_info, species_list, rxn_list = create_matrix(labels_path)
    species_idx      = {s: i for i, s in enumerate(species_list)}
    target_indices   = [species_idx[t] for t in targets if t in species_idx]
    non_target_mask  = np.ones(len(species_list), dtype=bool)
    non_target_mask[target_indices] = False

    with open(study_path, 'rb') as f:
        study = pkl.load(f)

    n_runs = len(run_indices)
    cmap   = plt.cm.plasma
    colors = [cmap(i / max(n_runs - 1, 1)) for i in range(n_runs)]

    n_cols = min(3, n_runs)
    n_rows = int(np.ceil(n_runs / n_cols))
    fig_hist, axes_hist = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    axes_hist = np.atleast_1d(axes_hist).flatten()

    fig_frac, ax_frac = plt.subplots(figsize=(8, 6))
    eps_continuous = np.linspace(0.0, 1.0, 50)

    for plot_idx, run_idx in enumerate(run_indices):
        run     = study['runs'][run_idx]
        powplug = run['param_value']
        tec     = run['tec']
        color   = colors[plot_idx]

        print(f"\nRun {run_idx} (P = {powplug}):")

        timepoints = timepoint_fn(tec['T (s)'], tec['EN_Td'])
        dic_rates = {
            tp_idx: load_rate_vector(study_path, tp['time'], run_idx)
            for tp_idx, tp in timepoints.items()
        }

        # 1D vector: already-computed R_{T,S} value for each species
        R_TS = compute_RTS(nu, dic_rates, species_list, targets)
        R_TS[target_indices] = 1.0  # Always keep the targets
        R_TS_others = R_TS[non_target_mask]

        # --- Histogramm pro Run ---
        ax_h = axes_hist[plot_idx]
        ax_h.hist(R_TS_others, bins=30, color=color, alpha=0.8, edgecolor='black')
        ax_h.axvline(epsilon_list[0], color='red', linestyle='--', linewidth=1.2)
        ax_h.set_title(f'Run {run_idx} (P={powplug})')
        ax_h.set_xlabel('R_{T,S} (generalisiert)')
        ax_h.set_ylabel('Anzahl Spezies')
        ax_h.set_xlim(0, 1)
        ax_h.grid(True, alpha=0.3)

        # --- Fraction-vs-epsilon ---
        kept_counts = [(R_TS >= e).sum() for e in eps_continuous]
        ax_frac.plot(eps_continuous, kept_counts, marker='o', markersize=3,
                     color=color, label=f'{powplug} W')

        for eps in epsilon_list:
            kept_species = [s for s, r in zip(species_list, R_TS) if r >= eps]
            print(f"  ε={eps:.3f}: {len(kept_species)} Spezies behalten")

    for ax in axes_hist[n_runs:]:
        ax.axis('off')
    fig_hist.suptitle(f'Generalisierte R_{{T,S}}-Verteilung | Targets: {targets}', fontsize=12)
    fig_hist.tight_layout()
    fig_hist.savefig(save_dir / f'R_TS_histogram_per_run_{len(targets)}targets.png', dpi=300)
    plt.close(fig_hist)

    for eps in epsilon_list:
        ax_frac.axvline(eps, color='gray', linestyle=':', linewidth=1.2, alpha=0.8)
        ax_frac.text(eps + 0.003, ax_frac.get_ylim()[1] * 0.98, f'ε={eps}',
                     color='gray', fontsize=8, va='top')
    ax_frac.set_xlim(-0.05, 1.05)
    ax_frac.set_xlabel('ε (Schwellwert)')
    ax_frac.set_ylabel('Behaltene Spezies')
    ax_frac.grid(True, alpha=0.3)
    ax_frac.legend(fontsize=12)
    fig_frac.tight_layout()
    fig_frac.savefig(save_dir / f'R_TS_fraction_vs_epsilon_{len(targets)}targets.png', dpi=300)
    plt.close(fig_frac)

    print("\nAll plots saved.")