from dataclasses import dataclass
import os
from pathlib import Path
from typing import Optional
from scipy.interpolate import interp1d
from dataclasses import dataclass

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import math
import re

@dataclass
class DeviationAnalyzer:
    study: dict
    ref_param_value: float
    value_col: str

def calculate_SEI(nam):
    pow_list = nam.get('POWPLUG', [0])          # W/cm³ (Liste)
    if not isinstance(pow_list, list):
        pow_list = [pow_list]
    P_peak = max(pow_list)                    # Peak power density [W/cm³]    
    
    t_res   = nam.get('XPOW', 15e-9)
    if not isinstance(t_res, list):
        t_res = [t_res]
    t_res = t_res[2] - t_res[1]  # Verweilzeit [s]
    vol     = nam.get('VOL_PLASMA', 1.0)        # plasma volume [cm³]
    T_gas   = nam.get('TGAS0', 300.0)           # gas temperature [K]
    press   = nam.get('PRESS0', 1.0)            # pressure [atm]
    
    # Pressure in Pa (positiv = atm, negativ = Torr laut .nam-Kommentar)
    P_Pa = press * 101325 if press > 0 else abs(press) * 133.322
    
    # Gasdichte
    k_B = 1.380649e-23
    eV  = 1.602176634e-19
    n0  = P_Pa / (k_B * T_gas)   # [m⁻³]
    n0_cm3 = n0 * 1e-6            # [cm⁻³]
    
    # Total energie in plasma
    E_total_J = P_peak * vol * t_res   # [J]

    # SCCM → L/s
    flow_sccm = nam.get('FLOW_SCCMTOT', 0)
    flow_cm3_s = flow_sccm / 60.0
    flow_liter_min = flow_sccm * 1e-3  # [L/min]
    flow_liter_sec = flow_liter_min / 60   # [L/s]
    
    # Particle flux density [cm⁻³ s⁻¹]
    Q_n = (flow_cm3_s * n0_cm3) / vol
    
    
    # SEI per Molekül
    power_W = P_peak * vol  # [W] = [J/s]
    SEI_eV = (power_W * 6.24e18 * 24.6) / (flow_liter_sec * 6.022e23)  # [eV/molekül]

    return SEI_eV

def calculate_energy_efficiency(n_CO2_0, n_CO2_final, SEI_eV):
    X = 1 - (n_CO2_final) / n_CO2_0 if n_CO2_0 != 0 else 0
    n = X * (2.9/SEI_eV)

    if n_CO2_0 == n_CO2_final:
        n = 'CO2_no_change'
    elif n_CO2_final > n_CO2_0:
        n = 'CO2_increase'
    
    return n

# calculate the ratio of diffrent species to CO2
def calculate_ratio(study, species: list) -> None:
    """
    Berechnet Molekül/CO2 Ratios und hängt sie direkt an den tec DataFrame an.
    """
    for idx, run in study['runs'].items():
        df: pd.DataFrame = run['tec']
        new_cols = {}

        for molecule in species:
            if molecule in df.columns and 'CO2' in df.columns:
                new_cols[f'{molecule}_to_CO2_ratio'] = df[molecule] / df['CO2']
            else:
                print(f"  ⚠ Run {idx}: {molecule} nicht gefunden, überspringe.")

        if new_cols:
            run['tec'] = pd.concat([df, pd.DataFrame(new_cols, index=df.index)], axis=1)


def parse_boltz_EN_table(boltz_file):
    """
    Liest E/N [V·cm²] und Te [eV] Lookup-Tabelle aus global_kin_boltz.out.
    Gibt EN_Td [Td] und Te_eV [eV] als Arrays zurück.
    """
    EN_list = []
    Te_list = []

    with open(boltz_file, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if 'ELECTRIC FIELD/NUMBER DENSITY=' in line:
            EN_vcm2 = float(line.split('=')[1].strip().split()[0])
            EN_list.append(EN_vcm2)
            # Suche nächste Te-Zeile
        if 'EQUIVALENT ELECTRON TEMPERATURE=' in lines[i]:
            Te = float(lines[i].split('=')[1].strip().split()[0])
            Te_list.append(Te)

   
    EN_arr = np.array(EN_list)
    Te_arr = np.array(Te_list)
    if len(EN_arr) == 0 or len(Te_arr) == 0:
        raise ValueError("Keine E/N oder Te Daten in der Boltzmann-Datei gefunden.")

    # Umrechnung: V·cm² → Td  (1 Td = 1e-17 V·cm²)
    EN_Td = EN_arr / 1e-17

    return EN_Td, Te_arr


def add_EN_to_df(df, boltz_file):
    """
    Fügt E/N [Td] als neue Spalte zu self.df hinzu,
    interpoliert aus der Boltzmann-Lookup-Tabelle via TE.
    """
    EN_Td, Te_table = parse_boltz_EN_table(boltz_file)

    # Te muss monoton steigend sein für np.interp
    # Sortiere nach Te (für den Fall dass Einträge nicht sortiert sind)
    sort_idx = np.argsort(Te_table)
    Te_sorted = Te_table[sort_idx]
    EN_sorted = EN_Td[sort_idx]

    Te_sim = df['TE'].values
    EN_sim = np.interp(Te_sim, Te_sorted, EN_sorted)

    # pd.concat statt direkter Zuweisung → kein Fragmentieren
    df = pd.concat([df, pd.Series(EN_sim, name='EN_Td', index=df.index)], axis=1)
    return df

def add_EN_to_study(study: dict) -> None:
    """
    Fügt E/N [Td] als Spalte 'EN_Td' an jeden tec DataFrame an.
    Liest die Boltzmann-Lookup-Tabelle aus dem jeweiligen Run-Ordner.
    """
    for idx, run in study['runs'].items():
        #print(run['path'])
        boltz_path = run['path'] / 'global_kin_boltz.out'

        if not boltz_path.exists():
            #print(boltz_path)
            print(f"  ⚠ Run {idx}: global_kin_boltz.out nicht gefunden, überspringe.")
            continue

        run['tec'] = add_EN_to_df(run['tec'], boltz_path)


def plot_xpow_size(xpow_values, tend=2e-4):
    """
    Plottet die Dreieckspulse für verschiedene XPOW Flankenbreiten.
    
    :param xpow_values: Liste von Flankenbreiten (Peak-Zeiten)
    :param tend: Endzeit der Simulation (feste rechte Grenze)
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    for flank in xpow_values:
        peak = flank
        end  = 2 * flank
        # Dreieckspuls: 0 → peak → end → tend
        t = [0, peak, end, tend]
        p = [0, 1,    0,   0   ]
        ax.plot(t, p, label=f'XPOW={flank*1e9:.0f} ns')

    ax.set_xlabel('Zeit (s)')
    ax.set_ylabel('Normierte Leistung [-]')
    ax.set_title('XPOW Dreieckspulse')
    ax.set_xscale('log')
    ax.set_xlim(1e-9, tend)
    ax.legend()
    ax.grid(True)
    plt.tight_layout()

def plot_xpow_from_study(study, path_run: Path):
    """
    Plottet die Dreieckspulse aller Runs aus einer geladenen Study.
    
    :param study: geladene Study-Struktur (dict)
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    for run_idx, run in study['runs'].items():
        nam = run['nam']
        if nam is None:
            print(f"  WARNUNG: Run {run_idx} hat keine .nam Datei, übersprungen.")
            continue

        xpow    = nam['XPOW']     # [x0, peak, end, tend, tcycle]
        powplug = nam['POWPLUG']  # [0.0, peak_power, 0.0, 0.0]

        peak    = xpow[1]
        power   = max(powplug)

        ax.plot(xpow, powplug, label=f'Run {run_idx}: Flanke={peak*1e9:.0f} ns, P={power:.2e} W/cm³')

    ax.set_xlabel('Zeit (s)')
    ax.set_ylabel('Leistung (W/cm³)')
    ax.set_title(f'XPOW Pulse – {study["param_name"]} Studie ({study["date"]})')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-9, max(run['nam']['XPOW'][3] for run in study['runs'].values() if run['nam'] is not None))
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(path_run / 'XPOW_vs_Time.png', dpi=300)

def plot_rates(run, path_labels, is_rate_study=False, plot_rates_flag=True):
    """
    Plottet die Reaktionsraten aus einem Run.
    Jede Rate bekommt ein eigenes Diagramm; je 4 Diagramme werden
    in einer 2×2-Abbildung zusammengefasst.
 
    :param run: einzelner Run (dict)
    :param path_labels: Pfad zur Pickle-Datei mit Reaktionslabels
    :param is_rate_study: Flag, ob es sich um eine Rate-Studie handelt
    """
    run_dir = Path(str(run['path'])) / 'Rates'   # str() erzwingt saubere Konvertierung
    run_dir.mkdir(parents=True, exist_ok=True)
    nam = run['nam']
    df = run['tec']
    if df is None:
        print(f"  WARNUNG: Run {run} hat keine tec Daten, übersprungen.")
        return
 
    with open(path_labels, 'rb') as f:
        labels = pickle.load(f)
 
    rate_cols = [c for c in df.columns if c.startswith("RATE_")]
    t = df['T (s)']
 
    # Je 4 Raten auf eine 2×2-Seite
    page_size = 4
    n_pages = math.ceil(len(rate_cols) / page_size)

    if is_rate_study:
        analyse_rate_contributions2(run_dir, run['tec'], nam, labels, species='CO2', reac_contribution=nam['ISIDE_PLOT'], threshold=0.00001)
    else:
        analyse_rate_contributions(run_dir, run['tec'], nam, labels, species='CO2', reac_contribution=nam['ISIDE_PLOT'], threshold=0.00001)
 
    if plot_rates_flag:
        for page_idx in range(n_pages):
            cols_on_page = rate_cols[page_idx * page_size : (page_idx + 1) * page_size]
 
            fig, axes = plt.subplots(2, 2, figsize=(14, 9))
            fig.suptitle(f'Reaction Rates – Run {run["param_value"]}  '
                         f'(Page {page_idx + 1}/{n_pages})', fontsize=13)
            axes_flat = axes.flatten()
 
            for ax_idx, col in enumerate(cols_on_page):
                ax = axes_flat[ax_idx]
                ax.plot(t, df[col])
                ax.set_xscale('log')
                ax.set_yscale('log')
                ax.set_xlabel('Time (s)')
                ax.set_ylabel('Rate [m⁻³ s⁻¹]')
                ax.set_title(labels.get(col, col), fontsize=9)
                ax.grid(True, which='both', ls='--', alpha=0.5)
 
            # Leere Subplots ausblenden (letzte Seite ggf. nicht voll)
            for ax_idx in range(len(cols_on_page), page_size):
                axes_flat[ax_idx].set_visible(False)
 
            plt.tight_layout()
            fname = run_dir / f'Rates_Run_{run["param_value"]}_p{page_idx + 1:03d}.png'
            plt.savefig(fname, dpi=300)
            plt.close(fig)

def analyse_rate_contributions(run_dir, df, nam, rxn_map, species, reac_contribution, threshold=0.00001):
    """
    Berechnet den relativen Beitrag jeder Rate zur Produktion/Vernichtung
    einer Spezies, integriert über die Zeit.
    threshold: Mindestanteil (z.B. 0.01 = 1%) um als relevant zu gelten
    """
    import numpy as np
    from scipy.integrate import trapezoid

    #print(run_dir)

    t = df['T (s)'].values
    rate_cols = [c for c in df.columns if c.startswith("RATE_")]
    #print(rxn_map)
    integrals = {}
    for col in rate_cols:
        #idx = int(col.split("_")[-1])
        rxn_str = rxn_map.get(col, (""))
        #print(f"  Rate {col}: Reaction={rxn_str}")
        # if species not in rxn_str:
        #     continue
        integral = trapezoid(np.abs(df[col].values), t)
        integrals[col] = (integral, rxn_str)
        #print(f"  {col}: Integral={integral:.3e}, Reaction={rxn_str}")
    #print(integrals)
    total = sum(v[0] for v in integrals.values())

    results = {
        col: {
            'integral': val,
            'fraction': val / total,
            'reaction': rxn
        }
        for col, (val, rxn) in integrals.items()
        #if val / total >= threshold
    }

    #speichere ergebnisse in txt datei
    with open(run_dir / 'rate_contributions.txt', 'w') as f:
        f.write(f"Rate Contributions for {nam['CRATE_PLOT']}; (Part in reaction:{reac_contribution}) (threshold={threshold*100:.1f}%):\n")
        f.write(f"{'Rate Column':<15} {'Integral':<15} {'Fraction':<10} {'Reaction'}\n")
        f.write("-" * 80 + "\n")
        for col, info in sorted(results.items(), key=lambda x: -x[1]['fraction']):
            f.write(f"{col:<15} {info['integral']:<15.3e} {info['fraction']:<10.3%} {info['reaction']}\n")

    # Sortiert nach Beitrag
    return dict(sorted(results.items(), key=lambda x: -x[1]['fraction']))

def analyse_rate_contributions2(run_dir, df, nam, rxn_map, species, reac_contribution, threshold):
    """
    Berechnet den relativen Beitrag jeder Rate zur Produktion/Vernichtung
    einer Spezies, integriert über die Zeit.
    threshold: Mindestanteil (z.B. 0.01 = 1%) um als relevant zu gelten
    """
    import numpy as np
    from scipy.integrate import trapezoid
    import re

    def get_edukte(rxn_str):
        rxn_clean = rxn_str.split('[')[0]
        parts = rxn_clean.split('>')
        return parts[0].upper() if parts else ''

    def species_is_educt(rxn_str, species):
        species_clean = species.strip('"\'')
        lhs = get_edukte(rxn_str)
        pattern = re.escape(species_clean.upper()) + r'[A-Z0-9]*'
        return bool(re.search(pattern, lhs))

    def is_co2_loss(rxn_str, species):
        """
        True nur wenn die Gesamtzahl der species-Moleküle abnimmt.
        species kommt als Edukt vor aber NICHT als Produkt.
        """
        species_clean = species.strip('"\'').upper()
        rxn_clean = rxn_str.split('[')[0]
        parts = rxn_clean.split('>')
        if len(parts) < 2:
            return False
        lhs = parts[0].upper()
        rhs = parts[1].upper()
        educt_pattern = re.escape(species_clean) + r'[A-Z0-9]*'
        if not re.search(educt_pattern, lhs):
            return False
        if re.search(educt_pattern, rhs):
            return False
        return True

    def classify_reaction(rxn_str, species):
        """
        Übernimmt die Klassifizierung 1:1 aus analyse_rates.py, statt einer
        eigenen (davon abweichenden) Logik zu folgen. Import bewusst HIER
        innerhalb der Funktion (lazy import), nicht am Dateikopf von
        functions_analysis.py - sonst entsteht ein Zirkelimport, da
        analyse_rates.py seinerseits `from functions_analysis import
        load_study_cache` importiert.

        Hinweis: `species` wird nur noch der Signatur zuliebe entgegen-
        genommen (Aufrufstellen unten übergeben es weiterhin), aber nicht
        mehr verwendet - die Koelman/Kozák-Klassifizierung hängt rein vom
        Reaktions-Label ab (X/E/N/I/V), nicht von der gerade betrachteten
        Zielspezies.
        """
        from functions.analyse_rates import classify_species_process
        return classify_species_process(rxn_str)

    t = df['T (s)'].values
    rate_cols = [c for c in df.columns if c.startswith("RATE_")]

    # --- Integrale für alle Reaktionen ---
    integrals = {}
    for col in rate_cols:
        rxn_str = rxn_map.get(col, "")
        integral = trapezoid(np.abs(df[col].values), t)
        integrals[col] = (integral, rxn_str)

    total_all = sum(v[0] for v in integrals.values())

    # --- Einzelreaktionen klassifizieren ---
    results = {}
    for col, (val, rxn) in integrals.items():
        loss = is_co2_loss(rxn, species)
        results[col] = {
            'integral':          val,
            'fraction_all':      val / total_all,      # Anteil am Gesamtaufkommen
            'fraction_co2_loss': 0.0,                  # wird unten befüllt
            'reaction':          rxn,
            'group':             classify_reaction(rxn, species) if loss else '—',
            'is_loss':           loss,
            'relevant':          species_is_educt(rxn, species),
        }

    # --- CO2-Verlust Nenner berechnen ---
    total_co2_loss = sum(
        info['integral']
        for info in results.values()
        if info['is_loss']
    )

    # fraction_co2_loss befüllen
    count_loss_reactions = sum(1 for info in results.values() if info['is_loss'])
    for info in results.values():
        if info['is_loss'] and total_co2_loss > 0:
            info['fraction_co2_loss'] = info['integral'] / total_co2_loss

    # --- Gruppenanalyse (nur echte CO2-Verlustprozesse) ---
    # GROUP_NAMES kommt jetzt ebenfalls aus analyse_rates.py, damit hier
    # exakt dieselben fünf Koelman/Kozák-Kategorien (X/E/N/I/V) verwendet
    # werden wie dort - vorher waren es vier eigene, davon abweichende
    # Gruppennamen (u.a. wurden X und E hier zu "Elektronenstoß"
    # zusammengefasst, was die beiden physikalisch unterschiedlichen
    # Prozesse Ionisation/Anregung und Anlagerung/Rekombination vermischt
    # hätte).
    from functions.analyse_rates import GROUP_NAMES
    group_names = list(GROUP_NAMES)
    group_totals    = {g: 0.0 for g in group_names}
    group_fractions = {g: 0.0 for g in group_names}

    for info in results.values():
        if info['is_loss'] and info['group'] in group_totals:
            group_totals[info['group']] += info['integral']

    if total_co2_loss > 0:
        group_fractions = {
            g: group_totals[g] / total_co2_loss
            for g in group_names
        }

    # --- Output-Datei ---
    with open(run_dir / 'rate_contributions.txt', 'w') as f:
        count_rates = len(results)

        # Kopfzeile
        f.write(f"Rate Contributions for {species}     "
                f"(Part in reaction: {reac_contribution})     "
                f"(threshold={threshold*100:.3f}%)    "
                f"Loss reactions: {count_loss_reactions}     "
                f"Total number of rate columns: {count_rates}\n")
        f.write("=" * 110 + "\n\n")

        # Gruppenübersicht
        f.write("GRUPPENANALYSE — Anteil am CO2-Gesamtverlust (Kozak/Koelmann)\n")
        f.write(f"  Gesamt-CO2-Verlust (integriert): {total_co2_loss:.3e}\n\n")
        f.write(f"  {'Gruppe':<50} {'Integral':<15} {'Anteil'}\n")
        f.write("  " + "-" * 55 + "\n")
        for g in sorted(group_names, key=lambda g: -group_fractions[g]):
            f.write(f"  {g:<50} {group_totals[g]:<15.3e} {group_fractions[g]:.2%}\n")

        f.write("\n" + "=" * 110 + "\n\n")

        # Einzelreaktionen — nur CO2-Verlustprozesse über threshold
        f.write("CO2-VERLUSTPROZESSE (Einzelreaktionen)\n")
        f.write(f"  {'Rate Column':<15} {'Gruppe':<50} "
                f"{'Anteil CO2-Verlust':<22} {'Reaktion'}\n")
        f.write("  " + "-" * 105 + "\n")
        for col, info in sorted(results.items(),
                                key=lambda x: -x[1]['fraction_co2_loss']):
            if not info['is_loss']:
                continue
            if info['fraction_co2_loss'] < threshold:
                continue
            f.write(f"  {col:<15} {info['group']:<50} "
                    f"{info['fraction_co2_loss']:<22.3%} {info['reaction']}\n")

        f.write("\n" + "=" * 110 + "\n\n")

        # Alle Reaktionen (wie bisher, für Vollständigkeit)
        f.write("ALLE REAKTIONEN (Anteil am Gesamtaufkommen)\n")
        f.write(f"  {'Rate Column':<15} {'Integral':<15} {'Anteil (gesamt)':<18} "
                f"{'Reaktion':<35} {'Gruppe':<20} {'Loss'}\n")
        f.write("  " + "-" * 120 + "\n")
        for col, info in sorted(results.items(),
                                key=lambda x: -x[1]['fraction_all']):
            if info['fraction_all'] < threshold:
                continue
            f.write(f"  {col:<15} {info['integral']:<15.3e} "
                    f"{info['fraction_all']:<18.3%} {info['reaction']:<35} "
                    f"{info['group']:<20} {info['is_loss']}\n")

    return dict(sorted(results.items(),
                       key=lambda x: -x[1]['fraction_co2_loss'])), group_fractions

def load_study_cache(path: Path) -> dict:
    with open(path, 'rb') as f:
        study = pickle.load(f)
    # print(f"   Study geladen: {path}")
    return study

def check_file_exists(file_path: str) -> bool:
    """
    Check if a file exists at the given path.
    
    Args:
        file_path (str): Path to the file.
    
    Returns:
        bool: True if file exists, False otherwise.
    """
    # Validate input type
    if not isinstance(file_path, str) or not file_path.strip():
        raise ValueError("file_path must be a non-empty string.")

    # Check if the path exists and is a file
    return os.path.isfile(file_path)

# creat a function to show the progress of the analysis
def progress_bar(current, total, bar_length=40):
    """
    Displays a progress bar in the console.
    
    Args:
        current (int): Current progress count.
        total (int): Total count for completion.
        bar_length (int): Length of the progress bar in characters.
    """
    if total <= 0:
        raise ValueError("Total must be greater than zero.")
    
    fraction = current / total
    filled_length = int(bar_length * fraction)
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    percent = fraction * 100
    print(f'\r  Progress: |{bar}| {percent:.2f}% Complete', end='\r')
    
    if current == total:
        print()  # Move to the next line on completion


def calculat_E_diffrence(study):
    '''
    Berechnet die Differenze zwischen Elektronendichte vom vollständigen Mechanismus und dem reduzierten Mechanismus.
        Args:
            study (dict): Geladene Study-Struktur (dict)
    '''
    # 1) Referenz (voller Mechanismus) und reduzierte Runs laden
    for run in study['runs'].values():
        if run['param_value'] == 0.0:
            t_ref = run['tec']['T (s)'].values
            E_ref = run['tec']['E'].values
        else:
            t = run['tec']['T (s)'].values
            E = run['tec']['E'].values
            run['t'] = t
            run['E'] = E

    # 2) Gemeinsames Zeitgitter definieren: log-äquidistant über den relevanten Bereich
    t_grid = np.logspace(np.log10(t_ref.min()), np.log10(t_ref.max()), 2000)

    # 3) Referenz auf das Gitter interpolieren (linear in log-t, das ist für E(t) meist glatter)
    def interp_on_grid(t, E, t_grid):
        f = interp1d(np.log10(t), E, bounds_error=False, fill_value=np.nan)
        return f(np.log10(t_grid))

    E_ref_grid = interp_on_grid(t_ref, E_ref, t_grid)

    results = {}
    for eps, (t, E) in runs.items():
        E_grid = interp_on_grid(t, E, t_grid)

        abs_dev = E_grid - E_ref_grid

        # Relative Abweichung nur dort auswerten, wo E_ref über einer Rauschschwelle liegt
        # (z.B. 0.1 % des Peakwerts), sonst dividiert man durch ~0 im Afterglow-Ende
        threshold = 1e-3 * np.nanmax(E_ref_grid)
        mask = E_ref_grid > threshold
        rel_dev = np.full_like(E_grid, np.nan)
        rel_dev[mask] = abs_dev[mask] / E_ref_grid[mask]

        results[eps] = {
            't': t_grid,
            'abs_dev': abs_dev,
            'rel_dev': rel_dev,
            'max_rel_dev': np.nanmax(np.abs(rel_dev)),
            'nrmse': np.sqrt(np.nanmean((abs_dev[mask])**2)) / np.nanmax(E_ref_grid),
        }

    # 4) Skalare Kennzahl je Epsilon
    for eps, r in results.items():
        print(f"eps={eps:.1f}  max|rel. Abw.|={r['max_rel_dev']*100:6.2f} %   NRMSE={r['nrmse']*100:6.3f} %")

@dataclass
class RunDeviation:
    """Ergebnis der Abweichungsanalyse für einen einzelnen Run."""
    param_value: float   # epsilon
    sigma: float          # Reduktions-Schwellwert (in compare_pkl.py "delta" genannt - dieselbe Größe)
    t: np.ndarray
    values: np.ndarray
    abs_dev: np.ndarray
    rel_dev: np.ndarray
    nrmse_full: float
    nrmse_window: float
    max_rel_dev: float
 
 
class DeviationAnalyzer:
    """
    Analysiert die Abweichung einer Zeitreihe (Standard: Elektronendichte 'E')
    zwischen einem Referenz-Run (voller Mechanismus) und mehreren Vergleichs-
    Runs (z.B. reduzierte Mechanismen bei unterschiedlichen DRGEP-Epsilons).
 
    Parameters
    ----------
    study : dict
        study-dict wie oben beschrieben.
    ref_param_value : float, default 0.0
        param_value, der als Referenz (voller Mechanismus) gilt.
    value_col : str, default "E"
        Name der Spalte in study['runs'][k]['tec'], die verglichen werden soll.
        Kandidaten-Namen werden automatisch gesucht, falls die exakte
        Bezeichnung nicht passt (siehe _find_col).
    time_col : str, default "Time"
        Analog für die Zeitspalte.
    n_grid : int, default 2000
        Anzahl Punkte des log-äquidistanten Auswertungsgitters.
    window : tuple(float, float), default (1e-6, 1e-4)
        Zeitfenster (in s) für die fenstergewichtete NRMSE, z.B. der
        Afterglow-Bereich.
    rel_threshold : float, default 1e-3
        Relative Schwelle (bezogen auf das Maximum der Referenzkurve),
        unterhalb der die relative Abweichung nicht berechnet wird
        (verhindert Division durch ~0 am Ende des Afterglows).
    """
 
    TIME_CANDIDATES = ["T (s)", "time", "Time(s)", "t", "TIME"]
    VALUE_CANDIDATES_DEFAULT = {
        "E": ["E", "E [cm-3]", "ne", "n_e", "electron_density"],
    }
 
    def __init__(
        self,
        study: dict,
        ref_param_value: float = 0.0,
        ref_sigma: float = 0.0,
        value_col: str = "E",
        time_col: str = "T (s)",
        n_grid: int = 2000,
        window: tuple[float, float] = (1e-7, 1e-4),
        rel_threshold: float = 1e-4,
        dic_reduced_path: Optional[Path] = None,
        sigma_outer_key: Optional[str] = None,
    ):
        """
        ref_sigma : float, default 0.0
            Sigma-Wert (= "delta" in compare_pkl.py, Reduktions-Schwellwert)
            des Referenzruns. Zusammen mit ref_param_value (epsilon)
            identifiziert das den Baseline-Run eindeutig.
        dic_reduced_path : Path, optional
            Pfad zur dic_reduced_*.pkl mit dem epsilon/sigma-Mapping
            (siehe Erzeugungsskript / compare_pkl.py). Die Datei heisst
            pro Studie z.B. "dic_reduced_<powplug>.pkl". Wird nicht
            angegeben, sucht DeviationAnalyzer automatisch im selben
            Ordner wie die Run-Unterordner der Study nach genau einer
            dic_reduced_*.pkl.
        sigma_outer_key : str, optional
            Oberster Key in der dic_reduced_*.pkl (z.B. '2e+06'), falls
            dort mehrere Top-Level-Gruppen vorhanden sind. Bei genau
            einem Top-Level-Key wird dieser automatisch verwendet.
        """
        self.study = study
        self.ref_param_value = float(ref_param_value)
        self.ref_sigma = float(ref_sigma)
        self.value_col = value_col
        self.time_col = time_col
        self.n_grid = n_grid
        self.window = window
        self.rel_threshold = rel_threshold
        self.dic_reduced_path = Path(dic_reduced_path) if dic_reduced_path else None
        self.sigma_outer_key = sigma_outer_key
 
        self.t_grid: Optional[np.ndarray] = None
        self.t_ref: Optional[np.ndarray] = None
        self.E_ref: Optional[np.ndarray] = None
        self.E_ref_grid: Optional[np.ndarray] = None
        # run_idx -> {'epsilon': float, 'sigma': float}, befuellt von _load_sigma_mapping()
        self.run_meta: dict = {}
        # Key ist jetzt (epsilon, sigma) statt nur epsilon, damit Runs mit
        # gleichem epsilon aber unterschiedlichem sigma nicht mehr
        # gegenseitig in self.results ueberschrieben werden.
        self.results: dict[tuple[float, float], RunDeviation] = {}
        self.skipped: list[tuple[float, float]] = []
        # Wird in run() gesetzt: True, wenn mehr als ein Sigma-Wert unter
        # den Ergebnissen vorkommt (siehe fmt_key()).
        self.sigma_varies: bool = False
 
    # ------------------------------------------------------------------
    # Spalten-/Datenextraktion
    # ------------------------------------------------------------------
    def _find_col(self, df: pd.DataFrame, requested: str, fallback_key: str) -> str:
        if requested in df.columns:
            return requested
        candidates = self.VALUE_CANDIDATES_DEFAULT.get(fallback_key, []) if fallback_key != "time" \
            else self.TIME_CANDIDATES
        for c in candidates:
            if c in df.columns:
                return c
        hits = [c for c in df.columns if fallback_key.lower() in c.lower()]
        if hits:
            return hits[0]
        raise ValueError(
            f"Spalte '{requested}' nicht gefunden. Verfuegbare Spalten: {list(df.columns)}"
        )
 
    def _extract_tE(self, run_entry: dict) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        tec = run_entry.get("tec")
        if tec is None or not isinstance(tec, pd.DataFrame) or tec.empty:
            return None, None
 
        t_col = self._find_col(tec, self.time_col, "time")
        e_col = self._find_col(tec, self.value_col, self.value_col)
 
        t = tec[t_col].to_numpy(dtype=float)
        E = tec[e_col].to_numpy(dtype=float)
 
        order = np.argsort(t)
        t, E = t[order], E[order]
 
        # Doppelte Zeitpunkte (z.B. Restart-Artefakte) entfernen
        t, unique_idx = np.unique(t, return_index=True)
        E = E[unique_idx]
        return t, E
 
    def _interp_on_grid(self, t: np.ndarray, E: np.ndarray) -> np.ndarray:
        valid = t > 0
        f = interp1d(
            np.log10(t[valid]), E[valid],
            bounds_error=False, fill_value=np.nan,
        )
        return f(np.log10(self.t_grid))
 
    # ------------------------------------------------------------------
    # Sigma-Mapping (epsilon, sigma) pro Run-Index
    # ------------------------------------------------------------------
    def _find_dic_reduced_path(self) -> Path:
        """
        Sucht die dic_reduced_*.pkl automatisch im selben Ordner wie die
        Run-Unterordner der Study (analog zur Ordnerstruktur in
        compare_pkl.py: dic_reduced_*.pkl liegt neben den Run_*-Ordnern,
        aus denen study['runs'] aufgebaut wurde).
        """
        if self.dic_reduced_path is not None:
            return self.dic_reduced_path

        sample_run = next(iter(self.study["runs"].values()))
        pkl_dir = Path(sample_run["path"]).parent
        candidates = sorted(pkl_dir.glob("dic_reduced_*.pkl"))

        if not candidates:
            raise FileNotFoundError(
                f"Keine dic_reduced_*.pkl in '{pkl_dir}' gefunden. "
                f"Bitte dic_reduced_path beim Konstruktor explizit angeben."
            )
        if len(candidates) > 1:
            raise ValueError(
                f"Mehrere dic_reduced_*.pkl in '{pkl_dir}' gefunden: "
                f"{[c.name for c in candidates]}. Bitte dic_reduced_path "
                f"beim Konstruktor explizit angeben."
            )
        return candidates[0]

    def _load_sigma_mapping(self) -> None:
        """
        Baut self.run_meta[run_idx] = {'epsilon': eps, 'sigma': sigma} auf.

        Struktur der dic_reduced_*.pkl (siehe Erzeugungsskript):
            dic_reduced[powplug] = {
                'targets': [...],            # KEIN Run-Eintrag, wird uebersprungen
                epsilon_1: {sigma_a: {...}, sigma_b: {...}, ...},
                epsilon_2 (z.B. 0.0): {0.0: {...}},   # Sonderfall: nur 1 Eintrag
                ...
            }
        'sigma' entspricht 'delta' im Erzeugungsskript/compare_pkl.py -
        dieselbe Groesse, nur andere Bezeichnung. Der oberste Key
        ('powplug') ist pro dic_reduced-Datei einmalig vorhanden.
        Model nummeriert Runs ab 1 (wie in compare_pkl.py); die
        Reihenfolge folgt der Dict-Einfuegereihenfolge im Erzeugungsskript
        (epsilon_list, darin delta_list).
        """
        path = self._find_dic_reduced_path()
        with open(path, "rb") as f:
            dic_reduced_data = pickle.load(f)

        outer_key = self.sigma_outer_key
        if outer_key is None:
            outer_keys = list(dic_reduced_data.keys())
            if len(outer_keys) != 1:
                raise ValueError(
                    f"dic_reduced_data in '{path.name}' hat mehrere "
                    f"oberste Keys {outer_keys}; bitte sigma_outer_key "
                    f"beim Konstruktor angeben."
                )
            outer_key = outer_keys[0]

        run_meta = {}
        run_idx = 1
        for epsilon, eps_data in dic_reduced_data[outer_key].items():
            if epsilon == "targets" or not isinstance(eps_data, dict):
                # 'targets' ist die Liste der Zielspezies, kein Run-Eintrag.
                continue
            for sigma, _data in eps_data.items():
                run_meta[run_idx] = {
                    "epsilon": float(epsilon),
                    "sigma": float(sigma),
                }
                run_idx += 1

        self.run_meta = run_meta

    def _key_for(self, run_idx) -> tuple[float, float]:
        """
        Liefert den (epsilon, sigma)-Key fuer einen Run-Index. Faellt auf
        (param_value, 0.0) zurueck, falls kein Sigma-Mapping fuer diesen
        Run existiert (z.B. Baseline-Run mit epsilon=sigma=0.0, der u.U.
        gar nicht in dic_reduced_data auftaucht).
        """
        meta = self.run_meta.get(run_idx)
        if meta is not None:
            return (meta["epsilon"], meta["sigma"])
        param_value = self.study["runs"][run_idx].get("param_value", float("nan"))
        try:
            param_value = float(param_value)
        except (TypeError, ValueError):
            pass
        return (param_value, 0.0)

    # ------------------------------------------------------------------
    # Hauptanalyse
    # ------------------------------------------------------------------
    def run(self) -> "DeviationAnalyzer":
        self._load_sigma_mapping()
        runs = self.study["runs"]
        ref_target = (self.ref_param_value, self.ref_sigma)
 
        ref_key = next(
            (k for k in runs if self._key_for(k) == ref_target),
            None,
        )
        if ref_key is None:
            found_keys = sorted(self._key_for(k) for k in runs)
            study_param_values = sorted(
                v.get("param_value") for v in runs.values()
            )
            raise ValueError(
                f"Kein Run mit epsilon={self.ref_param_value:.3e}, "
                f"sigma={self.ref_sigma} gefunden.\n"
                f"  Anzahl Runs in study: {len(runs)}\n"
                f"  Anzahl Eintraege aus dic_reduced (run_meta): {len(self.run_meta)}\n"
                f"  (epsilon, sigma)-Paare aus dic_reduced: {found_keys}\n"
                f"  study['runs'][*]['param_value']: {study_param_values}"
            )
 
        self.t_ref, self.E_ref = self._extract_tE(runs[ref_key])
        if self.t_ref is None:
            raise ValueError(
                f"Referenzrun (epsilon={self.ref_param_value}, "
                f"sigma={self.ref_sigma}) hat kein gueltiges tec-DataFrame."
            )
 
        self.t_grid = np.logspace(
            np.log10(self.t_ref[self.t_ref > 0].min()),
            np.log10(self.t_ref.max()),
            self.n_grid,
        )
        self.E_ref_grid = self._interp_on_grid(self.t_ref, self.E_ref)
        threshold = self.rel_threshold * np.nanmax(self.E_ref_grid)
 
        self.results = {}
        self.skipped = []
 
        for key, entry in runs.items():
            eps, sigma = self._key_for(key)
            t, E = self._extract_tE(entry)
            if t is None:
                self.skipped.append((eps, sigma))
                continue
 
            E_grid = self._interp_on_grid(t, E)
            abs_dev = E_grid - self.E_ref_grid
            mask = self.E_ref_grid > threshold
 
            rel_dev = np.full_like(E_grid, np.nan)
            rel_dev[mask] = abs_dev[mask] / self.E_ref_grid[mask]
 
            window_mask = mask & (self.t_grid > self.window[0]) & (self.t_grid < self.window[1])
 
            nrmse_full = (
                np.sqrt(np.nanmean(abs_dev[mask] ** 2))  / np.nanmax(self.E_ref_grid)
                if mask.any() else np.nan
            )
            nrmse_window = (
                np.sqrt(np.nanmean(abs_dev[window_mask] ** 2)) / np.nanmax(self.E_ref_grid)
                if window_mask.any() else np.nan
            )
            max_rel_dev = np.nanmax(np.abs(rel_dev)) if mask.any() else np.nan
 
            self.results[(eps, sigma)] = RunDeviation(
                param_value=eps,
                sigma=sigma,
                t=self.t_grid,
                values=E_grid,
                abs_dev=abs_dev,
                rel_dev=rel_dev,
                nrmse_full=nrmse_full,
                nrmse_window=nrmse_window,
                max_rel_dev=max_rel_dev,
            )
 
        if self.skipped:
            print(f"Warnung: {len(self.skipped)} Runs uebersprungen (kein tec-DataFrame): {self.skipped}")

        # Nur True, wenn tatsaechlich mehr als ein Sigma-Wert unter den
        # Ergebnissen vorkommt - bei reinen Epsilon-Sweeps (sigma immer
        # 0.0) bleibt das False, damit summary()/plot() sigma nicht als
        # unnoetige, konstante Info mitschleppen.
        sigmas = {sigma for (_, sigma) in self.results}
        self.sigma_varies = len(sigmas) > 1

        return self

    def fmt_key(self, eps: float, sigma: float) -> str:
        """
        Formatiert einen (epsilon, sigma)-Key fuer Legenden/Tabellen. Sigma
        wird nur angezeigt, wenn es unter den Ergebnissen tatsaechlich
        variiert (self.sigma_varies) - sonst waere es bei einem reinen
        Epsilon-Sweep nur redundantes "sigma=0.0" in jeder Beschriftung.
        """
        if self.sigma_varies:
            return f"eps={eps}, sigma={sigma}"
        return f"eps={eps}"
 
    # ------------------------------------------------------------------
    # Ausgabe
    # ------------------------------------------------------------------
    def summary(self) -> pd.DataFrame:
        if not self.results:
            raise RuntimeError("Erst .run() aufrufen.")
 
        rows = []
        # self.results ist ueber (epsilon, sigma)-Tupel geschluesselt; die
        # Standard-Tupelsortierung sortiert automatisch erst nach epsilon,
        # dann nach sigma. 'sigma' erscheint als Spalte nur, wenn es
        # tatsaechlich variiert (siehe fmt_key()/sigma_varies).
        for (eps, sigma), r in sorted(self.results.items()):
            row = {self.study.get("param_name", "epsilon"): eps}
            if self.sigma_varies:
                row["sigma"] = sigma
            row.update({
                "NRMSE_gesamt_%": r.nrmse_full * 100,
                f"NRMSE_fenster_{self.window[0]:.0e}-{self.window[1]:.0e}s_%": r.nrmse_window * 100,
                "max_rel_abweichung_%": r.max_rel_dev * 100,
            })
            rows.append(row)
        df = pd.DataFrame(rows)
        print('\n')
        print(df.to_string(index=False, float_format=lambda x: f"{x:8.3f}"))
        print('\n')
        return df
 
    def plot(self, title: Optional[str] = None, figsize=(9, 8)):
        if not self.results:
            raise RuntimeError("Erst .run() aufrufen.")
 
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
 
        ax1.plot(self.t_ref, self.E_ref, "k-", lw=2,
                 label=f"{self.fmt_key(self.ref_param_value, self.ref_sigma)} (voll)")
        for (eps, sigma), r in sorted(self.results.items()):
            ax1.plot(r.t, r.values, label=self.fmt_key(eps, sigma))
        ax1.set_ylabel(f"{self.value_col}")
        ax1.legend(fontsize=8, ncol=2)
        ax1.set_title(title or f"{self.study.get('param_name', '')} – {self.value_col}(t) Abweichung")
 
        for (eps, sigma), r in sorted(self.results.items()):
            ax2.plot(r.t, r.rel_dev * 100, label=self.fmt_key(eps, sigma))
        ax2.axhline(0, color="k", lw=0.5)
        ax2.axvspan(self.window[0], self.window[1], color="gray", alpha=0.1,
                    label="Auswertefenster")
        ax2.set_xscale("log")
        ax2.set_ylabel("rel. Abweichung [%]")
        ax2.set_xlabel("Time (s)")

        plt.tight_layout()
        plt.show()
        return fig
 
    def get(self, epsilon: float, sigma: float = 0.0) -> RunDeviation:
        key = (epsilon, sigma)
        if key not in self.results:
            raise KeyError(
                f"(epsilon={epsilon}, sigma={sigma}) nicht in results. "
                f"Vorhanden: {sorted(self.results)}"
            )
        return self.results[key]

    def get_by_run_index(self, run_idx) -> RunDeviation:
        """
        Bequemer Zugriff ueber den originalen Study-Run-Index (also den Key
        in study['runs']) statt ueber (epsilon, sigma) direkt - nuetzlich
        fuer Code, der Ergebnisse mit anderen Run-indizierten Tabellen
        (z.B. run_table in create_run_report.py) verknuepfen will.
        """
        key = self._key_for(run_idx)
        if key not in self.results:
            raise KeyError(
                f"Run {run_idx} (abgeleiteter Key={key}) nicht in results. "
                f"Vorhanden: {sorted(self.results)}"
            )
        return self.results[key]