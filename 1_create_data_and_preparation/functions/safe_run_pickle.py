"""
Model Parameter Study - Multi-Run Loader
=============================================
Lädt alle .tec Dateien aus den Run-Unterordnern einer Parameterstudie.
Basiert auf der GKAnalysis-Logik aus read_tec.py.

06/2025 N. Meier: creation
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
import pickle
from functions.gk_rate_mapper import build_reaction_lookup, get_rate_label

BASE_DIR = Path(__file__).parent
RUN = 'POWPLUG_20260715_1849'

# Spezies die beim Einlesen ignoriert werden
EXCLUDE_GROUP = ['MIS']

# Alle CO2-Spezies (für Konversionsberechnung)
CO2_GROUP = [
    'CO2', 'CO2va', 'CO2vb', 'CO2vc', 'CO2vd',
    *[f'CO2v{i:02d}' for i in range(1, 22)],
    'CO2e1', 'CO2e2'
]


# ── .tec Parser (aus GKAnalysis übernommen) ────────────────────────────────────

def parse_tec_file(filepath: Path) -> pd.DataFrame:

    tec_info_part = []
    tec_data_part = []
    header_done = False

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '':
                header_done = True
            elif not header_done:
                #print(line)
                #line = line.replace('\x00', '')
                tec_info_part.append(line)
            else:
                tec_data_part.append(line)

    # ── Header: exakt wie in GKAnalysis ───────────────────────────────────────
    variables = []
    step_no = None
    for i, line in enumerate(tec_info_part):
        #print(f"  Header-Zeile {i}: {line}")
        if i == 0:
            pass  # Titel-Zeile
        elif i < len(tec_info_part) - 1:
            if '"' in line:
                var_name = line.split('"')[1].strip()
                if var_name:
                    variables.append(var_name)
        else:
            # letzte Header-Zeile: I = XXXX
            m = re.findall(r'\d+', line)
            if m:
                step_no = int(m[0])

    if not variables or step_no is None:
        raise ValueError(f"Header konnte nicht gelesen werden: {filepath}")

    # ── Daten: exakt wie in GKAnalysis ────────────────────────────────────────
    tec_data_full = []
    for line in tec_data_part:
        tec_data_full.extend(line.split(' '))
    tec_data_full = ' '.join(tec_data_full).split()  # leere Einträge entfernen

    # Fortran-Notation → Python float (beide Fälle)
    tec_data_full = [re.sub(r'([0-9])([+\-])(\d+)$', r'\1E\2\3', s) for s in tec_data_full]
    tec_data_full = [re.sub(r'[Dd]([+\-]?\d+)', r'E\1', s) for s in tec_data_full]

    # ── Fix für unterschiedliche Run-Längen ───────────────────────────────────
    n_vars = len(variables)
    actual_step_no = len(tec_data_full) // n_vars

    # ── Variablen spaltenweise aufteilen ─────────────────────────────────────
    var_data = {}
    for i, var in enumerate(variables):
        if var not in EXCLUDE_GROUP:
            var_data[var] = tec_data_full[i * actual_step_no : (i + 1) * actual_step_no]

    df = pd.DataFrame(var_data)
    df = df.astype('float64')
    df = _add_custom_parameters(df)

    return df


def _add_custom_parameters(df: pd.DataFrame) -> pd.DataFrame:
    new_cols = {}

    # CO2-Konversion
    present_co2 = [s for s in CO2_GROUP if s in df.columns]
    if present_co2:
        n_co2_zero = df[present_co2].iloc[0].sum()
        co2_total  = df[present_co2].sum(axis=1)
        new_cols['CONVERSION'] = (n_co2_zero - co2_total) / n_co2_zero

    # Vibrationstemperatur Tv1
    if 'CO2' in df.columns and 'CO2v01' in df.columns:
        new_cols['Tv1'] = (-4.08 - -3.78851) / (
            8.817333e-5 * np.log(df['CO2v01'] / df['CO2'])
        )

    # Alle neuen Spalten auf einmal hinzufügen → kein Fragmentierungs-Warning
    if new_cols:
        df = pd.concat([df, pd.DataFrame(new_cols, index=df.index)], axis=1)

    return df


# ── .nam Parser ────────────────────────────────────────────────────────────────

def parse_nam_file(filepath: Path) -> dict:
    """
    Liest Parameter aus einer model.nam Datei.
    Gibt Dict {Parametername: float} zurück.
    """
    params = {}
    
    # Regex: KEY = Wert(e), optional gefolgt von Komma oder !-Kommentar
    pattern = re.compile(
        r'^\s*([A-Z_][A-Z0-9_]*)\s*=\s*([^!$\n]+)',
        re.IGNORECASE
    )
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            # Kommentare und Namelist-Tags überspringen
            stripped = line.strip()
            if stripped.startswith('!') or stripped.startswith('$') or stripped.startswith('*'):
                continue
            
            match = pattern.match(line)
            if not match:
                continue
            
            key = match.group(1).upper()
            raw = match.group(2).strip().rstrip(',').strip()
            
            # Fortran D-Notation → Python E-Notation
            raw = re.sub(r'([0-9]\.?)D([+-]?[0-9])', r'\1E\2', raw, flags=re.IGNORECASE)
            
            # Mehrere Werte (Array)?
            values = [v.strip() for v in raw.split(',') if v.strip()]
            parsed = []
            for v in values:
                try:
                    parsed.append(float(v))
                except ValueError:
                    parsed.append(v)
            
            params[key] = parsed[0] if len(parsed) == 1 else parsed
    return params


# ── Ordnernamen-Parser ─────────────────────────────────────────────────────────

def parse_study_dir(dir_name: str) -> dict:
    """
    POWPLUG_20260626_1038        →  {param_name: 'POWPLUG', date, time}
    POWPLUG_XPOW_20260626_1038   →  {param_name: ['POWPLUG', 'XPOW'], date, time}
    """
    m = re.match(r'(.+?)_(\d{8})_(\d{4})$', dir_name)
    if not m:
        raise ValueError(f"Unerwartetes Studien-Ordner-Format: {dir_name}")

    param_part = m.group(1)
    params = param_part.split('_')

    # Einzelner Parameter oder mehrere
    param_name = params[0] if len(params) == 1 else params

    return {'param_name': param_name, 'date': m.group(2), 'time': m.group(3)}


def parse_run_dir(dir_name: str) -> dict | None:
    """
    Run_1_2000.0              →  {run_index: 1, param_value: 2000.0}
    Run_1_[20000000.0, 1.5e-08] →  {run_index: 1, param_value: [20000000.0, 1.5e-08]}
    Gibt None zurück wenn Format nicht passt.
    """
    # Einzelner Wert: Run_1_2000.0
    m = re.match(r'Run_(\d+)_([\d.eE+\-]+)$', dir_name)
    if m:
        return {'run_index': int(m.group(1)), 'param_value': float(m.group(2))}

    # Liste: Run_1_[20000000.0, 1.5e-08]
    m = re.match(r'Run_(\d+)_\[(.+)\]$', dir_name)
    if m:
        values = [float(v.strip()) for v in m.group(2).split(',')]
        return {'run_index': int(m.group(1)), 'param_value': values}

    return None


# ── Haupt-Ladefunktion ─────────────────────────────────────────────────────────

def load_study(study_dir: Path) -> dict:
    """
    Lädt alle Runs einer Parameterstudie.

    Rückgabe-Struktur:
    {
        'param_name': 'POWPLUG',
        'date':       '20260626',
        'time':       '1038',
        'runs': {
            1: {
                'param_value': 2000.0,
                'path':        Path(...),
                'tec':         pd.DataFrame,   # None wenn Datei fehlt
                'nam':         dict,            # None wenn Datei fehlt
            },
            ...
        }
    }
    """
    study = {**parse_study_dir(study_dir.name), 'runs': {}}

    run_dirs = sorted(
        (d for d in study_dir.iterdir() if d.is_dir()),
        key=lambda d: parse_run_dir(d.name)['run_index']
                      if parse_run_dir(d.name) else -1
    )

    for subdir in run_dirs:
        info = parse_run_dir(subdir.name)
        if info is None:
            print(f"  Überspringe '{subdir.name}' – kein Run-Ordner")
            continue

        idx = info['run_index']
        param_value = info['param_value']

        # Formatierung: Liste oder einzelner Wert
        if isinstance(param_value, list):
            param_value_fmt = [float(f"{v:.3e}") for v in param_value]
        else:
            param_value_fmt = f"{param_value:.3e}"

        print(f"  Lade Run {idx}  (param_value = {param_value_fmt})")

        tec_path = subdir / 'global_kin.tec'
        nam_path = subdir / 'global_kin.nam'

        study['runs'][idx] = {
            'param_value': param_value_fmt,
            'path':        subdir,
            'tec':         parse_tec_file(tec_path) if tec_path.exists() else None,
            'nam':         parse_nam_file(nam_path) if nam_path.exists() else None,
        }

        if not tec_path.exists():
            print(f"    ⚠ Keine global_kin.tec in {subdir.name}")

    print(f"\n→ {len(study['runs'])} Runs geladen  |  Parameter: '{study['param_name']}'")
    return study


# def study_to_dataframe(study: dict) -> pd.DataFrame:
#     """
#     Fasst alle Runs in einem DataFrame zusammen.
#     Jede Zeile = ein Zeitschritt, zusätzliche Spalten: 'run_index', 'param_value'.
#     Nützlich für schnelle Vergleichsplots über alle Runs.
#     """
#     frames = []
#     for idx, run in study['runs'].items():
#         if run['tec'] is None:
#             continue
#         df = run['tec'].copy()
#         df.insert(0, 'param_value', run['param_value'])
#         df.insert(0, 'run_index',   idx)
#         frames.append(df)

#     if not frames:
#         return pd.DataFrame()

#     combined = pd.concat(frames, ignore_index=True)
#     return combined

def save_study(study: dict, path: Path):
    with open(path, 'wb') as f:
        pickle.dump(study, f)
    print(f"Study gespeichert: {path}")


def load_save_study(dir_path: Path, dat_path):
    study = load_study(dir_path)

    rxn_map = build_reaction_lookup(dat_path)
    df = study['runs'][1]['tec']

    rate_cols = [c for c in df.columns if c.startswith("RATE_")]
    labels = {col: get_rate_label(rxn_map, int(col.split("_")[-1]), include_comment=True)
              for col in rate_cols}

    # Dateiname: ein Parameter oder mehrere
    param_name = study["param_name"]
    if isinstance(param_name, list):
        param_str = '_'.join(param_name)
    else:
        param_str = param_name

    base_name = f'{param_str}_{study["date"]}_{study["time"]}'

    with open(dir_path / f'{base_name}_reaction_labels.pkl', 'wb') as f:
        pickle.dump(labels, f)
    print(f"Reaktionslabels gespeichert: {dir_path / f'{base_name}_reaction_labels.pkl'}")

    save_study(study, dir_path / f'{base_name}.pkl')

# ── Einstiegspunkt ─────────────────────────────────────────────────────────────

if __name__ == '__main__':
    study = load_study(BASE_DIR / 'RUNS' / RUN)

    rxn_map = build_reaction_lookup(BASE_DIR / 'Model' / "global_kin.dat")

    df = study['runs'][1]['tec']

    # Beim Plotten: Spaltennamen automatisch in Reaktionsgleichungen umwandeln
    rate_cols = [c for c in df.columns if c.startswith("RATE_")]
    labels = {col: get_rate_label(rxn_map, int(col.split("_")[-1]), include_comment=True)
            for col in rate_cols}

    with open(BASE_DIR / 'RUNS' / RUN / f'{study["param_name"]}_{study["date"]}_{study["time"]}_reaction_labels.pkl', 'wb') as f:
        pickle.dump(labels, f)

    # safe study to pickle
    save_study(study, BASE_DIR / 'RUNS' / RUN / f'{study["param_name"]}_{study["date"]}_{study["time"]}.pkl')

