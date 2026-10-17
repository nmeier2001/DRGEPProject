'''
Hier finden sich die Funktionen, die die Parameter in den nam-Datei ändern.
Ebenso das Starten der Simulation und das Speichern der Output-Dateien in den jeweiligen Unterordnern.
'''

import math
import re
from pathlib import Path
import subprocess
import sys
import os
import shutil

BASE_DIR = Path(__file__).parent

OUTPUT_BASE_DIR = BASE_DIR / 'RUNS'
FILE_PATH = BASE_DIR / 'Model'

# Output files which will be saved after each run of GLOBAL_KIN simulation
OUTPUT_FILES  = [
    "global_kin.tec",
    "global_kin.out",
    "global_kin.pdt",
    "global_kin_boltz.out",
    "global_kin.rates",  
    "global_kin.dat"
]

def format_fortran(value):
    """Formatiert einen Wert in Fortran D-Notation, z.B. 2.0D+08"""
    exp = int(f"{value:0e}".split('e')[1])
    mantissa = value / (10 ** exp)
    return f"{mantissa:.2f}D{exp:+d}"


def changing_parameters(nam_path, params: dict):
    """
    Change one or multiple parameters in input file of GLOBAL_KIN simulation

    :param nam_path: path to input file
    :param params: dict of {param_name: param_value}
                   e.g. {'POWPLUG': 2e7, 'XPOW': 15e-9}
    """

    with open(nam_path, 'r') as file:
        nam_file = file.read()

    if not Path(nam_path).exists():
        raise FileNotFoundError(f"Input file 'global_kin.nam' not found.")

    VALID_PARAMS = ['POWPLUG', 'TGAS0', 'TCYCLE', 'TEND', 'XPOW', 'CRATE_PLOT', 'ISIDE_PLOT', 'PRESS0']
    #print(params)

    for param_name, param_value in params.items():
        if param_name not in VALID_PARAMS:
            raise ValueError(f"Parameter '{param_name}' is not valid.")

        if param_name == 'POWPLUG':
            fortran_val = format_fortran(param_value)
            # Matcht ersten Wert (0.0), optional zweiten (0.0), dann den Peak
            nam_file = re.sub(
                r'(POWPLUG\b\s*=\s*0\.0\s*,(?:\s*0\.0\s*,)?\s*)([^\s,]+)(\s*,)',
                lambda m: m.group(1) + fortran_val + m.group(3),
                nam_file
            )
        elif param_name == 'XPOW':
            peak = param_value
            end  = 2 * param_value
            # Matcht ersten Wert (0. oder beliebig), optional zweiten (Startpunkt), dann Peak und Ende
            nam_file = re.sub(
                r'(XPOW\s*=\s*[^\s,]+\s*,(?:\s*[^\s,]+\s*,)?\s*)([^\s,]+)(\s*,\s*)([^\s,]+)(\s*,)',
                lambda m: m.group(1) + format_fortran(peak) + m.group(3) + format_fortran(end) + m.group(5),
                nam_file
            )
        elif param_name == 'CRATE_PLOT':
            new_values = ', '.join(str(v) for v in param_value) + ','
            nam_file = re.sub(
                r'(CRATE_PLOT\s*=\s*)([^\n!]+)',
                lambda m: m.group(1) + new_values,
                nam_file
            )
        elif param_name == 'ISIDE_PLOT':
            new_values = ', '.join(str(v) for v in param_value) + ','
            nam_file = re.sub(
                r'(ISIDE_PLOT\s*=\s*)([^\n!]+)',
                lambda m: m.group(1) + new_values,
                nam_file
            )
        else:
            nam_file = re.sub(
                rf'({re.escape(param_name)}\s*=\s*)([^\s,!]+)',
                lambda m: m.group(1) + str(param_value),
                nam_file
            )

        print(f"  {param_name} gesetzt: {param_value}")

    tend_match = re.search(r'TEND\s*=\s*([\d.E+\-D]+)', nam_file)
    tend_value = float(tend_match.group(1).replace('D', 'E').replace('d', 'e')) if tend_match else None

    with open(nam_path, 'w') as f:
        f.write(nam_file)

    return tend_value


def save_outputs(run_dir, NAM_FILE, param_val, run_index):
    """
    Copy output files from the GLOBAL_KIN simulation to a specified directory.

    :param run_dir: directory where the output files will be saved
    :param NAM_FILE: name of the .nam file used for the simulation
    :param param_val: value of the parameter being tested
    :param run_index: index of the current simulation run
    """
    run_name = f"Run_{run_index+1}_{param_val}"
    run_dir = os.path.join(run_dir, run_name)
    os.makedirs(run_dir, exist_ok=True)

    for fname in OUTPUT_FILES:
        src = os.path.join(FILE_PATH, fname)
        dst = os.path.join(run_dir, fname)
        if os.path.exists(src):
            shutil.copyfile(src, dst)
            print(f"  Gespeichert: {fname}")
        else:
            print(f"  WARNUNG: {fname} nicht gefunden, übersprungen.")

    shutil.copyfile(
        os.path.join(FILE_PATH, NAM_FILE),
        os.path.join(run_dir, NAM_FILE)
    )
    print(f"  Gespeichert: {NAM_FILE}")

def run_model(GK_DIR, GK_EXECUTABLE, tend, dtinit=None):
    """
    Starts the GLOBAL_KIN simulation and tracks its progress based on the TEND value.

    :param GK_DIR: directory where the GLOBAL_KIN executable is located
    :param GK_EXECUTABLE: name of the GLOBAL_KIN executable
    :param tend: TEND value for progress tracking
    :param dtinit: initial time step for logarithmic progress calculations
    :return: return code of the process and the full output as a string
    """

    clean = subprocess.run(
        "sh SH_clean.sh",
        cwd=GK_DIR,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    if clean.returncode != 0:
        print(f"  WARNUNG: SH_clean.sh fehlgeschlagen (Code {clean.returncode})")
        print(f"  STDERR: {clean.stderr[-200:]}")

    process = subprocess.Popen(
        GK_EXECUTABLE,
        cwd=GK_DIR,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,   # stderr in denselben Stream mergen
        text=True
    )

    pattern = re.compile(r'TIME=\s*([\d.E+\-]+)\s*s.*\(\s*(\d+)\)')
    log_start = math.log10(dtinit) if (dtinit is not None and dtinit > 0) else None
    output_lines = []

    for line in process.stdout:
        output_lines.append(line)

        match = pattern.search(line)
        if match:
            t_current = float(match.group(1))
            step      = int(match.group(2))

            # 1. Schutz gegen log10(0.0) beim Simulationsstart
            if t_current <= 0.0:
                continue

            # 2. Falls dtinit nicht übergeben wurde: Ersten positiven Zeitschritt erfassen
            if log_start is None:
                log_start = math.log10(t_current)

            # 3. Falls t_current kleiner als 10**log_start (dtinit) ist, t_current abfangen
            current_log_val = math.log10(max(t_current, 10**log_start))
            log_end = math.log10(tend)

            # 4. Fortschritt berechnen mit Abfängen von Fehlerfällen
            if log_end > log_start:
                progress = max(0.0, min(
                    (current_log_val - log_start) / (log_end - log_start) * 100,
                    100.0
                ))
            else:
                progress = 100.0

            bar_len = 40
            filled  = int(bar_len * progress / 100)
            bar     = '█' * filled + '░' * (bar_len - filled)

            sys.stdout.write(
                f'\r[{bar}] {progress:5.1f}%  '
                f't={t_current:.3e} s  '
                f'Schritt {step}'
            )
            sys.stdout.flush()

    process.wait()
    
    # Am Ende einmal sauber 100% erzwingen für die Ausgabe im Terminal
    sys.stdout.write(f'\r[{"█" * 40}] 100.0%  Simulation abgeschlossen.    \n')
    sys.stdout.flush()

    full_output = "".join(output_lines)
    returncode = 1 if "ERROR" in full_output[:300] else 0
    
    return returncode, full_output