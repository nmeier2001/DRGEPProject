"""
Start GLOBAL_KIN simulation with a single command
"""

import os
from model_change_parameter import changing_parameters, save_outputs, run_model
import time
from datetime import date

# BASE_DIR = Path(__file__).parent

class ModelRunner:
    """
    Start GLOBAL_KIN simulation with a single command

    """
    def __init__(self, BASE_DIR, param_name, param_values, run_folder = 0):
        self.BASE_DIR = BASE_DIR
        self.run_folder = run_folder
        self.param_name = param_name
        self.param_values = param_values

        self.GK_PATH = BASE_DIR / 'Model'
        self.EXE_FILE_NAME = './global_kin.exe'
        self.SH_FILE_NAME  = 'SH_clean.sh'
        self.NAM_FILE_PATH = BASE_DIR / 'Model' / 'global_kin.nam'

        self.run()

    def run(self):
        start_time = time.time()

        # Normalisierung: param_name und param_values immer als Listen behandeln
        if isinstance(self.param_name, str):
            param_names  = [self.param_name]
            param_values = [self.param_values]
        else:
            param_names  = self.param_name
            param_values = self.param_values

        # Liste von Dicts: [{'POWPLUG': 2e7, 'XPOW': 15e-9}, ...]
        param_runs = [
            dict(zip(param_names, values))
            for values in zip(*param_values)
        ]
        n_runs = len(param_runs)

        # Run-Verzeichnis
        run_name = f"{'_'.join(param_names)}_{date.today().strftime('%Y%m%d')}_{time.strftime('%H%M')}"
        print(f"Run directory: {run_name}")
        if self.run_folder:
            run_dir = self.BASE_DIR / 'RUNS' / self.run_folder / run_name
        else:
            run_dir = self.BASE_DIR / 'RUNS' / run_name
        os.makedirs(run_dir, exist_ok=True)
        print(f"Run directory created: {run_dir}")
        with open(run_dir / 'parameters.txt', 'w') as f:
            f.write(f"Parameter: {param_names}\n")
            f.write(f"Werte: {param_runs}\n\n")
            f.write("Rates:\nEdukt: \nProdukt: \n")

        for i, params in enumerate(param_runs):
            print(f"Run {i+1}/{n_runs}: {params}")

            tend = changing_parameters(self.NAM_FILE_PATH, params)

            print()
            print('Starting GLOBAL_KIN simulation')
            returncode, output = run_model(self.GK_PATH, self.EXE_FILE_NAME, tend)

            if returncode != 0:
                print(f"\n  FEHLER: Model beendet mit Fehler – überspringe Speichern.")
            else:
                print(f"\n  Model erfolgreich abgeschlossen.")
                save_outputs(run_dir, self.NAM_FILE_PATH.name, list(params.values()), i)

        end_time = time.time()
        print(f"  Laufzeit: {(end_time - start_time) / 60:.2f} Minuten")
        #os.system('sh ' + self.SH_FILE_NAME)
        #os.system('./' + self.EXE_FILE_NAME)