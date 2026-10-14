import os
from pathlib import Path
import re
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from functions.functions_analysis import plot_rates, load_study_cache, progress_bar
from functions.gk_rate_mapper import build_reaction_lookup

#BASE_DIR = Path(__file__).parent

def combine_runs(run_path: Path, combined_run_data: list):
    '''
    Combines multiple runs from subdirectories into a single dictionary and saves it as a pickle file.

    :param run_path: Path to the directory containing the subdirectories with the runs
    :param combined_run_data: List of subdirectory names to combine
    '''
    big_df: dict = {
        'param_name': "",
        'date':       combined_run_data[1],
        'time':       combined_run_data[2],
        'runs': {},
    }
    labels_dic: dict = {}
    loaded_runs = {
        'loaded_run': 0,
        'skiped_run': 0
    }

    counter = 0
    print('Lade Runs aus den Unterordnern...')
    for idx, subdir in enumerate((run_path).iterdir()):
        progress_bar(idx + 1, len(list((run_path).iterdir())))
        counter += 1
        #print(subdir)
        if subdir.is_dir():
            if subdir.name == 'Combined_Runs':
                loaded_runs['skiped_run'] += 1
                continue
            #print(f"Lade Run-Ordner: {subdir.name}")
            pkl_path = subdir / f'{subdir.name}.pkl'
            labels_path = subdir / f'{subdir.name}_reaction_labels.pkl'

            
            if pkl_path.exists() and labels_path.exists():
                study = load_study_cache(pkl_path)
                labels = load_study_cache(labels_path) if labels_path.exists() else None
                big_df['param_name'] = 'CombinedRuns'
                loaded_runs['loaded_run'] += 1

                # add study data to big_df
                for n, run in study['runs'].items():    
                    #print(f'n: {n}')
                    if run['tec'] is None or run['nam'] is None:
                        print(f"  ⚠ Run {idx} hat keine Daten (tec oder nam fehlt)")
                        continue
                    #print(f"  Lade Run {n}  (param_value = {run['param_value']})")

                    if loaded_runs['loaded_run'] == 1:
                        big_df['runs'][n] = { 
                            'param_value': run['param_value'],
                            'path':        subdir,
                            'tec': run['tec'],
                            'nam': run['nam'],
                        }
                        #print(f"Hier: {big_df['runs'][n]['path']}")
                    else:
                        if big_df['runs'][n]['param_value'] == run['param_value'] and len(big_df['runs'][n]['tec']['T (s)']) == len(run['tec']['T (s)']):
                            # checking if colums in tec are the same
                            if big_df['runs'][n]['tec'].columns.equals(run['tec'].columns):
                                continue
                            else:
                                # append new columns to big_df['runs'][n]['tec'], only diffrent colums
                                new_cols = run['tec'].columns.difference(big_df['runs'][n]['tec'].columns)
                                #print(f"  Run {n}: Neue Spalten in tec-Datei gefunden: {len(new_cols)}")
                                big_df['runs'][n]['tec'] = pd.concat([big_df['runs'][n]['tec'], run ['tec'][new_cols]], axis=1)

                for key, value in labels.items():
                    if value not in labels_dic.keys():
                        labels_dic[key] = value

    # count how many colums in tec start with rate
    rate_cols = [col for col in big_df['runs'][1]['tec'].columns if col.startswith('RATE')]

    print(f'  {counter} pkl-files geladen\n')
    print(f"  → {len(big_df['runs'])} Runs geladen  |  Parameter: '{big_df['param_name']}'")
    print(f'  → {len(rate_cols)} Spalten in tec-Datei geladen')
    print(f'  → {len(labels_dic)} Reaktionsgleichungen geladen\n')
    #print(labels_dic)

    for idx, run in big_df['runs'].items():
        run_folder = run_path / f'CombinedRuns_{big_df["date"]}_{big_df["time"]}' / f'Run_{idx+1}_{run["param_value"]}'
        os.makedirs(run_folder, exist_ok=True)
        run['path'] = run_folder
        # plot_rates(run, path_labels=run_path / f'{big_df["param_name"]}_{big_df["date"]}_{big_df["time"]}_reaction_labels.pkl')

    #save big_df and labels_dic as pickle file
    with open(run_path / f'CombinedRuns_{big_df["date"]}_{big_df["time"]}' / f'CombinedRuns_{big_df["date"]}_{big_df["time"]}.pkl', 'wb') as f:
        pickle.dump(big_df, f)
    with open(run_path / f'CombinedRuns_{big_df["date"]}_{big_df["time"]}' / f'CombinedRuns_{big_df["date"]}_{big_df["time"]}_reaction_labels.pkl', 'wb') as f:
        pickle.dump(labels_dic, f)