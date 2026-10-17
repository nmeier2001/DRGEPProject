from pathlib import Path

from functions.functions_analysis import add_EN_to_study, calculate_SEI, calculate_energy_efficiency, calculate_ratio, plot_xpow_size, plot_rates, load_study_cache, check_file_exists, DeviationAnalyzer
from functions.plot_overview import PlotOverview
from functions.safe_run_pickle import load_save_study
from functions.combine_pkl import combine_runs
from functions.safe_run_pickle import save_study

BASE_DIR = Path(__file__).parent

rate_study_path = BASE_DIR / "Runs" / 'Rate study 1 Puls 10000 Datapoints'
# runs_rate_study = rate_study_path.iterdir()
# print(f"Rate Study runs found in {rate_study_path}: {[run.name for run in runs_rate_study]}")

# runs = [run.name for run in runs_rate_study]
runs = ['Reduced_20260921_2328']
plot_rates_flag = True

for RUN in runs:
    if RUN.startswith('CombinedRuns'):
        combined_run = True
        run_path = rate_study_path
    elif len(runs) >= 2:
        print(f"⚠ {RUN} ist ein Rate Study Run.")
        combined_run = False
        run_path = rate_study_path
    else:
        print(f"⚠ {RUN} ist kein CombinedRun.")
        combined_run = False
        run_path = BASE_DIR / 'RUNS'

    study_pkl_path = run_path / RUN / f'{RUN}.pkl'

    #-----------------------------Study/Reaktionslabels ggf. erst erzeugen-----------------------------
    if not combined_run:
        if not check_file_exists(str(study_pkl_path)): # check if study pickle file exists
            # try:
            load_save_study(run_path / RUN, BASE_DIR / 'Model' / "global_kin.dat")
            # except Exception as e:
            #     print(f"Fehler beim Erzeugen der Study-Pickle-Datei: {e}")
            #     sys.exit(1)
        else:
            print(f"Study-Pickle-Datei existiert bereits: {study_pkl_path}")
    else:
        combined_run_data = RUN.split('_')
        combine_runs(run_path, combined_run_data)
    #-----------------------------Study laden-----------------------------
    study = load_study_cache(study_pkl_path)

    if RUN.startswith('Reduced'):
        da = DeviationAnalyzer(study, ref_param_value='0.000e+00', value_col="E")
        da.run()
        da.summary()
        da.plot()

    param_name = study["param_name"]
    if isinstance(param_name, list):
        param_str = '_'.join(param_name)
    else:
        param_str = param_name

    if column_exists := 'EN_Td' not in study['runs'][1]['tec'].columns:
        add_EN_to_study(study)

    SEI_list = []
    energy_efficiency_list = []
    e_density_list = []
    conversion_list = []
    time = study['runs'][1]['tec']['T (s)']

    #-----------------------------Analyse Runs-----------------------------
    for i, run in study['runs'].items():
        if run['tec'] is None or run['nam'] is None:
            print(f"  ⚠ Run {i} hat keine Daten (tec oder nam fehlt)")
            continue

        SEI_list.append(calculate_SEI(run['nam']))
        n_CO2_0 = run['tec']['CO2'][0]
        n_CO2_final = min(run['tec']['CO2'])

        #print(f'Run {i}: {run["tec"]["TGAS"].iloc[-1]}')

        print(f"Run {i}: n_CO2_0 = {n_CO2_0}, n_CO2_final = {n_CO2_final}")

        conversion_list.append(1 - ((n_CO2_final) / n_CO2_0))
        energy_efficiency = calculate_energy_efficiency(n_CO2_0, n_CO2_final, SEI_list[-1])
        energy_efficiency_list.append(energy_efficiency)

        e_density = run['tec']['E']
        e_density_list.append(e_density)

        plot_rates(run, path_labels=run_path / RUN / f'{param_str}_{study["date"]}_{study["time"]}_reaction_labels.pkl', is_rate_study=combined_run, plot_rates_flag=plot_rates_flag)

    species = ['O', 'O2', 'CO', 'C']
    calculate_ratio(study, species)

    #-----------------------------XPOW-----------------------------
    if study['param_name'] == 'XPOW':
        xpow_values = []
        for i, run in study['runs'].items():
            xpow_values.append(float(run['param_value']))

        plot_xpow_size(xpow_values)

    #plot_xpow_from_study(study, path_run=run_path / RUN)

    PlotOverview(study, path_save=run_path / RUN, conversion_list=conversion_list, SEI_list=SEI_list, energy_efficiency_list=energy_efficiency_list, species=species, is_rate_study=combined_run)

    # safe study to pickle
    save_study(study, run_path / RUN / f'{param_str}_{study["date"]}_{study["time"]}.pkl')

    # print('')
    # print('=== Analyseergebnisse ===')
    # print(f"SEI (eV): {SEI_list}")
    # print(f"Energieeffizienz: {energy_efficiency_list}")