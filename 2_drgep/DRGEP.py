from functions.create_matrix import create_matrix, load_rate_vector, select_species, filter_reactions, filter_reactions_with_sigma, compute_RTS, plot_RAB_fraction_vs_epsilon
from functions.select_timepoints_one_puls import select_timepoints_plateau
import pickle as pkl
from functions.change_dat import ChangeDat
import matplotlib.pyplot as plt
from functions.compare_reduced_dat import create_data_reduced_dic
from pathlib import Path
from functions.setup_reduced_run_folders import prepare_reduced_run_folders

# Set up font size for plots
plt.rcParams.update({'font.size': 14})

BASE_DIR = Path(__file__).parent

# Load the reaction labels and study data, input for DRGEP
labels_path = BASE_DIR / "CombinedRuns_20260715_0828_reaction_labels.pkl"
study_path = BASE_DIR / "CombinedRuns_20260715_0828.pkl"

# DRGEP parameters
epsilon_list = [0.01]  # Threshold for DRGEP, Spezies reduction; for analysis epsilon has to start with 0.0
sigma_list = [0.00] # Threshold for DRGEP, Reaction reduction
targets = ['CO2', 'CO'] # List of species to be retained in the reduced mechanism

# Index of the run to be processed
run_idx = 5 

# Set to True to plot the fraction of kept species vs epsilon; recalculates R_T_S for the given targets
plot_kept_species_vs_epsilon = True 

# path where the reduced dat files are stored, used for preparing the run folders
DAT_SOURCE = BASE_DIR / "dat-file" / "modified_dat"
# target path where the run folders for the reduced mechanisms will be created
TARGET_RUN_DIR = Path(r"C:\Users\NiklasMeier\Desktop\1_create_data_and_preparation\Runs\Reduced_20260921_2328")

with open(study_path, 'rb') as f:
        study = pkl.load(f)

nu, rxn_info, species_list, rxn_list = create_matrix(labels_path)
print(nu.shape, len(species_list), len(rxn_list), '\n')  # Ausgabe der Dimensionen der Stöchiometriematrix und der Listen

# find the timepoints based on E/N
timepoints = select_timepoints_plateau(study['runs'][run_idx]['tec']['T (s)'], study['runs'][run_idx]['tec']['EN_Td'])
powplug = study['runs'][run_idx]['param_value']


dic_rates = {}
for idx, time in timepoints.items():
    rate_vector, rate_cols = load_rate_vector(study_path, time['time'], run_idx)
    dic_rates[idx] = (rate_vector, rate_cols)

# dic for saving parameters and results
dic_reduced = {}
dic_reduced[powplug] = {}
dic_reduced[powplug]['targets'] = targets

# compute R_T_S for the given targets
R_T_S = compute_RTS(nu, dic_rates, species_list, targets)

for epsilon in epsilon_list:
    dic_reduced[powplug][epsilon] = {}

    # species selection based on R_T_S and epsilon
    kept_species = select_species(
        R_T_S, species_list,
        targets=targets,
        epsilon=epsilon
    )

    # needed for the edge case when epsilon = 0.0, where no reactions are eliminated
    if epsilon == 0.0:
        skeletal_rxn_info = filter_reactions(rxn_info, kept_species)

        dic_reduced[powplug][epsilon][0.0] = create_data_reduced_dic(
            kept_species=kept_species,
            eliminated_species=[''],
            number_reactions=len(skeletal_rxn_info),
            powplug=powplug,
            epsilon=epsilon
        )

        skeletal_mechanism_path = BASE_DIR / f"skeletal_mechanism_{epsilon}_0.0.pkl"
        with open(skeletal_mechanism_path, 'wb') as f:
            pkl.dump(skeletal_rxn_info, f)

        ChangeDat(
            dat_file_path=BASE_DIR / "Model" / "model_complete.dat",
            skeletal_mechanism_path=skeletal_mechanism_path,
            epsilon=epsilon,
            sigma=0.0
        )

    else:
        # Sigma-Loop only for epsilon > 0.0
        for sigma in sigma_list:
            dic_reduced[powplug][epsilon][sigma] = {}

            skeletal_rxn_info = filter_reactions_with_sigma(
                nu, rxn_info, kept_species, species_list,
                dic_rates, sigma=sigma
            )

            dic_reduced[powplug][epsilon][sigma] = create_data_reduced_dic(
                kept_species=kept_species,
                eliminated_species=[s for s in species_list if s not in kept_species],
                number_reactions=len(skeletal_rxn_info),
                powplug=powplug,
                epsilon=epsilon
            )

            skeletal_mechanism_path = BASE_DIR / f"skeletal_mechanism_{epsilon}_{sigma}.pkl"
            with open(skeletal_mechanism_path, 'wb') as f:
                pkl.dump(skeletal_rxn_info, f)

            ChangeDat(
                dat_file_path=BASE_DIR / "dat-file" / "model_complete.dat",
                skeletal_mechanism_path=skeletal_mechanism_path,
                epsilon=epsilon,
                sigma=sigma
            )

prepare_reduced_run_folders(
        source_dat_dir=DAT_SOURCE,
        target_base_dir=TARGET_RUN_DIR
    )

if plot_kept_species_vs_epsilon:
    plot_RAB_fraction_vs_epsilon(
        study_path=study_path,
        labels_path=labels_path,
        run_indices=[1, 2, 3, 4, 5],
        timepoint_fn=select_timepoints_plateau,
        targets = targets,
        epsilon_list=[0.1, 0.5, 0.98],
        save_dir=BASE_DIR / "plots" / "RAB_fractions"
    )
