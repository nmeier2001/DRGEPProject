# 2_drgep – DRGEP mechanism reduction

Step 3 of the workflow (see main README). This folder reduces the complete
Model mechanism with the Directed Relation Graph with Error Propagation
(DRGEP) method and writes reduced `.dat` files that can be simulated in
step 4.

## Contents

| File / folder | Purpose |
|---------------|---------|
| `DRGEP.py` | main script: species and reaction reduction, writes reduced `.dat` files |
| `functions/` | helper modules (see below) |
| `dat-file/Model_complete.dat` | complete mechanism (input) |
| `dat-file/modified_dat/` | reduced mechanisms (output) |
| `plots/` | plots, e.g. `plots/RAB_fractions/` |
| `CombinedRuns_<date>_<time>.pkl` | combined study from step 2 (input) |
| `CombinedRuns_<date>_<time>_reaction_labels.pkl` | reaction labels of the rate columns (input) |

## Input
- The combined study pickle and its `_reaction_labels.pkl` from step 2,
  copied into this folder.
- The complete mechanism as `dat-file/Model_complete.dat`. It must be
  the same `.dat` file that was used for the rate simulation in step 1,
  otherwise the rate indices do not match the reactions.

## Configuration
All settings are at the top of `DRGEP.py`:

| Variable | Meaning |
|----------|---------|
| `labels_path`, `study_path` | input pickles from step 2 |
| `epsilon_list` | species threshold ε. Species whose interaction coefficient with all targets is below ε are removed. `0.0` = no species removed (reference) |
| `delta_list` | reaction threshold δ, only used for ε > 0. [TODO: exact criterion of `filter_reactions_with_delta`] |
| `targets` | target species that are always kept, e.g. `['CO2', 'CO']` |
| `run_idx` | index of the run in the study (1-based) whose rates are used, i.e. the POWPLUG value the reduction is based on |
| `plot_kept_species_vs_epsilon` | `True`: additionally plot the fraction of kept species vs. ε |

Orientation from the thesis (CO₂ mechanism, POWPLUG = 2e6): up to ε ≈ 0.50
the error stays small (50 species, about 2824 reactions, NRMSE ≈ 0.45 %);
above that the error increases sharply. δ = 0.1 reduces the number of
reactions strongly (2824 → 173) at the same species count, but with a
higher error.

## Run
The reduction and the creation of the run folders are two separate steps:

```bash
# 1. reduce the mechanism -> dat-file/modified_dat/*.dat
python DRGEP.py

# 2. create the run folders for Model from all reduced .dat files
python functions/setup_reduced_run_folders.py
```

Run step 2 only after all desired ε/δ combinations have been reduced, and
make sure `dat-file/modified_dat/` contains only the `.dat` files that
should be simulated, because every `.dat` file in this folder gets its own
run folder.

### Creating the run folders (`setup_reduced_run_folders.py`)
Settings at the bottom of the script:

| Variable | Meaning |
|----------|---------|
| `SOURCE_DAT_DIR` | folder with the reduced `.dat` files, default `2_drgep/dat-file/modified_dat` |
| `TARGET_RUN_DIR` | folder for the run folders, default `3_model/RUNS/Reduced_<date>_<time>` |

The script
1. collects all `.dat` files in `SOURCE_DAT_DIR`,
2. puts the reference mechanism (ε = 0, δ = 0, file name ending in
   `_0.0_0.0.dat`) first,
3. sorts all other files **numerically** by the numbers in the file name
   (ε, then δ),
4. creates one folder `Run_<n>_<suffix>` per file in `TARGET_RUN_DIR` and
   copies the mechanism into it as `model.dat`.

Example:
```
3_model/RUNS/Reduced_20261005_1219/
├── Run_1_0.0_0.0/model.dat     reference (complete mechanism)
├── Run_2_0.1_0.0/model.dat
├── Run_3_0.5_0.0/model.dat
└── Run_4_0.5_0.1/model.dat
```
The folder `Reduced_<date>_<time>` is then simulated in step 4
(`3_model/model_reduced_run.py`, `REDUCED_FOLDER`).

## How the code works

1. **Stoichiometric matrix** – `create_matrix(labels_path)` builds the matrix
   `nu` (species × reactions) together with `rxn_info`, `species_list`
   and `rxn_list` from the reaction labels.
2. **Time points** – `select_timepoints_plateau()` selects the time points
   used for the reduction from the reduced electric field E/N (`EN_Td`)
   of the chosen run. [TODO: criterion, number of points]
3. **Rate vectors** – `load_rate_vector()` reads the reaction rates of all
   reactions at each selected time point.
4. **Interaction coefficients** – `compute_RTS()` computes the DRGEP
   coefficients R_T,S between each target T and every species S.
   Standard DRGEP (Pepiot-Desjardins & Pitsch, 2008): the direct
   interaction coefficient r_AB measures how strongly the production and
   consumption of A depend on reactions involving B; R_TS is the maximum
   product of r along all paths from T to S. [TODO: check against
   implementation, e.g. maximum over time points and targets]
5. **Species selection** – `select_species()` keeps all species with
   R_TS ≥ ε plus the targets.
6. **Reaction selection**
   - ε = 0: `filter_reactions()` keeps all reactions that contain only
     kept species.
   - ε > 0: `filter_reactions_with_delta()` additionally removes reactions
     below the threshold δ.
7. **Writing the mechanism** – the kept reactions are saved as
   `skeletal_mechanism_<ε>_<δ>.pkl`, and `ChangeDat` writes the reduced
   `.dat` file into `dat-file/modified_dat/`. [TODO: file naming]
8. **Run folders** – done separately with
   `functions/setup_reduced_run_folders.py`, see "Run" above.
9. **Optional plot** – `plot_RAB_fraction_vs_epsilon()` plots the fraction
   of kept species vs. ε for several runs. Run indices and ε values for this
   plot are set directly in the function call at the end of the script,
   independent of `epsilon_list`.

## Output
- `dat-file/modified_dat/*.dat` – reduced mechanisms
- `skeletal_mechanism_<ε>_<δ>.pkl` – kept reactions per parameter set
- run folders in `3_model/RUNS/Reduced_<date>_<time>/` (from `setup_reduced_run_folders.py`)
- `plots/RAB_fractions/` – plots of kept species fraction vs. ε

Results per parameter set are collected in
`dic_reduced[powplug][epsilon][delta]` with `kept_species`,
`eliminated_species` and the number of kept reactions.

## Modules in `functions/`

| Module | Content |
|--------|---------|
| `create_matrix.py` | `create_matrix`, `load_rate_vector`, `compute_RTS`, `select_species`, `filter_reactions`, `filter_reactions_with_delta`, `plot_RAB_fraction_vs_epsilon`: core of the DRGEP method |
| `select_timepoints_one_puls.py` | `select_timepoints_plateau`: selection of time points from E/N for a single pulse |
| `change_dat.py` | `ChangeDat`: writes a reduced `.dat` file from the complete one and a list of kept reactions |
| `compare_reduced_dat.py` | `create_data_reduced_dic`: summary of one reduction (kept/eliminated species, number of reactions) |
| `setup_reduced_run_folders.py` | `prepare_reduced_run_folders(source_dat_dir, target_base_dir)`: creates Model run folders for all reduced `.dat` files; executed separately |

## Assumptions and pitfalls
- The rate indices in the pickle refer to the complete `.dat` file. Never
  use a different or already reduced `.dat` as input.
- The reduction is based on one run (`run_idx`), i.e. one POWPLUG value.
  Check that `run_idx` points to the intended value; transfer to other
  energies must be validated in step 5.
- Near the stability threshold, small changes in ε cause large error
  changes, so ε should be varied in small steps there.
- The run number `n` in `Run_<n>_...` determines the order in which
  `4_evaluation/evaluate_run.py` assigns the `dic_reduced` entries to the
  runs. This only works if `epsilon_list` and `delta_list` in `DRGEP.py`
  are sorted in ascending order and `modified_dat/` contains exactly the
  mechanisms of the saved `dic_reduced`.
- Run folders and `.dat` files must be sorted numerically by ε/δ, not
  lexicographically.
