# Mechanism Reduction of CO₂ Plasma Chemistry with DRGEP
Master's thesis, Niklas Meier, University of Leipzig / EnaDyne, 2026
Contact: [TODO e-mail]

This repository contains the workflow used to reduce a Model reaction
mechanism (CO₂ plasma chemistry) with the DRGEP method and to evaluate the
reduced mechanism against the full one.

## Requirements
- Python 3.13.14, packages in `requirements.txt`
- GNU Fortran (gfortran) 15.2.0 on Ubuntu [TODO: Ubuntu version, native or WSL2]
- Model is included in `3_global_kin/Model/`. This repository is for
  internal use on the EnaDyne server only; Model must not be redistributed
  outside EnaDyne.

## Installation
```bash
python -m venv mec_reduction
source mec_reduction/bin/activate      # Windows: mec_reduction\Scripts\activate
pip install -r requirements.txt
```
Model: [TODO compile instructions, where the executable must be placed]

## Folder structure
```
1_create_data_and_preparation/   simulate all reaction rates, combine runs
2_drgep/                         DRGEP reduction of the mechanism
3_global_kin/                    Model and run folders
4_evaluation/                    evaluation of the reduced mechanism
```
Note: `param_study_analysis.py` and the files in `functions/` exist in both
`1_...` and `4_...` and are identical. The script combines several Model
runs into one pickle file. It is used both for preparing the DRGEP input and
for analysing the reduced runs.

## Workflow

### 1. Simulate all reaction rates
Model cannot output the rates of all reactions in a single run: the number
of plot variables per run is limited by `IPLOTVARDIM` (about 2000).
`1_create_data_and_preparation/Model_run.py` therefore splits the
mechanism into several runs and starts them automatically.

### 2. Combine the runs
Run `param_study_analysis.py` to merge all runs into one pickle file.
- Set the run folder in [TODO variable].
- Output: `CombinedRuns_<date>_<time>.pkl` and `..._reaction_labels.pkl`
- An example for CO₂ chemistry with one pulse is included in `Runs/[TODO name]`.

### 3. DRGEP reduction
Copy both pickle files into `2_drgep/` and run `DRGEP.py`.
- The complete mechanism must be located at `dat-file/global_kin_complete.dat`.
- Set the reduction parameters ε and δ in [TODO file/variable].
- Output: reduced `.dat` files in `dat-file/modified_dat/`.

### 4. Simulate the reduced mechanism
Copy the reduced `.dat` into `3_global_kin/Model_reduced_run/` [TODO check]
and start the simulation with [TODO script]. The rates do not need to be
simulated again, which saves a lot of time.
- `Model_change_parameter/` is for rate simulation (step 1)
- `Model_reduced_run/` is for runs with the reduced mechanism (this step)

### 5. Evaluate
1. Run `4_evaluation/param_study_analysis.py` to combine the reduced runs.
2. Run `4_evaluation/evaluate_run.py` to compare against the full mechanism.
- Output: [TODO plots/tables, where]

## Data structures
- `study["runs"][idx]` with keys `tec` (DataFrame), `nam` (dict), `param_value`
- `dic_reduced[powplug][epsilon][delta]` with `kept_species`,
  `eliminated_species`, `kept_reactions`

## Known issues
- [TODO, e.g. sorting of run labels, species count discrepancies]
- Under Windows/WSL2, case sensitivity of folder names can cause
  `FileNotFoundError` despite the file existing.
- Pickle files depend on the pandas/Python version and can be truncated
  by interrupted writes.