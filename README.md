# Mechanism Reduction of CO₂ Plasma Chemistry with DRGEP
Master's thesis, Niklas Meier, University of Leipzig / EnaDyne, 2026
Contact: nm@enadyne.de

This repository contains the implementation of the DRGEP method
(Directed Relation Graph with Error Propagation) that was used to reduce a
reaction mechanism for CO₂ plasma chemistry, together with the scripts used
to prepare the input data and to evaluate the reduced mechanisms.

## Scope of this repository

The reaction rates and the validation runs in the thesis were calculated
with a 0D plasma kinetics model. **All parts of this simulation model have
been removed**: source code, executable, input files and all scripts that
only control the model. Wherever the model is referred to in the code, it is
called *model* (e.g. `model_change_parameter.py`, `model_input.dat`).

The repository therefore contains only
- the DRGEP reduction (`2_drgep/`),
- the functions that prepare model results as input for DRGEP
  (`1_create_data_and_preparation/`),
- the functions that evaluate the reduced mechanisms against the full one
  (`3_evaluation/`).

The code is intended mainly as a reference for how DRGEP was implemented.
It can be used with any kinetics model, as long as the model results are
converted into the input format described in `2_drgep/README.md`.
`2_drgep/dat-file/model_input.dat` is included as an example of the
mechanism format and as an example input for the reduction.

## Requirements
- Python 3.13.14, packages in `requirements.txt`

## Installation
```bash
python -m venv mec_reduction
source mec_reduction/bin/activate      # Windows: mec_reduction\Scripts\activate
pip install -r requirements.txt
```

## Folder structure
```
1_create_data_and_preparation/   model results -> study pickle (DRGEP input)
2_drgep/                         DRGEP reduction of the mechanism
3_evaluation/                    comparison of reduced and full mechanism
README.md
requirements.txt
.gitignore
```
Each folder has its own README with details on the code.
`functions/` in `1_...` and `3_...` contain the same modules.

## Workflow

```
 model runs (full mechanism, all reaction rates)          [not included]
        │
        ▼
 1_create_data_and_preparation/param_study_analysis.py
        │  study pickle + reaction labels pickle
        ▼
 2_drgep/DRGEP.py
        │  reduced mechanism files (.dat)
        ▼
 2_drgep/functions/setup_reduced_run_folders.py
        │  one run folder per reduced mechanism
        ▼
 model runs (reduced mechanisms)                          [not included]
        │
        ▼
 3_evaluation/param_study_analysis.py  ->  3_evaluation/evaluate_run.py
           report: errors, species/reaction count, recommendation
```

1. **Model runs with the full mechanism.** DRGEP needs the rate of every
   reaction at the selected time points. In the thesis, the model could not
   write all rates in one run, so the output was split into several runs
   that each cover a subset of reactions.
2. **Preparation** (`1_create_data_and_preparation/`). The runs are read and
   combined into one study pickle and one reaction labels pickle.
3. **Reduction** (`2_drgep/`). `DRGEP.py` removes species with the threshold
   ε and reactions with the threshold σ and writes one reduced mechanism
   file per (ε, σ). `functions/setup_reduced_run_folders.py` then creates
   one run folder per reduced mechanism.
4. **Model runs with the reduced mechanisms**, under the same conditions as
   the full mechanism. Reaction rates are not needed in this step.
5. **Evaluation** (`3_evaluation/`). All reduced runs are compared with the
   reference run (ε = 0, σ = 0): error in CO₂ and CO, NRMSE of the electron
   density, number of species and reactions.

## Input for DRGEP (summary)
Details and examples in `2_drgep/README.md`.

| File | Content |
|------|---------|
| mechanism file (`.dat`) | species part, `*` line, one reaction per line: `A + B > C + D : <rate parameters> % <label>` |
| study pickle (`.pkl`) | `study["runs"][idx]["tec"]`: DataFrame with `T (s)`, `EN_Td` and one column `RATE_…_<n>` per reaction |
| reaction labels pickle (`_reaction_labels.pkl`) | dict `{rate column: "A+B>C+D  [ID]"}` in the same order as the rate columns |

Reaction *n* is the *n*-th active reaction line of the mechanism file.
The unit of the rates does not matter, because DRGEP only uses ratios.

## Thresholds
- ε (epsilon): species threshold
- σ (sigma): reaction threshold

In parts of the code σ is still called `delta`.

## Data structures
- Study: `study["runs"][idx]` with keys `tec` (DataFrame with the
  time-resolved results), `nam` (dict with the model input parameters),
  `param_value`
- Reduction results: `dic_reduced[powplug][epsilon][sigma]` with
  `kept_species`, `eliminated_species`, `kept_reactions`
  (`powplug` = power density of the run used for the reduction)

## Known issues

- The assignment of the DRGEP results (`dic_reduced`) to the reduced runs in
  `evaluate_run.py` is based on the order of the runs. Run folders must be
  sorted numerically by ε/σ, not lexicographically
  (`setup_reduced_run_folders.py` does this).
- `dic_reduced` must be saved as `dic_reduced_*.pkl` in the folder of the
  reduced runs, otherwise the reduction tables in the report are missing.
- Pickle files depend on the pandas/Python version and can be truncated
  by interrupted writes; delete and recreate them if loading fails.
- Under Windows/WSL2, case sensitivity of folder names can cause
  `FileNotFoundError` despite the file existing.
