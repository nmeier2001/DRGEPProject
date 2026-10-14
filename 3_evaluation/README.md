# 4_evaluation – Evaluation of the reduced mechanisms

Step 5 of the workflow (see main README). This folder compares the Model
runs with the reduced mechanisms (step 4) against the complete mechanism and
creates a Markdown report with tables, error plots and a recommendation for
the best reduced mechanism.

## Contents

| File / folder | Purpose |
|---------------|---------|
| `param_study_analysis.py` | reads the reduced runs, creates the study pickle, electron density deviation (`DeviationAnalyzer`), overview plots |
| `evaluate_run.py` | creates the run report (Markdown + CSV + plots) |
| `functions/` | same modules as in `1_create_data_and_preparation/functions` (see README there), plus `combine_pkl.py` |
| `Runs/` | [TODO: currently not used by the scripts, see "Paths" below] |

## Paths

Both scripts use the **current working directory** as base
(`BASE_DIR = Path.cwd()`) and read the runs from
`1_create_data_and_preparation/Runs/`. They must therefore be started from
the repository root:

```bash
cd <repository root>
python 4_evaluation/param_study_analysis.py
python 4_evaluation/evaluate_run.py <RUN_NAME>
```

The reduced runs are created in step 4 in `3_model/RUNS/Reduced_<date>_<time>/`.
Copy or move this folder to `1_create_data_and_preparation/Runs/` before
the evaluation. [TODO: check / unify]

## Input
A run folder `Reduced_<date>_<time>/` containing
- one subfolder per reduced mechanism with the Model output (step 4),
- **a reference run with the complete mechanism** (ε = 0, δ = 0). Without it
  no errors can be calculated,
- `dic_reduced_*.pkl` with the DRGEP results from step 3
  (`dic_reduced[powplug][epsilon][delta]`). Optional, but required for the
  reduction tables, the σ plots and the recommendation.

---

## Step 5a: `param_study_analysis.py`

Same script as in step 2, see `1_create_data_and_preparation/README.md`.
For the evaluation set:
```python
runs = ['Reduced_<date>_<time>']
```
For a run name starting with `Reduced`, the script creates the study pickle
`<RUN>.pkl` (if missing) and additionally runs `DeviationAnalyzer` on the
electron density `E` with the run `param_value = 0.000e+00` as reference.

## Step 5b: `evaluate_run.py`

### Run
```bash
python 4_evaluation/evaluate_run.py <RUN_NAME>
python 4_evaluation/evaluate_run.py <RUN_NAME> --reason "test eps=0.5 at 2e6"
python 4_evaluation/evaluate_run.py      # asks for run name and reason
```

| Argument | Meaning | Default |
|----------|---------|---------|
| `run` | name of the run folder in `Runs/` or absolute path | interactive |
| `--reason` | reason for the run, written into the report and `run_history.md` | interactive |
| `--co2-threshold` | max. mean CO₂ error for the recommendation [%] | 2.0 |
| `--nrmse-threshold` | max. NRMSE of the electron density for the recommendation [%] | 5.0 |

The study pickle (`param_study_analysis.py`) must exist before.

### What the script does
1. Loads the study pickle of the run (all `*.pkl` except `dic_reduced*`;
   preferably `<RUN>.pkl`) and adds E/N if missing.
2. Loads all `dic_reduced*.pkl` and builds the reduction tables
   (species/reactions per ε, δ; eliminated species per ε).
3. Finds the reference run: the run whose `param_value` is ≈ 0.
4. Calculates the errors of every run against the reference:
   - CO₂ and CO: mean relative deviation over time,
     `mean((ref − run) / ref) · 100`
   - electron density: NRMSE over the whole time and in an evaluation
     window (`DeviationAnalyzer`)
5. Creates plots and CSV tables.
6. Recommends the run with the fewest reactions among all runs with
   CO₂ error ≤ `--co2-threshold` and NRMSE(E) ≤ `--nrmse-threshold`.
7. Flags runs with identical error values; this usually means that two ε/δ
   values resulted in the same reduced mechanism.

### Output (in the run folder)
```
<RUN>/
├── Report/
│   ├── Run_Report_<RUN>.md               report
│   ├── run_table.csv                     all values per run
│   └── run_history.md                    log of all report runs with reason
├── Error_Analysis/
│   ├── species_error.png                 CO2/CO error vs. ε
│   ├── species_error_vs_sigma.png        CO2/CO error vs. δ
│   ├── nrmse_electron_density.png        E(t) and relative deviation
│   ├── nrmse_vs_epsilon.png
│   └── nrmse_vs_sigma.png
└── Reduction_Overview/                   only with dic_reduced*.pkl
    ├── species_vs_reactions.png
    ├── error_vs_species.png
    ├── error_vs_reactions.png
    ├── reduction_table.csv
    └── eliminated_species_by_epsilon.csv
```
Plots are only created if the required data exist; the console output
shows which plots were skipped and why.

### Report structure
1. General information (run, date, varied parameter, reference run, reason)
2. Parameter values per run (POWPLUG, XPOW, CO₂ conversion, CO production)
3. Mechanism reduction (DRGEP): targets, ε/δ table, eliminated species
4. Error analysis against the reference
5. Recommendation
6. Notes and assumptions

The report uses relative image links and can be opened directly in
Obsidian or any Markdown viewer.

### Definitions
| Quantity | Definition |
|----------|------------|
| CO₂ conversion [%] | `(1 − min(CO2) / CO2(t=0)) · 100` |
| CO production [%] | `(max(CO) − CO(t=0)) / CO2(t=0) · 100` |
| CO₂/CO error [%] | mean relative deviation from the reference over all time steps |
| NRMSE_E_gesamt / _fenster [%] | NRMSE of the electron density over the whole time / evaluation window |

---

## Assumptions and pitfalls
- **Assignment of `dic_reduced` entries to runs is based on order.**
  `dic_reduced` contains no run index; the n-th (ε, δ) entry is assigned
  to the n-th run. If the run folders were sorted differently (e.g.
  lexicographically instead of numerically), ε and δ are assigned to the
  wrong runs. Check the control columns `epsilon (dic_reduced)` and
  `delta (dic_reduced)` in the report.
- The CO₂/CO error is a signed mean: positive and negative deviations at
  different times can cancel each other out. Always look at the plots too.
- If CO₂ starts at zero, conversion and CO production cannot be calculated
  and stay empty.
- Run the scripts from the repository root, otherwise `Path.cwd()` points
  to the wrong folder.
- The docstring of `evaluate_run.py` still refers to the old name
  `create_run_report.py` and to `compare_pkl.py`, which is not part of
  this repository.
