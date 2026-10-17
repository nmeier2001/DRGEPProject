# 3_evaluation – Evaluation of the reduced mechanisms

Compares the model runs with the reduced mechanisms against the run with the
full mechanism and creates a Markdown report with tables, error plots and a
recommendation for the best reduced mechanism.

## Contents

| File | Purpose |
|------|---------|
| `param_study_analysis.py` | reads the reduced runs into a study pickle, deviation of the electron density (`DeviationAnalyzer`), overview plots |
| `evaluate_run.py` | creates the report (Markdown, CSV, plots) |
| `functions/` | same modules as in `1_create_data_and_preparation/functions` (see README there) |

## Input

A folder `Reduced_<date>_<time>/` in `1_create_data_and_preparation/Runs/`

The study pickle of the reduced runs has the same structure as in
`2_drgep/README.md`. Rate columns are not needed; required are `T (s)`,
`CO2`, `CO` and the electron density `E`.

## Paths
Both scripts use the current working directory as base (`Path.cwd()`) and
must be started from the repository root:

```bash
python 3_evaluation/param_study_analysis.py
python 3_evaluation/evaluate_run.py <RUN_NAME>
```

## Step 1: `param_study_analysis.py`

Same script as in `1_create_data_and_preparation/`. Set
```python
runs = ['Reduced_<date>_<time>']
```
For a run name starting with `Reduced`, the script creates the study pickle
`<RUN>.pkl` (if missing) and compares the electron density `E` of all runs
with the reference run (`DeviationAnalyzer`).

## Step 2: `evaluate_run.py`

```bash
python 3_evaluation/evaluate_run.py <RUN_NAME>
python 3_evaluation/evaluate_run.py <RUN_NAME> --reason "test eps=0.5"
python 3_evaluation/evaluate_run.py      # asks for run name and reason
```

| Argument | Meaning | Default |
|----------|---------|---------|
| `run` | name of the run folder or absolute path | interactive |
| `--reason` | reason for the run, written into the report and `run_history.md` | interactive |
| `--co2-threshold` | max. mean CO₂ error for the recommendation [%] | 2.0 |
| `--nrmse-threshold` | max. NRMSE of the electron density for the recommendation [%] | 5.0 |

### What the script does
1. Loads the study pickle of the run and adds E/N if missing.
2. Loads all `dic_reduced*.pkl` and builds the reduction tables.
3. Finds the reference run (`param_value` ≈ 0).
4. Calculates the errors of every run against the reference:
   - CO₂ and CO: mean relative deviation over time,
     `mean((ref − run) / ref) · 100`
   - electron density: NRMSE over the whole time and in an evaluation
     window
5. Creates plots and CSV tables.
6. Recommends the run with the fewest reactions among all runs with
   CO₂ error ≤ `--co2-threshold` and NRMSE(E) ≤ `--nrmse-threshold`.
7. Flags runs with identical errors; this usually means that two ε/σ values
   resulted in the same reduced mechanism.

### Output (in the run folder)
```
<RUN>/
├── Report/
│   ├── Run_Report_<RUN>.md
│   ├── run_table.csv
│   └── run_history.md                    log with reason per report
├── Error_Analysis/
│   ├── species_error.png                 CO2/CO error vs. ε
│   ├── species_error_vs_sigma.png        CO2/CO error vs. σ
│   ├── nrmse_electron_density.png
│   ├── nrmse_vs_epsilon.png
│   └── nrmse_vs_sigma.png
└── Reduction_Overview/                   only with dic_reduced*.pkl
    ├── species_vs_reactions.png
    ├── error_vs_species.png
    ├── error_vs_reactions.png
    ├── reduction_table.csv
    └── eliminated_species_by_epsilon.csv
```
Plots are only created if the required data exist; the console shows which
plots were skipped.

### Definitions
| Quantity | Definition |
|----------|------------|
| CO₂ conversion [%] | `(1 − min(CO2) / CO2(t=0)) · 100` |
| CO production [%] | `(max(CO) − CO(t=0)) / CO2(t=0) · 100` |
| CO₂/CO error [%] | mean relative deviation from the reference over all time steps |
| NRMSE_E [%] | normalised root mean square error of the electron density, whole time / evaluation window |

In parts of the code, σ is still called `sigma` (e.g. column `sigma`).

## Assumptions and pitfalls
- **The `dic_reduced` entries are assigned to the runs by order**, because
  `dic_reduced` contains no run index. Check the columns
  `epsilon (dic_reduced)` and `sigma (dic_reduced)` in the report.
- The CO₂/CO error is a signed mean: positive and negative deviations can
  cancel each other out. Always check the plots as well.
- If CO₂ starts at zero, conversion and CO production stay empty.
- Start the scripts from the repository root, otherwise `Path.cwd()` points
  to the wrong folder.
