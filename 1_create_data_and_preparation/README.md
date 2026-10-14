# 1_create_data_and_preparation – Rate simulation and data preparation

Steps 1 and 2 of the workflow (see main README). This folder simulates the
rates of all reactions in the full mechanism with Model and combines the
results into one pickle file, which is the input for the DRGEP reduction in
`2_drgep/`.

`param_study_analysis.py` and `functions/` also exist in `4_evaluation/`
and are used there again to analyse the reduced runs (step 5).

## Contents

| File / folder | Purpose |
|---------------|---------|
| `model_run.py` | Step 1: starts all Model runs needed to output every reaction rate |
| `param_study_analysis.py` | Step 2: reads the runs, creates study pickles, combines runs, analyses and plots |
| `functions/` | Helper modules (see below) |
| `Runs/` | Output of Model, one folder per run |

---

## Step 1: `model_run.py`

### Why several runs are needed
Model only writes reaction rates for species listed in `CRATE_PLOT` in the
`.nam` file. Each entry outputs the rates of all reactions in which that
species occurs on the side given in `ISIDE_PLOT`
(`1` = [TODO: reactant side?], `2` = [TODO: product side?]).
The total number of plot variables per run is limited by `IPLOTVARDIM`
(about 2000). To get the rates of *all* reactions, the species/side
combinations are therefore split over several simulations.

### Configuration
All settings are at the top of the script:

| Variable | Meaning |
|----------|---------|
| `nam_path` | Path to the Model input file `model.nam` |
| `param_name` | List of `.nam` parameters to vary, e.g. `['POWPLUG']` |
| `param_value` | One list of values per entry in `param_name`, e.g. `[[2e4, 2e5, 2e6, 2e7, 2e8]]` |
| `simulations` | Dict of simulations, each with a `CRATE_PLOT` list (species) and an `ISIDE_PLOT` list (side per species) |

Rules for `simulations`:
- `CRATE_PLOT` and `ISIDE_PLOT` must have the same length; entry *i* of
  `ISIDE_PLOT` belongs to entry *i* of `CRATE_PLOT`.
- The number of rates per simulation must stay below `IPLOTVARDIM`.
- Together, all simulations must cover every reaction of the mechanism.
  [TODO: how the split was determined / how to check completeness]
- A species may appear in several simulations or on both sides. Note that the
  number of plot variables is not simply additive in that case, so check
  the limit with a test run.

The current split (13 simulations) is made for `model_complete.dat`
(CO₂ chemistry). If the mechanism changes, the split must be redone.

### What the script does
For every entry in `simulations`:
1. `changing_parameters()` writes `CRATE_PLOT` and `ISIDE_PLOT` into the `.nam` file.
2. `ModelRunner` starts Model once for every value in `param_value`
   and stores the results in `Runs/`.

Total number of Model runs = number of simulations × number of parameter
values (currently 13 × 5 = 65).

### Run
```bash
python model_run.py
```
Output: [TODO: folder naming in `Runs/`, files per run (.tec, .nam, ...)]

---

## Step 2: `param_study_analysis.py`

### Configuration
| Variable | Meaning |
|----------|---------|
| `rate_study_path` | Folder in `Runs/` that contains the runs of one rate study |
| `runs` | List of run folder names to process |
| `plot_rates_flag` | `True`: plot reaction rates (slow for many reactions) |

### Modes
The script decides from the run name what to do:

| Run name / situation | Mode | Action |
|----------------------|------|--------|
| `CombinedRuns_<date>_<time>` | combined run | `combine_runs()` merges the individual runs into one study |
| any name, `runs` contains ≥ 2 entries | rate study run | creates `<run>.pkl` with `load_save_study()` if it does not exist yet |
| `Reduced_*` (single entry) | reduced run | additionally compares the runs with `DeviationAnalyzer` (reference: `param_value = 0.000e+00`, [TODO: what the reference run is]) |
| any other single name | single run | looks for the run directly in `Runs/` |

Typical procedure for step 2 [TODO: check]:
1. Set `runs` to all run folders of the rate study and run the script.
   This creates one study pickle per run.
2. Set `runs = ['CombinedRuns_<date>_<time>']` and run it again.
   This creates the combined pickle for DRGEP.

### Analysis performed for every study
- reduced electric field `EN_Td` (added if missing)
- SEI, CO₂ conversion and energy efficiency per run
- reaction rate plots (`plot_rates`, if `plot_rates_flag = True`)
- species ratios for `O, O2, CO, C`
- pulse-width plot if the varied parameter is `XPOW`
- overview plots (`PlotOverview`), saved in the run folder

### Run
```bash
python param_study_analysis.py
```

### Output (in `Runs/<rate study>/<run>/`)
- `<param>_<date>_<time>.pkl` – study pickle
- `<param>_<date>_<time>_reaction_labels.pkl` – reaction labels for the rate columns
- overview plots [TODO: file names]

For step 3, copy the combined pickle and its `_reaction_labels.pkl` into `2_drgep/`.

---

## Data structure of a study pickle
```python
study = {
    "param_name": "POWPLUG",          # or list of names
    "date": "...", "time": "...",
    "runs": {
        1: {                          # run index starts at 1
            "tec": DataFrame,         # Model output, columns e.g. 'T (s)', 'CO2', 'E', rates
            "nam": dict,              # parsed .nam parameters
            "param_value": "2.000e+06",
        },
        ...
    },
}
```
[TODO: further keys]

## Modules in `functions/`

| Module | Content |
|--------|---------|
| `function_model_run.py` | `ModelRunner`: starts Model for each parameter value, manages run folders |
| `model_change_parameter.py` | `changing_parameters(nam_path, params)`: writes parameters into the `.nam` file |
| `safe_run_pickle.py` | `load_save_study()` reads a run folder (.tec/.nam) into a study and saves it, `save_study()` |
| `functions_analysis.py` | analysis functions: SEI, energy efficiency, E/N, ratios, rate plots, `load_study_cache`, `DeviationAnalyzer` |
| `gk_rate_mapper.py` | maps `RATE___xxx` columns of the `.tec` output to reactions of the `.dat` file; counts in how many reactions a species appears as educt/product (see below) |
| `plot_overview.py` | `PlotOverview`: overview plots of a study |

### Rate mapping (`gk_rate_mapper.py`)
Model names rate columns only by number (`RATE___859`). The number is the
sequential index of the reaction in the `.dat` file: counting starts at 1 after
the `*` line that separates species from reactions, and every non-empty,
non-comment (`!`) line containing `>` is one reaction.

| Function | Purpose |
|----------|---------|
| `build_reaction_lookup(dat_path)` | returns `{index: (reaction, comment)}`; the comment is the text after `%` (e.g. `"N01  Eremin1997,Nist2015"`) |
| `get_rate_label(rxn_map, index, max_len=40, include_comment=False)` | short label for plot legends, e.g. `CO2+M>CO+O+M  [N01]` |
| `print_mapping(rxn_map, indices=None)` | prints a table index → reaction → comment |
| `count_species_occurrences(dat_path, species_list)` | returns `{species: {'educt', 'product', 'total'}}`, the number of reactions per species and side. Useful to estimate the plot variables per `CRATE_PLOT` entry when splitting the simulations in step 1 |

```python
from functions.gk_rate_mapper import build_reaction_lookup, get_rate_label
rxn_map = build_reaction_lookup("model.dat")
get_rate_label(rxn_map, 859, include_comment=True)
```
Quick check from the command line:
```bash
python functions/gk_rate_mapper.py path/to/model.dat
```

Important: the mapping is only valid for the exact `.dat` file used in the
simulation. A reduced mechanism has a different numbering, so for reduced
runs the reduced `.dat` file must be passed, not the complete one.
Species matching in `count_species_occurrences` is exact (`CO2` does not
match `CO2v01` or `CO2^`); a species appearing twice on one side of a
reaction is counted once.

## Assumptions and pitfalls
- `.tec` files use fixed 10-character column names and Fortran number formats
  (D-notation, exponents without `E` such as `3.252-160`); parsing is
  handled in `safe_run_pickle.py` / `functions_analysis.py`.
- Run labels must be sorted numerically, not lexicographically.
- On Linux/WSL2, folder names are case-sensitive (`Runs` ≠ `RUNS`).
- Pickle files can be truncated by interrupted writes; delete and recreate
  them if loading fails.
- Conversion and energy efficiency use the CO₂ value at t = 0. If CO₂
  starts at zero, these values are not meaningful.
