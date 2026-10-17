# 1_create_data_and_preparation – Preparation of the DRGEP input

Reads the results of the model runs with the full mechanism and converts
them into the two input files for DRGEP:

- `CombinedRuns_<date>_<time>.pkl` – study pickle
- `CombinedRuns_<date>_<time>_reaction_labels.pkl` – reaction labels

The exact format of both files is described in `2_drgep/README.md`.

## Contents

| File / folder | Purpose |
|---------------|---------|
| `param_study_analysis.py` | reads model runs, creates study pickles, combines runs, analysis plots |
| `functions/` | helper modules (see below) |
| `Runs/` | model results, one folder per rate study |

## Why runs have to be combined
DRGEP needs the rate of every reaction. In the thesis, the number of output
variables per model run was limited to about 2000, so the rates of all
reactions could not be written in a single run. The rate output was split
over several runs with identical conditions, each covering a subset of
species and reaction sides. `combine_runs()` merges these runs into one
study that contains the rates of all reactions.

With a model that writes all rates in one run, combining is not necessary;
the study pickle of that run can be used directly.

## `param_study_analysis.py`

### Configuration
Settings at the top of the script:

| Variable | Meaning |
|----------|---------|
| `rate_study_path` | folder in `Runs/` that contains the runs of one rate study |
| `runs` | list of run folder names to process |
| `plot_rates_flag` | `True`: plot reaction rates (slow for many reactions) |

### Modes
The script decides from the run name what to do:

| Run name / situation | Action |
|----------------------|--------|
| `CombinedRuns_<date>_<time>` | `combine_runs()` merges the individual runs of the rate study into one study |
| ≥ 2 entries in `runs` | creates `<run>.pkl` for each run with `load_save_study()` if it does not exist yet |
| one entry, other name | looks for the run directly in `Runs/` |

### Procedure
1. Set `runs` to all run folders of the rate study and run the script.
   One study pickle is created per run.
2. Set `runs = ['CombinedRuns_<date>_<time>']` and run it again.
   The combined study pickle and the reaction labels pickle are created.
3. Copy both files into `2_drgep/` or set `study_path` and `labels_path`
   in `2_drgep/DRGEP.py`.

```bash
cd 1_create_data_and_preparation
python param_study_analysis.py
```

### Analysis performed for every study
- reduced electric field `EN_Td` (added if missing; required by DRGEP)
- SEI, CO₂ conversion and energy efficiency per run
- reaction rate plots (if `plot_rates_flag = True`)
- species ratios for `O, O2, CO, C`
- overview plots (`PlotOverview`), saved in the run folder

### Output (in `Runs/<rate study>/<run>/`)
- `<param>_<date>_<time>.pkl` – study pickle
- `<param>_<date>_<time>_reaction_labels.pkl` – reaction labels
- overview plots

## Reading model results (`load_save_study`)

`load_save_study(study_dir, dat_path)` reads all runs of one study folder,
creates the study pickle and the reaction labels pickle
(format of both: `2_drgep/README.md`).

### Expected folder structure
```
Runs/<rate study>/
└── <PARAM>_<YYYYMMDD>_<HHMM>/          e.g. POWPLUG_20260715_1849
    ├── Run_1_2000.0/                   Run_<index>_<value>
    │   ├── <model output>.tec
    │   └── <model input>.nam
    ├── Run_2_20000.0/
    └── ...
```
- Study folder: parameter name(s) separated by `_`, then date and time
  (`POWPLUG_XPOW_20260626_1038` for two parameters).
- Run folder: `Run_<index>_<value>` or `Run_<index>_[<v1>, <v2>]`.
  Other subfolders are skipped.
- The file names of the model output are set in `load_study()`
  (`tec_path`, `nam_path`).

### Time-resolved results (`.tec`, Tecplot format)
Schematic structure:
```
TITLE = "..."                 line 1: title, ignored
"T (s)"                       one line per variable, name in quotes
"E"
"CO2"
...
"RATE___859"
ZONE I = 10000                last header line: number of time steps
                              empty line = end of header
<values>                      all values of variable 1, then variable 2, ...
```
- Values are separated by spaces and written variable by variable (block
  format), not row by row.
- Fortran number formats are converted: `2.00D+06` and `3.252-160`
  (exponent without `E`).
- If a run is shorter than given in the header, the number of time steps is
  calculated from the number of values.
- The column `MIS` is skipped. `CONVERSION` (CO₂ conversion of all CO₂
  states) and `Tv1` (vibrational temperature) are added.

### Input parameters (`.nam`, Fortran namelist)
Lines `NAME = value` or `NAME = v1, v2, ...` are read into a dict; lines
starting with `!`, `$` or `*` are skipped. Several values are stored as a
list.

### Reaction labels
The labels are created from the `RATE_` columns of run 1 and the mechanism
file `dat_path` with `get_rate_label(..., include_comment=True)`. The
mechanism file must be the one the model used for these runs.

To use another model, it is usually easier to create both pickles directly
in the format of `2_drgep/README.md` than to adapt the parser.

## Modules in `functions/`

| Module | Content |
|--------|---------|
| `safe_run_pickle.py` | `load_save_study()` reads run folders into a study, `save_study()` saves it |
| `combine_pkl.py` | `combine_runs()` merges the runs of one rate study into one study |
| `functions_analysis.py` | analysis: SEI, energy efficiency, E/N (`add_EN_to_study`), ratios, rate plots, `load_study_cache`, `DeviationAnalyzer` |
| `gk_rate_mapper.py` | assigns rate columns to the reactions of the mechanism file (see below) |
| `model_change_parameter.py` | reading and writing model input parameters, Fortran number format (`format_fortran`) |
| `plot_overview.py` | `PlotOverview`: overview plots of a study |

### Rate mapping (`gk_rate_mapper.py`)
Rate columns are named only by number (`RATE___859`). The number is the
position of the reaction in the mechanism file (see `2_drgep/README.md`).

| Function | Purpose |
|----------|---------|
| `build_reaction_lookup(dat_path)` | `{n: (reaction, label)}` |
| `get_rate_label(rxn_map, n, max_len=40, include_comment=False)` | short label for plots, e.g. `CO2+M>CO+O+M  [N01]` |
| `print_mapping(rxn_map, indices=None)` | table number → reaction → label |
| `count_species_occurrences(dat_path, species_list)` | number of reactions per species as educt/product; useful for splitting the rate output over several runs |

```bash
python functions/gk_rate_mapper.py ../2_drgep/dat-file/model_input.dat
```

The mapping is only valid for the mechanism file that was used in the
model run. For reduced runs, the reduced mechanism file must be used.

## Assumptions and pitfalls
- Run labels must be sorted numerically, not lexicographically.
- Pickle files can be truncated by interrupted writes; delete and recreate
  them if loading fails.
