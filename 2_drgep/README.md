# 2_drgep – DRGEP mechanism reduction

Reduces a reaction mechanism with the Directed Relation Graph with Error
Propagation (DRGEP) method and writes the reduced mechanisms as new
mechanism files. This is the core of the repository.

Two thresholds are used:
- **ε (epsilon)** – species threshold (DRGEP interaction coefficient),
- **σ (sigma)** – reaction threshold (contribution of a reaction to the
  production or consumption of a kept species).

## Contents

| File / folder | Purpose |
|---------------|---------|
| `DRGEP.py` | main script: species and reaction reduction, writes reduced mechanism files |
| `functions/` | helper modules (see below) |
| `dat-file/model_input.dat` | complete mechanism (input, example) |
| `dat-file/modified_dat/` | reduced mechanisms (output) |
| `plots/` | plots of R_TS vs. ε |

---

## Input

DRGEP needs three files that belong together:

| # | File | Created by |
|---|------|------------|
| 1 | mechanism file `dat-file/model_input.dat` | given |
| 2 | study pickle `CombinedRuns_<date>_<time>.pkl` | `1_create_data_and_preparation/param_study_analysis.py` |
| 3 | reaction labels pickle `CombinedRuns_<date>_<time>_reaction_labels.pkl` | `1_create_data_and_preparation/param_study_analysis.py` |

File 1 is used only to write the reduced mechanisms (`ChangeDat`). The
reduction itself works only with files 2 and 3. To use results of another
model, create files 2 and 3 with the structure below.

### 1. Mechanism file (`model_input.dat`)

Plain text, two parts separated by a line containing only `*`. `!` at the
start of a line marks a comment. Inside a data line, `!` is a field
separator, not a comment.

```
line 1–3     comment lines with the column description
line 4–84    species part, one species per line
line 85      !
line 86      *
line 87–91   comment lines (column description of the reactions)
line 92–     reaction part, one reaction per line
```

**Species part.** One species per line, fields separated by fixed
characters:

```
E           :         0.0; -1]   0.0[     0.0$  0.000546#   1{   0}    !
CO2         :       -4.08;  0]  3.94[   195.2$      44.0#   0{   0}    !NIST webbook.
```

| Separator | Field |
|-----------|-------|
| before `:` | species name |
| `:` … `;` | enthalpy of formation (eV/molecule) |
| `;` … `]` | charge |
| `]` … `[` | Lennard-Jones parameter 1 (Å) |
| `[` … `$` | Lennard-Jones parameter 2 (K) |
| `$` … `#` | molecular weight (amu) |
| `#` … `{` | fraction of flux consumed on surfaces |
| `{` … `}` | branching ratio for consumed flux |
| after `!` | comment / source |

Lines 78–84 contain pseudo-species of the model (`TE`, `TGAS`, `EDEP`,
`PDEP`, `SPEED`, `POSITION`, `MIS`), `M` (line 77) is the third body.
The example file contains 80 active species.

**Reaction part.** One reaction per line:

```
E + CO2 > E + E + CO2^              : 0.0 [ 0.0 ] 0.0 ! -4028 $0&   0. # 0. % X02        LXCat, Phelps
CO + M > C + O + M                  : 4.61E-04 [ -2.86 ] 129000. !  20 $0& 999. # 0. %      Appleton1970 - CO rate
```

| Element | Meaning |
|---------|---------|
| `educts > products` | reaction equation; species separated by `+`, spaces are ignored; a species can appear several times (`E + E`) |
| `:` | end of the equation |
| `A [ n ] E_act` | Arrhenius parameters (pre-factor, temperature exponent, activation energy) |
| `! sp no $ se & dH # ELOSS` | further model parameters (e.g. reference to an electron impact cross section) |
| `%` | start of the label; the first word is the reaction ID (e.g. `X02`) or the source |

The rate parameters are not used by DRGEP; `ChangeDat` copies them
unchanged into the reduced mechanism.

**Reaction numbering.** All non-comment lines with `>` after the `*` line
are numbered 1, 2, 3, … in the order of the file. Commented reactions
(`!` at the start) are not counted. The example file contains 6071 active
reactions. This number *n* links the mechanism file to the rate column
`RATE_…_<n>` in the study pickle.

**Fixed line positions.** `ChangeDat` currently uses fixed line ranges:
species are checked in lines 5–77 (CO2 … M; the electron in line 4 and the
pseudo-species are never removed), reactions from line 92 to the end.
If a mechanism with a different layout is used, the ranges in
`change_species()` and `change_reactions()` of `functions/change_dat.py`
must be adapted.

### 2. Study pickle (`CombinedRuns_<date>_<time>.pkl`)

Python `dict`, saved with `pickle`:

```python
study = {
    "param_name": "POWPLUG",              # varied parameter (str or list of str)
    "date": "20260715",                   # YYYYMMDD
    "time": "0828",                       # HHMM
    "runs": {
        1: {                              # run index (int), starts at 1
            "param_value": "2.000e+06",   # str for one parameter, list of float for several
            "path": Path(...),            # run folder
            "nam":  {...},                # dict of model input parameters
            "tec":  pandas.DataFrame,     # time-resolved results
        },
        2: {...},
    },
}
```

`DRGEP.py` uses only `study["runs"][run_idx]["tec"]` and
`study["runs"][run_idx]["param_value"]`.

Required columns of `tec` (one row per time step, float64):

| Column | Content |
|--------|---------|
| `T (s)` | time in s, ascending |
| `EN_Td` | reduced electric field E/N in Td, must be > 0 (logarithm) |
| `RATE_…_<n>` | rate of reaction *n* of the mechanism file, e.g. `RATE___859`. The column name must start with `RATE_` and end with `_<n>` |

Further columns (species densities `CO2`, `CO`, `E`, …) are ignored by
DRGEP but needed for the evaluation.

- The unit of the rates does not matter, as long as it is the same for all
  reactions: DRGEP only uses ratios of rates. The sign is ignored (`abs`).
- The order of the `RATE_` columns must be the same as the order of the
  keys in the reaction labels pickle (see below).
- A reaction without rate column is not part of the reduction and is
  **commented out** in every reduced mechanism, also for ε = 0.

### 3. Reaction labels pickle (`..._reaction_labels.pkl`)

Python `dict`, one entry per `RATE_` column, in the same order:

```python
labels = {
    "<rate column of reaction 1>": "E+CO2>E+CO2  [X01]",
    "<rate column of reaction 2>": "E+CO2e2>E+CO2e2  [X01[e2]]",
    ...
    "RATE___859": "CO2+M>CO+O+M  [N01]",
    ...
}
```

| Part | Format |
|------|--------|
| key | name of the rate column in `tec` |
| value | `<equation>` + **two spaces** + `[<group>]` |
| `<equation>` | educts and products separated by `>`, species by `+`, without spaces |
| `<group>` | first word after `%` in the mechanism file |

The two spaces are required: `create_matrix()` splits the label with
`split('  ')` and fails if there are none or more than one.
The labels are created by `get_rate_label(rxn_map, n, include_comment=True)`
from the mechanism file (`functions/gk_rate_mapper.py` in `1_...`). This
function shortens equations to 40 characters; all equations in the example
file are shorter (max. 27).

`create_matrix()` builds from the labels
- `species_list`: all species that appear in a labelled reaction, sorted
  alphabetically (incl. `M`),
- `nu`: stoichiometric matrix (species × reactions); duplicate species are
  counted (`E + E` → 2).

The targets must be contained in `species_list`.

---

## Configuration

Settings at the top of `DRGEP.py`:

| Variable | Meaning |
|----------|---------|
| `labels_path`, `study_path` | paths to the two input pickles |
| `epsilon_list` | list of ε values. `0.0` = no species removed (reference) |
| `sigma_list` | list of σ values, only used for ε > 0 |
| `targets` | target species, e.g. `['CO2', 'CO']` |
| `run_idx` | index of the run in the study whose rates are used (e.g. the run at 2·10⁶ W) |
| `plot_kept_species_vs_epsilon` | `True`: additional plots of R_TS (see below) |

For the evaluation in `3_evaluation/`, `epsilon_list` must contain `0.0`.
This creates the reference mechanism, which is simulated as reference run.

## Run

```bash
cd 2_drgep
python DRGEP.py                                  # -> dat-file/modified_dat/*.dat
python functions/setup_reduced_run_folders.py    # -> one run folder per mechanism
```

## How the code works

### 1. Time points – `select_timepoints_plateau(time, EN)`
The rates are evaluated at characteristic time points of one pulse, found
from E/N:

| Key | Time point |
|-----|------------|
| `EN_peak` | absolute maximum of E/N |
| `plateau_1` | 16 time steps after the maximum |
| `2ed peak` | second detected change of E/N (transition between plateaus) |
| `fall` | third detected change (strongest decrease before the last plateau) |
| `afterglow` | 3000 time steps after `fall` |
| `quasi_steady` | last time step |

Changes of E/N are detected with a sliding window (5 steps) over log(E/N):
`score = mean(log EN before) − mean(log EN after)`, then
`scipy.signal.find_peaks` (height 0.03, distance 40 steps).
With `n_points > 0`, additional points are placed logarithmically between
neighbouring time points.

The step offsets (16, 3000) were chosen for one pulse with about 10 000 time
steps. For other time resolutions or pulse shapes they must be adapted, and
at least three changes of E/N must be detected.

### 2. Rate vectors – `load_rate_vector(study_path, time, run)`
Takes the row of `tec` closest to the time point and returns the values of
all `RATE_` columns.

### 3. Direct interaction coefficients – `compute_rAB`
For each time point:

```
r_AB = Σ_i |ν_A,i · ω_i| · δ_B,i  /  Σ_i |ν_A,i · ω_i|
```

ν = stoichiometric coefficient, ω = rate, δ_B,i = 1 if species B takes part
in reaction i. r_AB is the fraction of the turnover of A that depends on
reactions involving B.

### 4. Generalised coefficients – `compute_generalized_R`, `compute_RTS`
From each target T, the path with the largest product of r_AB to every
species S is searched (max-product / widest path, priority queue):

```
R_TS = max over all paths T → S of Π r_AB
```

R_TS of a species is the maximum over all targets and all time points.
Targets have R_TS = 1.

### 5. Species selection – `select_species`
All species with R_TS ≥ ε are kept.

### 6. Reaction selection
- **ε = 0:** `filter_reactions` keeps all reactions whose species are all
  kept, i.e. all labelled reactions.
- **ε > 0:** `filter_reactions_with_sigma`
  1. keeps all reactions whose species (except `M`) are all kept,
  2. for σ > 0 additionally: for each kept species, the share of each
     reaction in its production P and consumption C is calculated from the
     net stoichiometry. A reaction is kept if
     `max over species and time points of max(P, C) ≥ σ`.
     Reactions without net change (e.g. elastic `E + CO2 > E + CO2`) have
     no share and are always removed for σ > 0.

### 7. Writing the reduced mechanism – `ChangeDat`
The kept reactions are saved as `skeletal_mechanism_<ε>_<σ>.pkl` (dict like
`rxn_info`). `ChangeDat` copies the complete mechanism to
`dat-file/modified_dat/model_<ε>_<σ>.dat` and
- comments out (`! `) all species of lines 5–77 that do not occur in a kept
  reaction,
- comments out (`!`) all reactions that are not kept.

Reactions are identified by their equation (species sorted, spaces
removed), not by their number. Reactions with the same equation (52 in the
example file) are therefore always kept or removed together.

### 8. Optional plots – `plot_RAB_fraction_vs_epsilon`
For several runs (`run_indices`): histogram of R_TS and number of kept
species vs. ε, saved in `plots/RAB_fractions/`. Run indices and ε markers
are set in the function call at the end of `DRGEP.py`.

## Output

| File | Content |
|------|---------|
| `dat-file/modified_dat/model_<ε>_<σ>.dat` | reduced mechanism, same format as the input file |
| `skeletal_mechanism_<ε>_<σ>.pkl` | kept reactions (dict like `rxn_info`) |
| `plots/RAB_fractions/*.png` | R_TS histograms, kept species vs. ε |
| `dic_reduced` (in memory) | `dic_reduced[powplug][ε][σ]` with `kept_species`, `eliminated_species`, `kept_reactions`, `powplug`, `epsilon`, plus `dic_reduced[powplug]['targets']`. `3_evaluation/evaluate_run.py` needs it as `dic_reduced_*.pkl` |

## Creating the run folders (`functions/setup_reduced_run_folders.py`)

1. collects all `.dat` files in `dat-file/modified_dat/`,
2. puts the reference (`model_0.0_0.0.dat`) first,
3. sorts all other files **numerically** by ε, then σ,
4. creates one folder `Run_<n>_<suffix>` per file in the target folder and
   copies the mechanism into it.

```
Reduced_20261005_1219/
├── Run_1_0.0_0.0/      reference (all labelled reactions)
├── Run_2_0.1_0.0/
├── Run_3_0.5_0.0/
└── Run_4_0.5_0.1/
```
These folders are simulated with the model and evaluated in `3_evaluation/`.
The run number determines how `evaluate_run.py` assigns the `dic_reduced`
entries, so `epsilon_list` and `sigma_list` must be sorted in ascending
order.

## Modules in `functions/`

| Module | Content |
|--------|---------|
| `create_matrix.py` | `create_matrix`, `load_rate_vector`, `compute_rAB`, `compute_generalized_R`, `compute_RTS`, `select_species`, `filter_reactions`, `filter_reactions_with_sigma`, `plot_RAB_fraction_vs_epsilon` |
| `select_timepoints_one_puls.py` | `sliding_window`, `select_timepoints_plateau` |
| `change_dat.py` | `ChangeDat`: writes the reduced mechanism file |
| `compare_reduced_dat.py` | `create_data_reduced_dic`: summary of one reduction |
| `setup_reduced_run_folders.py` | `prepare_reduced_run_folders`: run folders for the model |

## Assumptions and pitfalls
- The reduction is based on one run (`run_idx`). Transfer to other
  conditions has to be checked in `3_evaluation/`.
- Mechanism files and run folders must be sorted numerically, not
  lexicographically (`_` sorts after digits in ASCII).
