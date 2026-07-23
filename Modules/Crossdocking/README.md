# Cross-docking modules

This directory contains the reusable stages for cross-docking selected AHC
Glide poses against the opposite myosin-state receptor grid.

## Layout

```text
Crossdocking/
├── double_arm_docking.py     # run both state-cross-docking arms
├── single_arm_docking.py     # run one population against one grid
├── single_case_docking.py    # standalone one-ligand Glide test
├── s1_Prep/
│   ├── fetch_target.py       # filtered CSV -> target manifest
│   ├── unpack_p201.py        # selected SDFGZ -> SDF
│   └── prepare_ligand.py     # optional LigPrep route for new ligands
├── s2_Dock/
│   ├── make_glide_input.py   # create Glide input
│   ├── run_glide.py          # run Schrodinger Glide
│   └── extract_scores.py     # extract output pose scores
├── s3_analysis/
│   ├── validation.py         # double-arm input and score QC
│   ├── selectivity.py        # paired state metrics and ranking
│   └── run_analysis.py       # analysis CLI
└── x_patching/               # development code and test fixtures
```

## AHC pose cross-docking route

```text
filtered/selected scores CSV
    -> s1_Prep/fetch_target.py
    -> target manifest CSV
    -> s1_Prep/unpack_p201.py
    -> s2_Dock/make_glide_input.py
    -> s2_Dock/run_glide.py
    -> s2_Dock/extract_scores.py
```

The source `*_lib.sdfgz` files are read in place. `fetch_target.py` records
their paths in a manifest and does not copy them.

## New-ligand route

`s1_Prep/prepare_ligand.py` is retained for SMILES/SDF inputs that have
not already passed through LigPrep. Existing AHC `*_lib.sdfgz` poses use the
unpack route and should not be sent through LigPrep again by default.

## Standalone single-case docking

`single_case_docking.py` accepts one prepared `.sdf` or `.sdfgz` ligand and one Glide
grid. An SDFGZ input is unpacked automatically. The script runs without Dask
and does not run LigPrep.

```bash
python single_case_docking.py \
  --ligand /path/to/ligand.sdfgz \
  --grid /path/to/receptor_grid.zip \
  --output-dir /path/to/single_dock_result \
  --precision SP
```

Use a new output directory for each run. Results include the Glide pose file,
`docking_scores.csv`, `single_docking.in`, and `glide_console.log`.

## Single-arm pipeline

`single_arm_docking.py` connects target fetching, SDFGZ unpacking, Glide input generation,
docking, and score extraction. Each ligand is isolated under `jobs/`, while
`target_manifest.csv` tracks progress and `combined_docking_scores.csv`
collects all completed poses.

```bash
python single_arm_docking.py \
  --input /path/to/filtered_PPS.csv \
  --run-dir /path/to/PPS_AHC_run \
  --run-name PPS \
  --grid /path/to/PR_grid.zip \
  --output-dir /path/to/PPS_to_PR_crossdock \
  --ligand-prep-mode off \
  --precision SP
```

Use `--limit 10` for a small test. If source files are known to be incomplete,
`--allow-unavailable` skips the missing entries. After an interrupted run, use
the same arguments with `--resume` to reuse completed and partial outputs.

### Ligand preparation mode

`--ligand-prep-mode off` is the default and preserves the original workflow:
the selected AHC `*_lib.sdfgz` file is unpacked and docked directly.

`--ligand-prep-mode on` starts instead from the `smiles` column of the filtered
CSV, writes one `.smi` per selected instance, runs Schrodinger LigPrep, and
docks every prepared record produced for that instance. Use `--smiles-col` if
the CSV uses a different column name.

Both modes retain every extracted Glide pose in
`combined_docking_scores.csv`. The lowest (most negative)
`r_i_docking_score` across all prepared variants and poses for each selected
instance is written to `best_docking_scores.csv`. Each job also contains
`docking_scores.csv` and `best_docking_score.csv`.

## Double-arm cross-docking

The complete experiment contains two independent arms:

```text
PPS AHC population -> filtering -> docking to PR receptor grid
PR AHC population  -> filtering -> docking to PPS receptor grid
```

`double_arm_docking.py` runs both arms sequentially by calling the established
single-arm pipeline. It does not rename or copy the original `*_lib.sdfgz`
files. Instead, output folders and combined tables identify the ligand and
receptor states explicitly.

```bash
python double_arm_docking.py \
  --pps-input /path/to/filtered_PPS.csv \
  --pps-run-dir /path/to/PPS_AHC_run \
  --pr-input /path/to/filtered_PR.csv \
  --pr-run-dir /path/to/PR_AHC_run \
  --pps-grid /path/to/PPS_grid.zip \
  --pr-grid /path/to/PR_grid.zip \
  --output-dir /path/to/double_arm_crossdock \
  --ligand-prep-mode off \
  --precision SP
```

For a small test, `--limit 10` applies the limit independently to each arm.
The options `--allow-unavailable`, `--resume`, and `--fail-fast` are passed to
both single-arm runs. `--ligand-prep-mode on|off` and `--smiles-col` are also
applied identically to both arms.

```text
double_arm_crossdock/
├── L_PPS_to_R_PR/
├── L_PR_to_R_PPS/
├── double_arm_manifest.csv
├── double_arm_summary.csv
├── double_arm_docking_scores.csv
└── double_arm_best_docking_scores.csv
```

Combined tables retain the original `molecule_id` and add `ligand_state`,
`receptor_state`, `crossdock_arm`, and a collision-safe `crossdock_id` such as
`L_PPS_to_R_PR__2_3-2`.

## Double-arm selectivity analysis

After both docking arms are complete, build one paired own-state/opposite-state
row per molecule with:

```bash
python s3_analysis/run_analysis.py \
  --input-dir /path/to/double_arm_crossdock \
  --output-dir /path/to/double_arm_crossdock/analysis \
  --selectivity-threshold 2.0 \
  --strong-score-threshold -8.0
```

The core metric is:

```text
selectivity_margin = opposite_state_score - own_state_score
```

Because more negative docking scores are better, a positive margin indicates
preference for the molecule's original/own state. Outputs are
`qc_report.csv`, `paired_state_scores.csv`, `selectivity_ranking.csv`, and
`selectivity_summary.csv`.



# Testing status
260722 : Single arm testing complete, log files


# Future updates 

Change Input format to bash .in file \
Change Ligprep log out dir

