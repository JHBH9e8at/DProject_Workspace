# Cross-docking modules

This directory contains the reusable stages for cross-docking selected AHC
Glide poses against the opposite myosin-state receptor grid.

## Layout

```text
Crossdocking/
├── runner.py                 # sequential terminal pipeline runner
├── single_dock.py            # standalone one-ligand Glide test
├── s1_Prep/
│   ├── fetch_target.py       # filtered CSV -> target manifest
│   ├── unpack_p201.py        # selected SDFGZ -> SDF
│   └── prepare_ligand.py     # optional LigPrep route for new ligands
├── s2_Dock/
│   ├── make_glide_input.py   # create Glide input
│   ├── run_glide.py          # run Schrodinger Glide
│   └── extract_scores.py     # extract output pose scores
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

## Standalone single docking

`single_dock.py` accepts one prepared `.sdf` or `.sdfgz` ligand and one Glide
grid. An SDFGZ input is unpacked automatically. The script runs without Dask
and does not run LigPrep.

```bash
python single_dock.py \
  --ligand /path/to/ligand.sdfgz \
  --grid /path/to/receptor_grid.zip \
  --output-dir /path/to/single_dock_result \
  --precision SP
```

Use a new output directory for each run. Results include the Glide pose file,
`docking_scores.csv`, `single_docking.in`, and `glide_console.log`.

## Sequential pipeline runner

`runner.py` connects target fetching, SDFGZ unpacking, Glide input generation,
docking, and score extraction. Each ligand is isolated under `jobs/`, while
`target_manifest.csv` tracks progress and `combined_docking_scores.csv`
collects all completed poses.

```bash
python runner.py \
  --input /path/to/filtered_PPS.csv \
  --run-dir /path/to/PPS_AHC_run \
  --run-name PPS \
  --grid /path/to/PR_grid.zip \
  --output-dir /path/to/PPS_to_PR_crossdock \
  --precision SP
```

### additional features 
Use `--limit 10` for a small test. If source files are known to be incomplete,
`--allow-unavailable` skips the missing entries. After an interrupted run, use
the same arguments with `--resume` to reuse completed and partial outputs.
`--fail-fast` allows to stop the run at first error
