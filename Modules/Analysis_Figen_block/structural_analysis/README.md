# Structural analysis

Run the enabled structural blocks through one entry point:

```powershell
python run_structural_analysis.py --config structural_analysis.in
```

`pocket_volume` and `redocking` are independently selectable with ON/OFF, but
at least one must be enabled. Implementation files under `blocks/` are internal.

## Pocket volume

`blocks/pocket_volume.py` contains pure parsing, sphere selection, and Monte Carlo
volume calculations. It does not create output files.

`blocks/pocket_volume_runner.py` is the internal output boundary. It writes to:

```text
work_results/analysis/pocket_volume/<run_id>/
|-- inputs/
|-- intermediate/
|-- tables/
|   |-- <output_prefix>.csv
|   `-- <output_prefix>.json
|-- logs/
`-- run_manifest.json
```

Example from the `work_space` root:

```text
python -m Modules.Analysis_Figen_block.structural_analysis.pocket_volume_runner \
  --fpocket-dir <fpocket-result-directory> \
  --target 2,38,51 \
  --ligand <reference-ligand.pdb> \
  --target-label PPS \
  --output-prefix PPS_reference_3A
```

If `--run-id` is omitted, a timestamped run ID is generated. Existing runs
are rejected unless `--resume` is explicit. A resumed run allocates matching
CSV/JSON suffixes rather than overwriting existing scientific results.

The legacy implementation and its old result files are preserved under
`Analysis_Figen_block/bins/legacy_modules/Pocekt_volume_analysis/`. Migration
tests execute both implementations on the same synthetic fixture and require
identical CSV rows and JSON payloads.

## Redocking

When `redocking=ON`, all four PR/PPS docking-score and RMSD CSV paths are
required. Tables are written below `work_results/analysis/redocking/`; PNG/SVG
pairs are written below `work_results/Figs/redocking_validation/`. `score_ylim`
and `rmsd_ylim` each accept two comma-separated bounds.

Relative paths are resolved from `structural_analysis.in`. `resume=OFF` rejects
an existing run; `resume=ON` versions related tables and figures together.
