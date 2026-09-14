# Analysis and figure workflow

This directory provides one top-level runner for the analysis and figure modules.
Use it when several modules should share one result directory and one run ID.

```bash
python Modules/Analysis_Figen_block/run_analysis.py \
  --config Modules/Analysis_Figen_block/analysis.in
```

The top-level runner does not contain scientific calculations itself. It reads
`analysis.in`, selects the enabled modules, and calls each module's own runner.
Detailed scientific settings remain in the child `.in` files.

## What the top-level options do

| Option | Required | Accepted value | Purpose |
|---|---:|---|---|
| `results_root` | Yes | Path | Root directory for all generated results. Its final directory name must be `work_results`. This value overrides `results_root` in every enabled child config. |
| `run_id` | Yes | Text | Identifier shared by all enabled modules. It groups outputs from the same analysis run. This value overrides `run_id` in every enabled child config. |
| `dry_run` | No | `ON` / `OFF` | With `ON`, validates the top-level configuration and prints the resolved execution plan without running any child module. Default: `OFF`. |
| `ahc_analysis` | No | `ON` / `OFF` | Enables completed AHC-result analysis. Default: `OFF`. |
| `ahc_config` | Conditional | Path to `.in` | Child configuration used when `ahc_analysis=ON`. It may be omitted only when that module is OFF. |
| `structural_analysis` | No | `ON` / `OFF` | Enables pocket-volume and/or redocking analysis according to `structural_config`. Default: `OFF`. |
| `structural_config` | Conditional | Path to `.in` | Child configuration used when `structural_analysis=ON`. |
| `interaction_result_analysis` | No | `ON` / `OFF` | Enables post-processing of existing interaction tables. It does not run PoseCheck or ProLIF. Default: `OFF`. |
| `interaction_result_config` | Conditional | Path to `.in` | Child configuration used when `interaction_result_analysis=ON`. |
| `figures` | No | `ON` / `OFF` | Enables standalone figure generation according to `figures_config`. Default: `OFF`. |
| `figures_config` | Conditional | Path to `.in` | Child configuration used when `figures=ON`. |

At least one of the four module switches must be `ON`. Option names are strict;
unknown keys and invalid ON/OFF values cause an error instead of being ignored.
Relative paths are resolved from the directory containing `analysis.in`.

## Module responsibilities

Modules always execute in the following order, with disabled modules skipped:

| Order | Module | Main purpose | Important prerequisite |
|---:|---|---|---|
| 1 | `ahc_analysis` | Analyse completed PR/PPS AHC score populations: physicochemical properties, trajectory summaries, chemical space, scaffolds, and diversity. | Completed PR and PPS `scores.csv` files. |
| 2 | `structural_analysis` | Calculate selected Fpocket volumes and/or summarise redocking scores and RMSD. | Fpocket outputs and reference ligand for pocket volume; four redocking CSVs for redocking analysis. |
| 3 | `interaction_result_analysis` | Convert existing ProLIF/PoseCheck-derived tables into residue-hotspot and related comparison tables. | Cross-docking completed analysis must already have generated the required interaction tables. |
| 4 | `figures` | Regenerate selected AHC, chemical-space, physicochemical, and interaction figures from existing analysis tables. | Every input required by each enabled figure block must already exist. |

The detailed options for each module are documented beside its child config:

- `ahc_analysis/README.md` and `ahc_analysis/ahc_analysis.in`
- `structural_analysis/README.md` and `structural_analysis/structural_analysis.in`
- `interaction_result_analysis/README.md` and `interaction_result_analysis/interaction_result_analysis.in`
- `figures/README.md` and `figures/figures.in`

## Override rule

Only two values are propagated from the top-level config:

```text
top-level results_root -> overrides results_root in every enabled child config
top-level run_id       -> overrides run_id in every enabled child config
```

All other settings are read exclusively from each child `.in` file. For example,
turning `structural_analysis=ON` does not automatically enable both pocket-volume
and redocking calculations; their individual switches must be set in
`structural_analysis.in`.

This design allows a child module to run independently with its own `run_id`,
while the same config can participate in a combined run under the top-level
`run_id`.

## Recommended preflight and full run

First, check the paths and selected modules without starting calculations:

```ini
results_root=/home/andy/proj/work_results
run_id=PR_PPS_analysis
dry_run=ON
```

```bash
python Modules/Analysis_Figen_block/run_analysis.py \
  --config Modules/Analysis_Figen_block/analysis.in
```

The runner prints a plan similar to:

```text
Analysis run_id: PR_PPS_analysis
Results root: /home/andy/proj/work_results
[1/3] structural_analysis: .../structural_analysis.in
[2/3] interaction_result_analysis: .../interaction_result_analysis.in
[3/3] figures: .../figures.in
```

After confirming the plan, change `dry_run=OFF` and run the same command again.
The modules run sequentially. If one module raises an error, execution stops and
later modules are not started.

## Example module selections

Run structural analysis only:

```ini
ahc_analysis=OFF
structural_analysis=ON
interaction_result_analysis=OFF
figures=OFF
```

Post-process existing interactions and then generate figures:

```ini
ahc_analysis=OFF
structural_analysis=OFF
interaction_result_analysis=ON
figures=ON
```

Run every module:

```ini
ahc_analysis=ON
structural_analysis=ON
interaction_result_analysis=ON
figures=ON
```

Enabling every module does not guarantee that every internal analysis block is
enabled. The ON/OFF switches inside all four child `.in` files must also be
reviewed.

## Current production configuration

The supplied `analysis.in` currently enables:

```text
ahc_analysis=OFF
structural_analysis=ON
interaction_result_analysis=ON
figures=ON
dry_run=OFF
```

This configuration assumes that completed AHC analysis and cross-docking
completed analysis have already produced the tables consumed by later modules.
In particular, `interaction_result_analysis` does not rerun PoseCheck or ProLIF.

## Outputs and run identity

Results are written below the configured `work_results` root. Numerical tables,
metadata, and manifests are normally placed below `work_results/analysis/`, and
figures below `work_results/Figs/`. Each child module creates its own feature
subdirectory and uses the shared top-level `run_id` to identify this invocation.

Do not judge a calculation only from its PNG files. Review the associated CSV,
JSON, and manifest files for exact inputs, counts, settings, and output paths.

## Directory convention

Each executable subdirectory exposes a runner, a configuration file, a README,
and an internal `blocks/` directory. Files in `blocks/` are implementation
details and are not intended to be invoked individually. `common/` contains
shared configuration, path, output, and manifest utilities. Historical or
duplicated implementations under `bins/` are not part of the canonical workflow.
