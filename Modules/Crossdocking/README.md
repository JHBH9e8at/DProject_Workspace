# Cross-docking workflow

This module tests whether molecules generated against the PPS or PR receptor
state retain a docking preference for their original state when docked against
the opposite receptor.

## User-facing commands

Run the sequential Glide calculation from the repository root:

```bash
python Modules/Crossdocking/run_crossdocking.py \
  --config Modules/Crossdocking/crossdocking.in
```

Analyse an existing completed calculation:

```bash
python -m Modules.Crossdocking.CompletedAnalysis.crossdocking_analysis \
  --config Modules/Crossdocking/crossdocking.in
```

Scripts under `blocks/` are implementation units. Files under `bins/` are
legacy or development material and are not part of the production workflow.

## Experimental design

```text
L_PPS_to_R_PR: PPS-generated ligands → PR receptor grid
L_PR_to_R_PPS: PR-generated ligands  → PPS receptor grid
```

The original AHC score is the own-state score. The new Glide result is the
opposite-state score:

```text
selectivity_margin = opposite_state_score − own_state_score
```

Glide scores are more favourable when more negative. A positive margin
therefore supports preference for the molecule's original generation state.
This is a computational prioritisation signal, not experimental proof of
selectivity.

## Population filtering

Each AHC population is processed independently:

1. retain rows marked valid;
2. retain rows marked unique;
3. parse and canonicalise SMILES with RDKit;
4. require a finite, nonzero own-state docking score;
5. require a usable best-variant identifier;
6. retain the most negative score for each canonical molecule;
7. rank by own-state score;
8. select the configured top percentage;
9. annotate PAINS and BRENK alerts without excluding matches.

Production selection uses:

```text
N_selected = ceil(N_usable × best_dscore_top_per / 100)
```

`testmode_top_n` overrides percentage selection when it contains a positive
integer. Leave it blank for production runs.

## Required configuration

| Key | Meaning |
|---|---|
| `pps_input`, `pr_input` | Completed PPS/PR AHC score tables |
| `pps_run_dir`, `pr_run_dir` | AHC run directories used to locate selected poses |
| `pps_grid`, `pr_grid` | Prepared Glide receptor-grid ZIP files |
| `output_dir` | Cross-docking calculation output directory |
| `results_root` | Validated project `work_results` directory |
| `run_id` | Identifier shared with completed analysis |

Important options include:

| Key | Current value | Meaning |
|---|---:|---|
| `best_dscore_top_per` | `1` | Best own-state percentage selected after filtering |
| `testmode_top_n` | blank | Optional testing override |
| `precision` | `SP` | Glide precision; HTVS and XP are not validated here |
| `ligand_prep_mode` | `on` | Run LigPrep from SMILES (`on`) or reuse AHC pose (`off`) |
| `allow_unavailable` | `off` | Whether unresolved targets may remain unavailable |
| `resume` | `off` | Whether eligible existing job outputs may be reused |
| `fail_fast` | `on` | Stop when a docking arm fails |

## Ligand preparation modes

With `ligand_prep_mode=off`, the selected AHC `*_lib.sdfgz` variant is reused.
This preserves the original chemical variant for the receptor-state comparison.

With `ligand_prep_mode=on`, LigPrep starts from filtered SMILES and may generate
multiple protonation, tautomeric, or stereochemical variants. Glide docks all
generated records and the molecule-level result retains the most negative
score. Molecules with more enumerated variants consequently receive more
opportunities to obtain an extreme score. Interpret best-only rankings together
with `prepared_variant_count` and `docked_pose_count`.

## Execution and multiprocessing

Glide jobs always run sequentially to avoid licence-slot exhaustion and
unstable concurrent external jobs. The `multiprocessing`,
`posecheck_workers`, `posecheck_chunk_size`, and `prolif_workers` settings apply
only to completed interaction analysis.

With `posecheck_chunk_size=auto`:

```text
chunk_size = max(1, ceil(valid_pose_count / posecheck_workers))
```

## Calculation outputs

```text
<output_dir>/
├── filtering/
│   ├── PPS_filtered.csv
│   ├── PPS_filtered_filterlog.txt
│   ├── PR_filtered.csv
│   └── PR_filtered_filterlog.txt
├── crossdocking/
│   ├── L_PPS_to_R_PR/
│   ├── L_PR_to_R_PPS/
│   ├── double_arm_manifest.csv
│   ├── double_arm_summary.csv
│   ├── double_arm_docking_scores.csv
│   └── double_arm_best_docking_scores.csv
├── pipeline_run_summary.json
└── crossdocking_calculation_manifest.json
```

Inspect `pipeline_run_summary.json` first after an interruption. A completed
calculation must also contain the hashed calculation manifest consumed by the
score-analysis workflow.

## Completed analysis

The completed-analysis runner can perform:

- score validation, pairing, selectivity classification, and population overlap;
- PoseCheck pose-quality analysis;
- ProLIF analysis of cross-docked poses;
- ProLIF analysis of the bundled PPS and PR reference complexes.

Reference analysis generates:

```text
work_results/analysis/reference_interactions/<run_id>/reference_interactions.csv
```

The user does not prepare this file manually. Cross-docked interaction tables
are written below:

```text
work_results/analysis/crossdocking_interactions/<run_id>/
```

## Minimum review

For every run, confirm:

1. `pipeline_run_summary.json` reports `status=completed`;
2. filtering counts and PAINS/BRENK annotations are plausible;
3. both arms have the expected completed target counts;
4. calculation and completed-analysis manifests exist;
5. selectivity QC contains no unresolved failure;
6. selected poses are geometrically plausible;
7. LigPrep variant counts are reviewed when preparation is enabled.

## Limitations

- Docking-score differences are model-derived signals, not binding free energies.
- Top-percent candidates are conditioned on upstream AHC sampling and filtering.
- PAINS and BRENK matches are annotations, not hard exclusions.
- LigPrep enumeration can affect best-score comparisons.
- Population overlap and paired cross-docking answer different questions.
- Pose quality and interactions should be reviewed before prioritising candidates.
