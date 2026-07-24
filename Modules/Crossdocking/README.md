#### info! 
If you wish to see the –help for running please jump down to the **Software architecture** section below ! 
# Cross-Docking Workflow
This package evaluates whether molecules generated and optimized against the two receptor states under the same AHC protocol may exhibit selectivity for a specific state
## Biological question
The workflow compares two receptor states, **PPS** and **PR**. Separate AHC generative runs produce populations optimized against each state:
- the PPS population was generated and scored against the PPS receptor;
- the PR population was generated and scored against the PR receptor.
A strong score in the receptor state used during generation is not, by itself, evidence of state selectivity. A molecule may bind strongly to both states. The relevant question is therefore:
> Does a molecule that docks strongly to its own generation state dock substantially less strongly to the opposite receptor state?
## Working hypothesis
If state-specific binding environments are meaningfully different, the best-scoring molecules from a state-directed AHC population should show a systematic preference for their own state.
For a molecule selected from state \(A\), the pipeline defines:

```text
selectivity_margin = opposite_state_score - own_state_score
```
Glide docking scores are interpreted such that more negative values indicate
more favorable predicted binding. Therefore:
- a **+ve margin** supports preference for the molecule's own state;
- a margin near zero indicates similar docking strength in both states;
- a **-ve margin** indicates stronger predicted docking to the opposite state.

This is a **computational prioritization hypothesis**, not direct proof of biochemical selectivity. Docking scores should be interpreted together with pose quality, interactions, chemical liabilities, and ultimately experimental measurements.

## Experimental design
The experiment contains two symmetric cross-docking arms:
```text
PPS AHC population ── select strong PPS binders ── dock into PR grid
PR  AHC population ── select strong PR binders  ── dock into PPS grid
```
The original AHC docking score supplies the **own-state score**. Cross-docking supplies the **opposite-state score**. The two values are paired by molecule and used to calculate the selectivity margin.
### Arm definitions
| Arm | Ligand source | Own state | Cross-docking receptor |
|---|---|---|---|
| `L_PPS_to_R_PR` | PPS AHC population | PPS | PR |
| `L_PR_to_R_PPS` | PR AHC population | PR | PPS |

## Population filtering and candidate selection
Before cross-docking, each raw AHC score table is filtered independently. The
filtering block:

1. retains rows marked as valid;
2. retains rows marked as unique;
3. parses and canonicalizes SMILES with RDKit;
4. requires a finite, nonzero own-state docking score;
5. requires a usable best-variant identifier;
6. removes canonical duplicates, retaining the most favorable own-state score;
7. sorts molecules by own-state docking score;
8. selects the configured top percentage by own-state docking score;
9. optionally replaces percentage selection with a test-mode top-N override;
10. annotates PAINS and BRENK structural alerts without automatically removing
   flagged molecules.
`best_dscore_top_per` is applied **after filtering, canonical deduplication, and score ranking**. For example, `best_dscore_top_per=1` selects the best 1% of each usable population. The selected count is calculated independently for PPS and PR with `ceil(population_size × percentage / 100)`, so equally defined percentiles can contain different numbers of molecules when population sizes differ.
When `testmode_top_n` has a value, it overrides percentage selection for both populations. For example, `testmode_top_n=10` selects the ten best usable own-state binders from each population. Leave it blank for production percentage selection.
PAINS and BRENK columns are warning annotations. They are intended for downstream review and do not currently act as hard exclusion rules.

## Selectivity interpretation
The default analysis thresholds are:
```text
selectivity_threshold = 2.0
strong_score_threshold = -8.0
```
The score-analysis block pairs own-state and opposite-state results and assigns selectivity classes. In practical terms:
- `selective_own_state` identifies molecules whose own-state preference meets
  the configured selectivity threshold;
- `nonselective_strong` identifies molecules with strong docking but
  insufficient separation between the two receptor states.
The complete classification logic is implemented in `analysis/score/selectivity.py`. Thresholds are operational prioritization criteria and should be chosen with awareness of the docking protocol and the desired false-positive tolerance.

## Raw-population overlap
When both `pps_scores` and `pr_scores` are supplied, the score block also compares the complete raw AHC populations after canonicalization.
Simultaneously occurring molecules in both the PPS and PR states are analysed in the in an attempt to answer the same question in a different angle:
- **Cross-docking selectivity** compares receptor-state scores for selected
Molecules.
- **Population overlap** measures whether the independent PPS and PR AHC runs
  generated the same canonical chemical structures.
Low population overlap supports population-level separation, but it does not
by itself prove receptor-state selectivity. The paired cross-docking results
remain the direct test used by this workflow.
## Interaction analysis scope
`analysis/interaction/` provides a separate, single-complex interaction block
for:
- protein–ligand interaction fingerprints;
- receptor-state residue mapping;
- pose-quality checks;
- optional clash and strain thresholds.
This block is **not currently invoked automatically** by `crossdocking_runner.py`. It requires a prepared receptor structure and selected docked poses and is intended for focused follow-up analysis of prioritized complexes.

# Module Architecture
## Layout

```text
Crossdocking/
├── Crossdocking_main_readme.md       # scientific and software overview
├── crossdocking_runner.py            # complete config-driven workflow
├── testing_crossdock.in              # validated limit=10 example configuration
├── t_crossdock.in                    # additional configuration template/example
├── prep/
│   ├── __init__.py
│   ├── filter_population.py          # raw AHC CSV -> ranked filtered population
│   ├── fetch_targets.py              # filtered CSV -> target manifest
│   ├── unpack_sdfgz.py               # selected AHC SDFGZ -> SDF
│   ├── prepare_ligands.py            # optional LigPrep preparation
│   └── prep_runner.py                # independently runnable preparation block
├── docking/
│   ├── __init__.py
│   ├── glide_input.py                # create Glide input
│   ├── glide_execution.py            # execute Schrodinger Glide
│   ├── scores.py                     # extract docking pose scores
│   └── docking_runner.py             # independently runnable docking block
├── pipelines/
│   ├── __init__.py
│   ├── single_case.py                # standalone one-ligand docking test
│   ├── single_arm.py                 # one population -> opposite-state grid
│   └── double_arm.py                 # run and combine both cross-docking arms
├── analysis/
│   ├── score/
│   │   ├── __init__.py
│   │   ├── identity.py               # canonical identity and aggregation
│   │   ├── population_overlap.py     # raw-population overlap analysis
│   │   ├── validation.py             # double-arm result QC
│   │   ├── selectivity.py            # paired metrics, classes, and ranking
│   │   ├── run_analysis.py           # score-analysis implementation entry
│   │   └── score_runner.py           # independently runnable score block
│   └── interaction/
│       ├── __init__.py
│       ├── input_models.py           # complex input validation/materialization
│       ├── fingerprint.py            # ProLIF interaction fingerprints
│       ├── pose_quality.py           # PoseCheck-based quality assessment
│       ├── receptor_mapping.py       # PPS/PR residue correspondence
│       ├── single_complex.py         # single-complex implementation
│       └── interaction_runner.py     # independently runnable interaction block
├── tempresults/                      # local copied results from a test run
└── x_patching/                       # development code and test fixtures
```

`x_patching/` is intentionally outside the production workflow and may be
managed or removed separately.

## End-to-end execution flow

`crossdocking_runner.py` is the top-level runner. It accepts one flat
`key=value` configuration file and performs:

```text
raw PPS/PR score tables
        │
        ▼
population filtering and top-N selection
        │
        ▼
two-arm cross-docking
        │
        ├── PPS ligands -> PR grid
        └── PR ligands  -> PPS grid
        │
        ▼
double-arm QC and score reconciliation
        │
        ▼
paired selectivity ranking
        │
        └── optional raw-population overlap
```

The interaction-analysis block is intentionally separate from this flow.

## Running the full workflow

Run from the repository root so that package imports resolve consistently:

```bash
cd /path/to/Crossdockingtemp
python crossdocking_runner.py --config /path/to/testing_crossdock.in
```

The configuration parser accepts:

```text
key=value
```

As susual the `#`ed lines are ignored.
Any Duplicate keys, missing required keys, invalid switches, or unsupported docking precision values 
are intended flag and stall the workflow, but please double check the inputs with **cuation** as the code has not been fully reviewd.

### Required configuration keys
| Key | Meaning |
|---|---|
| `pps_input` | Raw PPS AHC score CSV |
| `pps_run_dir` | PPS AHC result directory containing selected ligand files |
| `pr_input` | Raw PR AHC score CSV |
| `pr_run_dir` | PR AHC result directory containing selected ligand files |
| `pps_grid` | Prepared PPS Glide receptor-grid ZIP |
| `pr_grid` | Prepared PR Glide receptor-grid ZIP |
| `output_dir` | Root directory for all workflow outputs |

### Common optional keys

| Key | Default | Meaning |
|---|---:|---|
| `pps_run_name` | `PPS` | PPS column/file prefix |
| `pr_run_name` | `PR` | PR column/file prefix |
| `step_col` | `step` | AHC step column |
| `pps_variant_col` | inferred | PPS best-variant column |
| `pr_variant_col` | inferred | PR best-variant column |
| `smiles_col` | `smiles` | SMILES column |
| `limit` | all | Top molecules selected per population after filtering |
| `precision` | `SP` | Glide precision: `HTVS`, `SP`, or `XP` |
| `ligand_prep_mode` | `off` | Reuse AHC pose files (`off`) or prepare from SMILES (`on`) |
| `allow_unavailable` | `off` | Permit missing/unusable targets in a manifest |
| `resume` | `off` | Reuse eligible existing job results |
| `fail_fast` | `off` | Stop immediately after an arm failure |
| `selectivity_threshold` | `2.0` | Minimum own-state score advantage |
| `strong_score_threshold` | `-8.0` | Strong-docking classification threshold |
| `pps_scores`, `pr_scores` | unset | Raw CSV pair enabling population-overlap analysis |
| `pps_score_col` | `PPS_r_i_docking_score` | PPS raw own-state score column |
| `pr_score_col` | `PR_r_i_docking_score` | PR raw own-state score column |
| `valid_col` | `valid` | Raw validity column used by overlap analysis |

`pps_scores` and `pr_scores` must either both be supplied or both be omitted.

**Note** : Precision options HTVS, XP are not tested and `may cause errors` !


## Ligand preparation modes

### `ligand_prep_mode=off`

The workflow reuses the selected `*_lib.sdfgz` ligand pose from the AHC run.
The best-variant identifier and AHC directory structure must resolve to a
nonempty source file.

This mode preserves the chemical variant selected during the original AHC run
and docks that variant against the opposite receptor state. It is therefore the
preferred mode for the primary controlled cross-docking experiment when the
scientific question is whether the **same AHC-selected ligand variant** retains
its docking preference across PPS and PR.

This mode was used in the validated test described below. Each molecule in that
test supplied one prepared SDF record to its cross-docking job.

### `ligand_prep_mode=on`

The workflow starts from the filtered CSV SMILES and runs LigPrep before docking. The current LigPrep command uses Epik at pH 7.0 ± 1.0 and permits up
to eight stereoisomers. Depending on the input molecule, LigPrep may produce different numbers of protonation states, tautomers, stereoisomers, and  other prepared variants.

This variable expansion introduces a potential **unequal-enumeration bias**. Glide docks all generated records, while the current molecule-level result retains the most negative `r_i_docking_score` across all returned variants and poses. A molecule with more generated variants therefore receives more opportunities to produce an extreme favorable score than a molecule with only one variant. Best-only scores from molecules with substantially different variant counts are not strictly equivalent comparisons.

For this reason, `ligand_prep_mode=on` should not be treated as directly interchangeable with the validated `off` protocol. It is most appropriate for:

**[the text regarding this section is underrevision]**
- preparing a ligand when the original AHC pose file is unavailable;
- deliberately exploring plausible protonation, tautomeric, or stereochemical
  states;
- sensitivity analysis following the primary fixed-variant comparison.

For a rigorous LigPrep-enabled comparison:

1. use identical LigPrep settings for both receptor-state arms;
2. dock the same prepared variant set for a molecule against both grids;
3. retain `prepared_variant_count`, pose counts, and variant-level scores;
4. report the best score together with the number of variants and poses tried;
5. avoid interpreting unadjusted best-only rankings without checking whether
   enumeration counts differ systematically between molecules or populations;
6. consider a fixed variant cap, a predefined state-selection rule, or a
   variant-count-matched sensitivity analysis.

The pipeline currently records `prepared_variant_count`, but its primary
molecule-level ranking still uses the best docking score. Results generated
with `ligand_prep_mode=on` must therefore be interpreted with this limitation
until an enumeration-aware aggregation policy is implemented.

## Independently runnable blocks

The top-level runner is the normal production entry point. Individual block
runners are also provided for debugging, partial reruns, and method
development:

```bash
# Population preparation
python -m prep.prep_runner --help

# One prepared ligand against one grid
python -m docking.docking_runner --help

# One filtered population against one receptor grid
python -m pipelines.single_arm --help

# Both cross-docking arms
python -m pipelines.double_arm --help

# Score QC, selectivity, and optional population overlap
python -m analysis.score.score_runner --help

# One selected receptor–pose complex
python -m analysis.interaction.interaction_runner --help
```

These block runners receive all information required for their own scope via
arguments. They do not replace `crossdocking_runner.py` as the full workflow
entry point.

## Output structure

```text
output_dir/
├── filtering/
│   ├── PPS_filtered.csv
│   ├── PPS_filtered_filterlog.txt
│   ├── PR_filtered.csv
│   └── PR_filtered_filterlog.txt
├── crossdocking/
│   ├── L_PPS_to_R_PR/
│   │   ├── target_manifest.csv
│   │   ├── combined_docking_scores.csv
│   │   ├── best_docking_scores.csv
│   │   └── jobs/
│   ├── L_PR_to_R_PPS/
│   │   ├── target_manifest.csv
│   │   ├── combined_docking_scores.csv
│   │   ├── best_docking_scores.csv
│   │   └── jobs/
│   ├── double_arm_manifest.csv
│   ├── double_arm_summary.csv
│   ├── double_arm_docking_scores.csv
│   └── double_arm_best_docking_scores.csv
├── analysis/
│   ├── selectivity/
│   │   ├── qc_report.csv
│   │   ├── paired_state_scores.csv
│   │   ├── selectivity_ranking.csv
│   │   └── selectivity_summary.csv
│   └── population_overlap/
│       ├── population_overlap_summary.csv
│       ├── overlapping_molecules.csv
│       ├── excluded_instances.csv
│       └── state-specific canonicalized and aggregated tables
└── pipeline_run_summary.json
```

`pipeline_run_summary.json` records start/end timestamps, active block, resolved
output paths, completion status, configuration path, and any terminal error.
It is the first file to inspect when diagnosing an interrupted run.

## Validated integration test

The complete workflow was validated on **2026-07-24** using
`testing_crossdock.in`.

### Test conditions

```text
PPS raw scores:
/home/andy/proj/Results/temptesting/Testingarea1/Rawoutputs/PPS_scores.csv

PR raw scores:
/home/andy/proj/Results/temptesting/Testingarea1/Rawoutputs/PR_306_raw.csv

PPS AHC run:
/home/andy/proj/Results/AHC_Resultsbin/2026_06_02_SMILES-RNN_PPS_tester_run_2

PR AHC run:
/home/andy/proj/701/output/2026_07_07_SMILES-RNN_PR_run

PPS grid:
/home/andy/proj/DProject_Workspace/AHC_Related/ref_structures/grids/PPS/PPSgrid.zip

PR grid:
/home/andy/proj/DProject_Workspace/AHC_Related/ref_structures/grids/PR/PR_Auto_Grid.zip

limit=10
precision=SP
ligand_prep_mode=off
allow_unavailable=off
resume=off
fail_fast=on
selectivity_threshold=2.0
strong_score_threshold=-8.0
```

### Filtering results

| Population | Raw rows | Usable unique canonical molecules | Selected |
|---|---:|---:|---:|
| PPS | 22,784 | 20,420 | 10 |
| PR | 19,584 | 17,035 | 10 |

The selected PR set contained seven BRENK annotations and no PAINS matches.
The selected PPS set contained neither PAINS nor BRENK matches. These alerts
were retained as annotations.

### Docking and QC results

| Arm | Completed | Failed | Unavailable |
|---|---:|---:|---:|
| `L_PPS_to_R_PR` | 10 | 0 | 0 |
| `L_PR_to_R_PPS` | 10 | 0 | 0 |

All 12 selectivity QC checks passed, including:

- required-column validation;
- unique cross-docking identifiers;
- correct arm/state mapping;
- completed-manifest to best-score reconciliation;
- finite and negative docking scores;
- consistency between manifests and best-score tables.

### Selectivity results

- all 10 PPS-derived molecules were classified as `selective_own_state`;
- 7 of 10 PR-derived molecules were classified as `selective_own_state`;
- 3 PR-derived molecules were classified as `nonselective_strong`;
- the best observed PPS-selective margin was `8.99358` for molecule
  `301_44-1` (`PPS=-15.8853`, `PR=-6.89172`);
- the raw canonical populations shared 3 molecules.

The test completed with:

```text
status=completed
arm_failures=0
analyzed_molecules=20
```

Local copies of the validated outputs are stored under `tempresults/`.

## Recommended run review

For every production run, review at minimum:

1. `pipeline_run_summary.json` — confirm `status=completed`;
2. both filtering logs — confirm the filtering counts and structural alerts;
3. `double_arm_summary.csv` — confirm completed/failed/unavailable counts;
4. `analysis/selectivity/qc_report.csv` — resolve any non-PASS item;
5. `analysis/selectivity/selectivity_ranking.csv` — inspect paired scores and
   prioritized candidates;
6. selected docking poses — verify plausible geometry before biological
   interpretation;
7. interaction analysis for the final prioritized complexes.

## Important limitations

- Docking-score differences are model-derived prioritization signals, not measured binding free energies.
- The top-N sets are conditioned on the upstream AHC search and filtering rules; they do not represent exhaustive chemical space.
- Reusing AHC poses with `ligand_prep_mode=off` assumes that the referenced ligand file is the intended chemical variant.
- Structural alerts are currently annotations, not exclusion criteria.
- Population overlap and paired cross-docking answer related but distinct questions and should not be conflated.
- Interaction analysis is not yet part of the automatic end-to-end runner.
