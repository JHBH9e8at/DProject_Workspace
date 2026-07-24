# Analysis Module

The `analysis` module converts double-arm cross-docking results into state-selectivity evidence and provides an independent structure-based follow-up workflow for selected receptor–ligand complexes. 

The module contains two blocks:
```text
analysis/
├── score/          # population identity, score QC, and selectivity
└── interaction/    # pose quality and residue-level interactions
```

`analysis/score` runs automatically as part of the complete cross-docking workflow. `analysis/interaction` is currently an independently runnable single-complex block and is not automatically invoked by `crossdocking_runner.py`.

---

# Part I. Scientific Reasoning

## 1. Analysis questions

The analysis module addresses three related but distinct questions.

### Question A: Do selected molecules prefer their own receptor state?

This is the primary cross-docking question. A molecule selected from one AHC population has:
- an **own-state score** from the original AHC run;
- an **opposite-state score** from the new cross-docking run.
Pairing these two values tests whether a molecule retains preference for the receptor state that originally selected it.

### Question B: Did the independent PPS and PR searches generate the same molecules?

Raw-population overlap measures convergence between the two complete AHC populations. It identifies canonical molecular structures generated in both independent searches and retains generation frequency and score-distribution information.

This is a population-level identity analysis. It does not replace cross-docking because the PPS and PR scores originate from separate AHC runs rather than a controlled redocking experiment.

### Question C: Are the prioritized docking poses structurally credible?

Docking-score separation alone is not sufficient for mechanistic interpretation. The interaction block can evaluate selected complexes for:
- protein–ligand contacts;
- residue-level interaction fingerprints;
- receptor-state residue correspondence;
- steric clashes;
- strain and other available pose-quality metrics.

This provides supporting structural evidence after score QC and selectivity ranking.

## 2. Double-arm selectivity

The two cross-docking arms are:

```text
L_PPS_to_R_PR: PPS-selected ligand -> PR receptor grid
L_PR_to_R_PPS: PR-selected ligand  -> PPS receptor grid
```

Glide docking scores are more favorable when they are more negative. 

For each completed molecule:

```text
selectivity_margin = opposite_state_score - own_state_score
```

Interpretation:
- +ve margin: own state has the more favorable score;
- -ve margin: opposite state has the more favorable score;
- larger absolute margin: greater separation between the two scores.

With the current defaults:

```text
selectivity_threshold = 2.0
strong_score_threshold = -8.0
```

the classification rules are:

| Class | Rule | Interpretation |
|---|---|---|
| `selective_own_state` | margin ≥ 2.0 and own score ≤ -8.0 | strong own-state docking with sufficient separation |
| `own_state_preferred_weak` | margin ≥ 2.0 and own score > -8.0 | separated scores, but weak own-state docking |
| `reverse_selective` | margin ≤ -2.0 and opposite score ≤ -8.0 | strong preference for the opposite state |
| `reverse_preferred_weak` | margin ≤ -2.0 and opposite score > -8.0 | opposite-state preference with weak docking |
| `nonselective_strong` | \|margin\| < 2.0 and both scores ≤ -8.0 | strong docking to both states without sufficient separation |
| `inconclusive` | all remaining cases | insufficient evidence under the configured thresholds |

Thresholds control categorical labels only. Continuous own-state scores, opposite-state scores, and margins are always retained.

These differences are relative docking-score comparisons. They must not be interpreted as experimentally measured binding free-energy differences.

## 3. Ranking policy

Candidates are ranked in the following class order:

```text
selective_own_state
own_state_preferred_weak
nonselective_strong
inconclusive
reverse_preferred_weak
reverse_selective
```

Within a class, candidates are sorted by:

1. larger selectivity margin;
2. more favorable own-state docking score.

The class-priority ranking is designed to prioritize the intended state-selective candidates. It should not obscure scientifically interesting reverse-selective or strong dual-binding molecules, which remain in the output tables with their continuous metrics.

## 4. Raw-population overlap

The raw PPS and PR score tables are processed independently:

```text
raw population
    -> technical validity checks
    -> isomeric SMILES canonicalization
    -> molecule-level aggregation
    -> PPS/PR canonical identity join
```

The identity definition:

- preserves formal charge;
- preserves defined stereochemistry;
- ignores SMILES atom-order differences;
- does not neutralize molecules;
- does not perform tautomer normalization.

Therefore, different charge states, defined stereoisomers, and different tautomers may remain distinct molecular identities.

Technically unusable rows are excluded. Current exclusion categories include:

- invalid validity flag;
- missing or unparsable SMILES;
- missing or nonfinite docking score;
- zero or positive docking score.

For a canonical molecule found in both populations:

```text
population_margin = PR_score - PPS_score
```

- positive margin: PPS score is more favorable;
- negative margin: PR score is more favorable.

Best, median, and mean margins are retained because repeated generation can produce a distribution of scores. Agreement between best and median margins is also reported.

Important: these are original scores from independent AHC experiments. A population-overlap preference is supporting evidence only and is not equivalent to a matched cross-docking selectivity measurement.

## 5. Interaction and pose-quality analysis

The interaction block operates on one prepared receptor and one docked SDF/SDFGZ pose file at a time.

It supports:
- `posecheck`: pose-quality analysis only;
- `prolif`: interaction fingerprints only;
- `both`: both analyses.

ProLIF results can be annotated with a PPS/PR common-residue mapping so that state-specific receptor residue identifiers can be compared in a shared reference system. Secondary contacts such as hydrophobic and van der Waals contacts are optional.

Pose-quality thresholds are optional. When supplied, they can flag results based on:
- total clash count;
- clashes per ligand heavy atom;
- strain energy.

The interaction block does not currently select candidates automatically from `selectivity_ranking.csv`. The user must provide the prepared receptor and pose file for the complex to analyze.

## 6. LigPrep-aware interpretation

`selectivity.py` retains:

- `ligand_prep_mode`;
- `prepared_variant_count`;
- `docked_pose_count`;
- the best cross-docking variant.

This metadata is important when `ligand_prep_mode=on`. Molecules that generate different numbers of LigPrep variants receive different numbers of docking opportunities, while the current molecule-level result uses the most favorable score. Best-only rankings can therefore contain unequal-enumeration bias.

When analyzing LigPrep-enabled runs:

- compare prepared-variant and pose counts;
- retain variant-level results;
- avoid interpreting best-only rankings without checking enumeration balance;
- prefer docking the same prepared variant set against both receptor grids.

---

# Part II. Module Architecture

## 7. Layout

```text
analysis/
├── README.md
├── __init__.py
├── score/
│   ├── __init__.py
│   ├── identity.py               # canonicalize and aggregate raw populations
│   ├── population_overlap.py     # compare PPS and PR raw populations
│   ├── validation.py             # QC double-arm docking outputs
│   ├── selectivity.py            # pair, classify, summarize, and rank scores
│   ├── run_analysis.py           # double-arm QC/selectivity implementation
│   └── score_runner.py           # complete score-block runner
└── interaction/
    ├── __init__.py
    ├── input_models.py           # validate inputs and materialize SDFGZ
    ├── fingerprint.py            # ProLIF interaction fingerprints
    ├── pose_quality.py           # PoseCheck analysis and thresholds
    ├── receptor_mapping.py       # map PPS/PR residue identifiers
    ├── single_complex.py         # single-complex implementation
    └── interaction_runner.py     # interaction-block runner
```

## 8. Score-block data flow

```text
crossdocking/double_arm_manifest.csv ──────────────┐
                                                   ├─> validation.py
crossdocking/double_arm_best_docking_scores.csv ──┘         │
                                                             ▼
                                                       selectivity.py
                                                             │
                                                             ▼
                                                  analysis/selectivity/

raw PPS scores.csv ──┐
                     ├─> identity.py -> population_overlap.py
raw PR scores.csv ───┘                       │
                                            ▼
                               analysis/population_overlap/
```

`score_runner.py` always runs double-arm QC and selectivity. Raw-population overlap runs only when both `--pps-scores` and `--pr-scores` are supplied.

## 9. Running the complete score block

Run from the Crossdocking repository root:

```bash
python -m analysis.score.score_runner \
  --input-dir /path/to/crossdocking \
  --output-dir /path/to/analysis \
  --selectivity-threshold 2.0 \
  --strong-score-threshold -8.0 \
  --pps-scores /path/to/PPS_scores.csv \
  --pr-scores /path/to/PR_scores.csv \
  --pps-score-col PPS_r_i_docking_score \
  --pr-score-col PR_r_i_docking_score \
  --smiles-col smiles \
  --valid-col valid
```

Required inputs:

| Argument | Meaning |
|---|---|
| `--input-dir` | Double-arm cross-docking output directory |
| `--output-dir` | Root analysis output directory |

Optional inputs:

| Argument | Default | Meaning |
|---|---:|---|
| `--selectivity-threshold` | `2.0` | Minimum absolute score separation for a state-preference class |
| `--strong-score-threshold` | `-8.0` | Strong docking threshold |
| `--pps-scores` | unset | Raw PPS score CSV for population overlap |
| `--pr-scores` | unset | Raw PR score CSV for population overlap |
| `--pps-score-col` | `PPS_r_i_docking_score` | PPS raw score column |
| `--pr-score-col` | `PR_r_i_docking_score` | PR raw score column |
| `--smiles-col` | `smiles` | Raw SMILES column |
| `--valid-col` | `valid` | Raw validity column |

`--pps-scores` and `--pr-scores` must be supplied together or both omitted.

## 10. Double-arm QC

Before selectivity pairing, `validation.py` checks:

1. required manifest columns;
2. required best-score columns;
3. unique `crossdock_id` values in both tables;
4. correct arm/state mapping;
5. presence of both expected arms;
6. reconciliation of completed manifest rows with best-score rows;
7. finite opposite-state docking scores;
8. negative opposite-state docking scores;
9. presence of PPS and PR own-state score columns;
10. finite own-state scores;
11. consistency between manifest and best-score values when both are present.

Results are written to `selectivity/qc_report.csv`. Any `ERROR` prevents selectivity analysis from proceeding. A completed run should contain only `PASS` rows.

## 11. Selectivity outputs

```text
analysis/selectivity/
├── qc_report.csv
├── paired_state_scores.csv
├── selectivity_ranking.csv
└── selectivity_summary.csv
```

### `qc_report.csv`

Machine-readable validation report. Each row contains:

- `severity`;
- `check`;
- `affected_count`;
- explanatory `message`.

### `paired_state_scores.csv`

One row per successfully paired molecule. It contains the own-state score, opposite-state score, margin, preferred state, class, preparation metadata, best cross-docking variant, and available ligand-efficiency metrics.

This is the primary table for statistical or custom downstream analysis.

### `selectivity_ranking.csv`

The paired table ordered by classification priority, selectivity margin, and own-state docking strength. `selectivity_rank` is added as a convenience for candidate review.

### `selectivity_summary.csv`

Aggregated by cross-docking arm, own state, and selectivity class. It reports:

- molecule count;
- mean own-state score;
- mean opposite-state score;
- mean selectivity margin;
- median selectivity margin.

## 12. Population-overlap outputs

```text
analysis/population_overlap/
├── population_overlap_summary.csv
├── overlapping_molecules.csv
├── PPS_population_aggregated.csv
├── PR_population_aggregated.csv
├── PPS_internal_duplicates.csv
├── PR_internal_duplicates.csv
├── PPS_canonicalized_instances.csv
├── PR_canonicalized_instances.csv
└── excluded_instances.csv
```

### `population_overlap_summary.csv`

Key/value overview of input counts, valid and excluded instances, unique canonical molecules, overlap size, union size, overlap fractions, Jaccard similarity, shared-molecule score margins, and overlap preference classes.

### `overlapping_molecules.csv`

One row per canonical molecule found in both populations. It contains PPS and PR occurrence counts, score distributions, best source rows, generation steps, variant identifiers, cross-population margins, and overlap classification.

### `PPS_population_aggregated.csv` and `PR_population_aggregated.csv`

One row per unique canonical molecule within each population. These files contain occurrence counts, best/median/mean/worst scores, score dispersion, generation-step range, and the source information for the best occurrence.

### `PPS_internal_duplicates.csv` and `PR_internal_duplicates.csv`

Subsets of the aggregated tables containing molecules with `occurrence_count > 1`. They identify structures rediscovered within one AHC run.

### `PPS_canonicalized_instances.csv` and `PR_canonicalized_instances.csv`

Audit tables at the original instance level. They map each accepted raw SMILES and source row to its canonical isomeric SMILES and retain score and generation metadata.

### `excluded_instances.csv`

Rows rejected during raw-population preparation, including population, source row, exclusion reason, original SMILES, and available score. Use this file to audit coverage losses.

## 13. Running interaction analysis

Run one prepared receptor–pose complex from the repository root:

```bash
python -m analysis.interaction.interaction_runner \
  --protein /path/to/prepared_receptor.pdb \
  --poses /path/to/docked_poses.sdfgz \
  --receptor-state PPS \
  --output-dir /path/to/interaction_result \
  --analysis both \
  --residue-map /path/to/PPS_PR_residue_map.csv \
  --include-secondary-interactions \
  --max-clashes 10 \
  --max-clashes-per-heavy-atom 0.5 \
  --max-strain-energy 20
```

Required arguments:

| Argument | Meaning |
|---|---|
| `--protein` | Prepared receptor PDB or MOL2 |
| `--poses` | Docked ligand SDF or SDFGZ |
| `--receptor-state` | `PPS` or `PR` |
| `--output-dir` | Result directory for this complex |

Optional arguments:

| Argument | Meaning |
|---|---|
| `--analysis` | `posecheck`, `prolif`, or `both` |
| `--residue-map` | Common PPS/PR residue-mapping CSV |
| `--include-secondary-interactions` | Include hydrophobic and van der Waals contacts |
| `--max-clashes` | Optional total-clash threshold |
| `--max-clashes-per-heavy-atom` | Optional normalized-clash threshold |
| `--max-strain-energy` | Optional strain-energy threshold |

The `posecheck` and `prolif` Python dependencies must be available when their
respective analyses are requested.

## 14. Interaction outputs

Depending on `--analysis`, the output directory contains:

```text
interaction_result/
├── intermediate/
│   └── <input_pose_name>.sdf
├── pose_metadata.csv
├── pose_quality.csv          # posecheck or both
├── interactions_long.csv     # prolif or both
└── run_summary.csv
```

### `pose_metadata.csv`

Stable pose indexing and available SDF metadata. Invalid SDF records stop the workflow because they would make pose-to-result mapping unreliable.

### `pose_quality.csv`

PoseCheck results with optional threshold-based quality flags.

### `interactions_long.csv`

Long-format residue interactions joined to pose metadata. When a residue map is provided, receptor-state residue identifiers are annotated with their commonmapping.

### `run_summary.csv`

Records which interaction components completed, their output row counts, the receptor state, and resolved input/intermediate paths.

## 15. Validated score-analysis test

The score block completed successfully in the full integration test performed on **2026-07-24** with ten selected molecules per arm.

```text
L_PPS_to_R_PR: completed 10, failed 0, unavailable 0
L_PR_to_R_PPS: completed 10, failed 0, unavailable 0
Analyzed molecules: 20
Selectivity QC checks: 12 PASS, 0 ERROR
```

Selectivity results:

- all 10 PPS-derived molecules were `selective_own_state`;
- 7 PR-derived molecules were `selective_own_state`;
- 3 PR-derived molecules were `nonselective_strong`;
- the largest margin was `8.99358` for PPS molecule `301_44-1`.

Raw-population results:

```text
PPS input instances:             22,784
PR input instances:              19,584
PPS valid instances:             21,015
PR valid instances:              18,475
PPS unique canonical molecules:  20,420
PR unique canonical molecules:   17,035
Shared canonical molecules:      3
Canonical union:                 37,452
```

The three shared molecules were `inconclusive` at the configured selectivity threshold of 2.0. This result indicates very low structural overlap between the two raw AHC populations but does not, by itself, establish receptor-state selectivity.

The interaction block was not executed as part of this integration test.

## 16. Recommended review order

For each completed analysis:

1. confirm that `qc_report.csv` contains no `ERROR`;
2. inspect `paired_state_scores.csv` for missing or unexpected metadata;
3. review score distributions and margins, not only categorical labels;
4. inspect `prepared_variant_count` and `docked_pose_count`, especially for
   LigPrep-enabled runs;
5. examine `selectivity_ranking.csv` for intended, reverse-selective, and
   strong dual-binding candidates;
6. use population overlap only as complementary convergence evidence;
7. inspect docking poses and run interaction analysis before mechanistic
   interpretation.

## 17. Limitations

- Docking scores are model-derived prioritization metrics, not measured binding
  affinities or free energies.
- Score thresholds are operational classification rules.
- Raw-population overlap compares independent AHC scores and is not controlled
  cross-docking.
- Canonical identity does not normalize tautomeric or charge states.
- Best-only scoring can favor molecules with more prepared variants or poses.
- The interaction block currently processes one receptor–pose file at a time.
- Interaction analysis is not automatically connected to candidate selection.
- Structural plausibility and experimental validation remain necessary.
