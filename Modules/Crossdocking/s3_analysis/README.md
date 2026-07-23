# Stage 3: state-selectivity analysis

## Summary

The `s3_analysis` module addresses state selectivity from two complementary
directions.

### Branch 1: pre-filter population overlap

The PPS- and PR-directed AHC runs are independent generative experiments, each
optimized against a different myosin receptor state. This branch identifies
canonical molecular identities that were generated in both raw populations
before property or diversity filtering.

Shared molecules indicate convergence of the two independent searches onto the
same chemical identity. Comparing their original PPS and PR docking scores
provides an initial estimate of state preference under the original AHC scoring
protocol. Generation frequency and score distributions are retained so that a
single unusually favourable occurrence is not interpreted in isolation.

This is primarily a population-level convergence and replication analysis. It
does not replace cross-docking because each score originates from its own
independent AHC run.

### Branch 2: double-arm cross-docking selectivity

This is the primary state-selectivity analysis. Filtered ligands generated
against one myosin state are docked into the opposite-state receptor grid:

```text
PPS population -> PR receptor grid
PR population  -> PPS receptor grid
```

For each ligand, its original own-state AHC docking score is paired with the
new opposite-state cross-docking score. This tests whether the ligand retains a
preference for the receptor state that originally selected it or also binds
similarly well to the alternative state.

The two branches are complementary: population overlap asks whether the
independent searches rediscovered the same molecules, whereas double-arm
cross-docking tests receptor-state preference for the selected candidates.

## Selectivity definition

Docking scores are more favourable when they are more negative.

For the double-arm cross-docking analysis:

```text
selectivity_margin = opposite_state_score - own_state_score
```

- Positive margin: the original/own receptor state is preferred.
- Negative margin: the opposite receptor state is preferred.
- Larger absolute margin: greater separation between the two state scores.

For a canonical molecule shared by the raw PPS and PR populations:

```text
population_margin = PR_score - PPS_score
```

- Positive margin: PPS has the more favourable score.
- Negative margin: PR has the more favourable score.

These differences are relative Glide docking-score comparisons under matched
protocols. They should not be interpreted as experimentally measured binding
free-energy differences.

## Code architecture

```text
s3_analysis/
├── population_overlap.py
│   └── molecule_identity.py
├── run_analysis.py
│   ├── validation.py
│   └── selectivity.py
├── __init__.py
└── README.md
```

```text
Branch 1: raw-population overlap

PPS raw scores.csv ─┐
                    ├─> molecule_identity.py ─> population_overlap.py
PR raw scores.csv ──┘

Branch 2: double-arm cross-docking

double_arm_manifest.csv ───────────┐
                                   ├─> validation.py ─> selectivity.py
double_arm_best_docking_scores.csv ┘                         │
                                                             └─> run_analysis.py
```

## Branch 1 usage: pre-filter population overlap

`molecule_identity.py` parses and canonicalizes the raw SMILES, applies
technical validity checks, and aggregates repeated occurrences within each
population. `population_overlap.py` joins the PPS and PR populations by
canonical isomeric SMILES and compares their score distributions.

The identity definition:

- preserves formal charge;
- preserves defined stereochemistry;
- ignores SMILES atom-order differences;
- does not neutralize molecules;
- does not perform tautomer normalization.

Technically unusable rows are excluded, including invalid or unparsable SMILES
and missing, zero, or positive docking scores.

```bash
python s3_analysis/population_overlap.py \
  --pps-scores /path/to/PPS/scores.csv \
  --pr-scores /path/to/PR/scores.csv \
  --output-dir /path/to/population_overlap \
  --selectivity-threshold 2.0
```

Principal outputs:

```text
population_overlap_summary.csv
overlapping_molecules.csv
PPS_population_aggregated.csv
PR_population_aggregated.csv
PPS_internal_duplicates.csv
PR_internal_duplicates.csv
PPS_canonicalized_instances.csv
PR_canonicalized_instances.csv
excluded_instances.csv
```

`overlapping_molecules.csv` contains the best, median, mean, standard
deviation, occurrence count, and available generation-step information for
each canonical molecule found in both populations.

## Branch 2 usage: double-arm cross-docking

`validation.py` checks the combined double-arm outputs before score pairing.
`selectivity.py` constructs per-ligand own-state/opposite-state comparisons,
classifies state preference, and ranks the candidates. `run_analysis.py`
coordinates validation, analysis, and output writing.

Inputs produced by `double_arm_docking.py`:

```text
double_arm_manifest.csv
double_arm_best_docking_scores.csv
```

Run:

```bash
python s3_analysis/run_analysis.py \
  --input-dir /path/to/double_arm_crossdock \
  --output-dir /path/to/double_arm_crossdock/analysis \
  --selectivity-threshold 2.0 \
  --strong-score-threshold -8.0
```

Thresholds control categorical classification only. The continuous docking
scores and margins are always retained.

Outputs:

```text
qc_report.csv
paired_state_scores.csv
selectivity_ranking.csv
selectivity_summary.csv
```

The current analysis does not yet include scaffold enrichment, interaction
fingerprints, or binding-pose inspection. These structure-based analyses can
be added after a complete double-arm result has passed QC.
