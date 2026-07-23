# Pose-quality and residue-interaction analysis plan

Date: 2026-07-23  
Project: Myosin PPS/PR state-selectivity analysis  
Status: Planned analysis; code not yet implemented

## Purpose

The preceding analysis stage identifies state-selective candidates from Glide
docking-score differences. Docking scores alone, however, do not explain why a
ligand appears to prefer the PPS or PR myosin state. The next analysis stage
will therefore examine the predicted protein-ligand interactions associated
with state preference while explicitly accounting for uncertainty in docking
pose quality.

The proposed workflow has two linked components:

1. a pose-quality assessment, used to determine whether the interaction trends
   remain stable after progressively removing less plausible docking poses;
2. a residue-level interaction fingerprint analysis, used to identify
   interactions enriched in PPS-selective, PR-selective, or nonselective
   molecules.

The resulting structural interpretation remains computational and
hypothesis-generating. It does not establish an experimentally observed binding
mode or experimentally calibrated binding confidence.

## Scientific hypotheses

### Primary hypothesis

Ligands predicted to prefer the PPS and PR myosin states will show different
patterns of residue-level interactions that reflect structural differences
between the two receptor states.

Examples of a state-associated pattern may include:

- formation or loss of a hydrogen bond to a state-dependent residue position;
- retention of an ionic anchor interaction in one state but not the other;
- state-specific aromatic stacking or cation-pi geometry;
- redistribution of hydrophobic contacts caused by pocket-shape differences;
- a recurrent combination of interactions rather than one isolated contact.

### Pose-quality robustness hypothesis

If an interaction pattern provides a meaningful structural explanation of
state preference, its direction and approximate effect size should remain
stable when increasingly strict pose-quality criteria are applied.

An interaction trend observed only among low-quality poses, or one that
disappears immediately after poses with severe clashes or excessive ligand
strain are removed, will be treated as weak evidence.

### Interaction-retention hypothesis

For a state-selective ligand, a characteristic interaction may be retained in
the preferred-state pose but lost, weakened, or geometrically disrupted when
the same ligand is docked into the alternative receptor state.

This will be examined using paired own-state and opposite-state poses from the
double-arm cross-docking experiment.

## Two-stage analytical strategy

### Stage A: hypothesis generation from original AHC poses

The original PPS- and PR-directed AHC docking outputs will be analysed first.
These poses represent the structures used during the independent
state-conditioned generative runs.

The objectives are to:

- construct residue-by-interaction fingerprints for each valid pose;
- compare interaction frequencies between the PPS and PR populations;
- identify interactions enriched in one state-directed population;
- evaluate whether enrichment remains stable across pose-quality thresholds;
- define a restricted set of candidate interaction hypotheses for subsequent
  testing.

This stage is exploratory. Results will be interpreted as population-level
interaction hypotheses rather than confirmation of state selectivity.

### Stage B: validation using double-arm cross-docking poses

The same pose-quality and interaction-extraction protocol will then be applied
to the double-arm cross-docking results:

```text
PPS-origin ligand: PPS original pose versus PR cross-docked pose
PR-origin ligand:  PR original pose versus PPS cross-docked pose
```

The interaction hypotheses selected in Stage A will be fixed before the
primary Stage B comparison. Validation will focus on:

- retention or loss of the hypothesised interaction;
- change in interaction geometry where available;
- association between interaction change and docking selectivity margin;
- consistency across multiple molecules and scaffolds;
- stability of the conclusion across pose-quality thresholds.

Interactions discovered only after examining the Stage B results will be
labelled exploratory rather than confirmatory.

The two stages are not experimentally independent because both rely on docking
poses and related scoring protocols. Where the dataset size permits,
scaffold-level separation between hypothesis-generation and validation
molecules should be considered to reduce chemical-series leakage.

## Pose Quality Index and sensitivity analysis

### Terminology

The planned metric will be called the **Pose Quality Index (PQI)** or
**Pose Plausibility Index**, not a calibrated confidence probability.

Without comparison against experimentally resolved poses, PQI cannot estimate
the probability that a docking pose is correct. It is an internally defined
measure of whether a pose passes selected geometric and energetic plausibility
checks.

### Proposed PQI components

The first implementation should retain the individual quality measurements
rather than immediately hiding them inside one weighted score.

Candidate measurements include:

- protein-ligand steric clash count;
- severe-clash count;
- clash count normalized by ligand heavy-atom count;
- ligand strain energy;
- strain normalized by heavy atoms or rotatable bonds;
- invalid ligand geometry or sanitization failure;
- minimum protein-ligand atomic distance;
- whether the ligand remains inside the intended binding pocket.

Some failures should be treated as hard exclusions, for example unreadable
structures, failed sanitization, or extreme atomic overlap. Other measurements
can be used as continuous or thresholded QC variables.

If a composite PQI is later introduced, its formula, normalization, weights,
and missing-value handling must be reported explicitly. The uncombined
component metrics will always be preserved.

### Threshold sensitivity analysis

Interaction results will be recalculated over progressively stricter
pose-quality subsets. Initial thresholds may be expressed using direct
physical cutoffs, retained fractions, or both:

```text
all technically valid poses
moderate-QC subset
strict-QC subset
very-strict-QC subset
```

At every threshold, the analysis will record:

- number and fraction of retained PPS and PR poses;
- number of retained molecules and scaffolds;
- Glide score and selectivity-margin distributions;
- molecular size, formal charge, and rotatable-bond distributions;
- residue-interaction frequencies and enrichment effect sizes.

A trend will be considered more robust if its direction is preserved and its
effect size remains reasonably stable as QC stringency increases. Statistical
significance alone will not be used as the stability criterion.

Thresholding can introduce composition bias. For example, larger or more
flexible ligands may fail strain or clash filters more frequently. Therefore,
state-specific retention rates and changes in molecular-property distributions
must be reported alongside the interaction trends.

## Tool selection

### PoseCheck: primary pose-quality tool

PoseCheck is the preferred starting tool for PQI-related measurements because
it is designed to evaluate 3D ligand poses using measures such as ligand strain
and protein-ligand steric clashes.

PoseCheck will be used for pose-quality assessment, not for defining
state-specific residue interactions. Its outputs will be retained as separate
QC variables before any composite index is considered.

The exact PoseCheck version and the availability of each metric in the server
environment must be verified before implementation. The current MolScore
environment previously reported a missing PoseCheck-related module, so a
separate compatible environment may be required.

Reference: <https://posecheck.readthedocs.io/en/latest/>

### ProLIF: primary interaction-fingerprint tool

ProLIF is selected as the main quantitative residue-interaction tool because it
is designed to generate interaction fingerprints from docking poses and to
export pose-by-residue-by-interaction results directly as structured Python and
Pandas objects.

This suits the planned batch analysis of hundreds of poses and supports:

- residue-level Boolean or count fingerprints;
- multiple docking poses;
- interaction atom indices;
- distance and angle metadata for supported interactions;
- frequency and enrichment calculations;
- later extension to molecular-dynamics trajectories if required.

The primary interaction set should initially emphasize chemically specific
interactions:

```text
HBDonor
HBAcceptor
Anionic
Cationic
PiStacking
CationPi
PiCation
XBDonor
XBAcceptor
```

Hydrophobic interactions should be analysed separately, and broad van der
Waals contacts should not dominate the main enrichment table.

References:

- <https://prolif.readthedocs.io/en/latest/>
- <https://prolif.readthedocs.io/en/latest/notebooks/docking.html>

### PLIP: secondary cross-check and reporting tool

PLIP will not be the primary batch-analysis engine. It is a rule-based
protein-ligand interaction profiler that provides detailed residue, atom,
distance, angle, and interaction reports for individual complexes.

PLIP is well suited to:

- independently checking selected ProLIF interactions;
- inspecting representative high-PQI poses;
- producing detailed reports for key candidate molecules;
- supporting manual structural review and report figures.

PLIP is less convenient than ProLIF for the principal population-level
analysis because its per-complex reports require additional parsing and
standardization before they become a common residue-by-interaction matrix.
Differences between PLIP and ProLIF do not automatically imply an error because
their interaction definitions and geometric cutoffs are not identical.

References:

- <https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md>
- <https://pmc.ncbi.nlm.nih.gov/articles/PMC4489249/>

### Why ProLIF will not define pose quality

ProLIF determines whether specified interaction-geometry rules are satisfied;
it is not a global docking-pose quality predictor. Using the number of detected
ProLIF interactions as a pose-quality filter would also create circularity if
those same interactions were later tested for enrichment.

Accordingly, the main pose-quality filters will be based on independent clash,
strain, and geometry-validity measurements. PLIP or ProLIF geometry may be
reported as interaction-level supporting information but will not be the
principal PQI definition.

### Why docking score alone is insufficient

Glide score remains the selectivity metric used in the preceding analysis
stage, but it should not serve simultaneously as the sole pose-quality metric.
Using score both to define state selectivity and to filter pose credibility
would make the structural interpretation overly dependent on the same scoring
function.

## Required input preparation

Interaction profiling requires both ligand coordinates and the corresponding
protein structure. A Glide receptor grid alone is insufficient.

The workflow therefore requires:

- the prepared PPS receptor used to generate the PPS grid;
- the prepared PR receptor used to generate the PR grid;
- original AHC docked ligand poses;
- double-arm cross-docked ligand poses;
- stable molecule, variant, pose, and receptor-state identifiers.

Protein protonation, histidine state, explicit hydrogen handling, ligand formal
charge, and atom typing should remain consistent with the docking preparation
where possible.

PPS and PR receptor residues must also be mapped onto a common numbering
scheme. Each extracted interaction should preserve:

```text
receptor_state
original_chain
original_residue_number
residue_name
common_aligned_residue_id
```

Without this mapping, the same structural position may be incorrectly treated
as two different residues, or unrelated residue numbers may be compared
directly.

## Planned statistical outputs

The base long-format interaction table should contain one row per detected or
evaluated residue interaction:

```text
pose_id
crossdock_id
molecule_id
ligand_origin_state
receptor_state
pose_rank
ligand_variant
docking_score
selectivity_margin
selectivity_class
PQI metrics
common_residue_id
original_residue_id
interaction_type
interaction_present
interaction_count
distance
angle
```

Derived analyses will include:

- state-specific interaction frequency;
- PPS-versus-PR frequency difference;
- odds ratio and Fisher's exact test where appropriate;
- multiple-testing correction;
- bootstrap confidence intervals for effect sizes;
- interaction retention/loss between paired receptor states;
- association between interaction change and selectivity margin;
- PQI-threshold sensitivity curves;
- candidate-level interaction summaries for structural inspection.

Repeated variants or poses from the same molecule must not be treated as fully
independent biological observations. Primary population statistics should be
aggregated at molecule level, while pose-level support can be retained as an
uncertainty or consensus measure.

## Planned code architecture

The proposed next module is:

```text
s4_interaction_analysis/
├── README.md
├── input_models.py
├── receptor_mapping.py
├── pose_quality.py
├── interaction_fingerprint.py
├── interaction_statistics.py
├── run_raw_interactions.py
├── run_crossdock_interactions.py
└── summarize_sensitivity.py
```

Planned responsibilities:

- `input_models.py`
  - locate and validate receptor and ligand-pose files;
  - standardize pose, variant, molecule, and state identifiers;
  - assemble protein-ligand complexes where required.

- `receptor_mapping.py`
  - define PPS-to-PR residue correspondence;
  - preserve original and common residue identifiers.

- `pose_quality.py`
  - run PoseCheck-compatible QC calculations;
  - record hard failures and continuous QC metrics;
  - apply explicit threshold policies without deleting the raw measurements.

- `interaction_fingerprint.py`
  - generate ProLIF fingerprints;
  - preserve interaction geometry and atom indices;
  - export both long-format and matrix-format tables.

- `interaction_statistics.py`
  - aggregate pose results at molecule level;
  - calculate interaction frequency, enrichment, effect sizes, and corrected
    statistical tests;
  - calculate paired interaction retention/loss.

- `run_raw_interactions.py`
  - execute Stage A hypothesis-generation analysis on original AHC poses.

- `run_crossdock_interactions.py`
  - execute Stage B analysis on double-arm cross-docking poses;
  - test the pre-specified Stage A hypotheses.

- `summarize_sensitivity.py`
  - rerun or summarize results across PQI thresholds;
  - report retention and population-composition changes;
  - identify interaction trends stable across thresholds.

## Implementation sequence

Implementation should proceed incrementally:

1. verify one PPS and one PR prepared receptor structure;
2. confirm residue numbering and construct the common residue map;
3. load one Glide SDF pose with its matching receptor;
4. confirm PoseCheck metrics on that complex;
5. generate one ProLIF fingerprint and inspect it manually;
6. cross-check selected interactions using PLIP and Maestro/PyMOL;
7. define stable intermediate CSV schemas;
8. implement batch extraction;
9. implement molecule-level aggregation and enrichment statistics;
10. add PQI-threshold sensitivity analysis;
11. apply the locked workflow to double-arm cross-docking results.

No population-scale conclusion should be drawn until the single-complex test
confirms correct protonation, residue identification, coordinate alignment, and
interaction detection.

## Interpretation limits

The planned analysis can identify docking-derived structural patterns
consistent with state selectivity. It cannot by itself prove the experimental
binding pose, binding affinity, or causal role of an individual residue.

The strongest resulting statement will therefore be of the form:

> A residue-level interaction pattern was enriched among docking-predicted
> state-selective ligands, remained stable across pose-quality thresholds, and
> was preferentially retained in the corresponding receptor state during
> paired cross-docking.

Experimental structures, mutagenesis, biochemical measurements, or
higher-level simulation would be required for stronger mechanistic claims.
