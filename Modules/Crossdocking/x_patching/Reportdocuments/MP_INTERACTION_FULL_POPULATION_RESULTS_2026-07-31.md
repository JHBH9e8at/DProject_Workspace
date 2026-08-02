# MP interaction full-population analysis notes

Date: 2026-07-31  
Project: Myosin PPS/PR state-directed molecular generation  
Status: Full-population PoseCheck and ProLIF extraction completed; initial comparative interpretation completed

## Purpose

This note records the first full-population structural-interaction analysis of
the independently generated PPS- and PR-directed AHC populations. The analysis
was designed to:

1. verify that the multiprocessing interaction run completed without missing
   or duplicated population records;
2. summarize PoseCheck clash and ligand-strain measurements;
3. identify the dominant ProLIF residue-interaction patterns in each
   receptor-state population;
4. determine whether the PPS and PR prepared receptors already use a common
   residue numbering system;
5. define the quality-control and enrichment analyses required before
   interaction patterns are interpreted as state-associated hypotheses.

The results are exploratory population-level docking analyses. They do not
establish experimental binding poses, affinities, or causal residue effects.

## Input populations

The interaction run used the same deduplicated full populations analysed in
the preceding PR/PPS optimization and diversity report.

| Population | Analysed poses | Unique molecules | Valid pose records |
|---|---:|---:|---:|
| PPS | 20,420 | 20,420 | 20,420 |
| PR | 19,688 | 19,688 | 19,688 |
| Total | 40,108 | 40,108 | 40,108 |

Primary result directory:

```text
Q:\coding_dir\701_Project\Resultsbin\MP_Interrun
```

Primary prepared receptors:

```text
Q:\coding_dir\701_Project\workspace\AHC_Related\ref_structures\Structures\
    PPS_Protein_no_conect.pdb
    PR_Protein_no_conect.pdb
```

## Multiprocessing execution

The full run used:

```text
PoseCheck workers:       30
PoseCheck chunk size:    250 poses
ProLIF workers:          30
Selection mode:          full population
Analysis mode:           PoseCheck and ProLIF
```

The run started at 2026-07-30 18:25:54 UTC and finished at
2026-07-30 22:57:37 UTC, for a total wall time of approximately 4 hours
32 minutes.

Component wall times:

| Population | Component | Output rows | Wall time |
|---|---|---:|---:|
| PPS | PoseCheck | 20,420 | 2 h 51 min |
| PPS | ProLIF | 109,622 | 1 min 31 s |
| PR | PoseCheck | 19,688 | 1 h 37 min |
| PR | ProLIF | 81,794 | 1 min 10 s |

PoseCheck was the dominant computational cost. ProLIF contributed only a
small fraction of the total wall time.

## Completion and data-integrity audit

The top-level run status was `completed`. All 40,108 manifest rows had
`interaction_status=completed`.

The following checks passed:

- missing source pose files: 0;
- non-empty interaction errors: 0;
- duplicate population-state/molecule identifiers: 0;
- duplicate pose indices within either population: 0;
- invalid pose-read records: 0;
- incomplete PoseCheck records: 0;
- unparsed ProLIF protein residue identifiers: 0.

The aggregate row counts reconciled with the PPS and PR population-level
outputs.

## Pose-quality measurements

### Docking and molecular-size context

| Statistic | PPS | PR |
|---|---:|---:|
| Mean docking score | -10.782 | -8.140 |
| Median docking score | -10.959 | -8.204 |
| Mean heavy-atom count | 35.45 | 28.02 |
| Median heavy-atom count | 36 | 28 |
| Mean rotatable-bond count | 8.37 | 8.13 |
| Median formal charge | +1 | +1 |

The more favourable PPS docking-score distribution is accompanied by a
substantially larger ligand population. Absolute clash and strain
measurements must therefore be interpreted together with size-normalized
measurements.

### Clash distributions

| Statistic | PPS clash count | PR clash count | PPS clash/heavy atom | PR clash/heavy atom |
|---|---:|---:|---:|---:|
| Mean | 7.42 | 5.87 | 0.208 | 0.210 |
| Median | 7 | 6 | 0.195 | 0.200 |
| 5th percentile | 2 | 2 | 0.063 | 0.059 |
| 95th percentile | 14 | 11 | 0.387 | 0.400 |
| Maximum | 35 | 21 | 1.111 | 0.731 |

PPS had a higher absolute clash count, but the heavy-atom-normalized
distributions were very similar. Much of the absolute difference is
therefore consistent with the larger PPS ligand size.

### Ligand-strain distributions

| Statistic | PPS strain energy | PR strain energy | PPS strain/rotatable bond | PR strain/rotatable bond |
|---|---:|---:|---:|---:|
| Mean | 34.97 | 25.86 | 4.27 | 3.36 |
| Median | 30.72 | 20.67 | 3.66 | 2.58 |
| 5th percentile | 8.81 | 5.95 | 1.37 | 0.96 |
| 95th percentile | 75.49 | 60.63 | 9.16 | 7.90 |
| Maximum | 500.03 | 724.57 | 117.62 | 159.93 |

`strain_per_rotatable_bond` was undefined for 15 PPS and 34 PR ligands with
zero rotatable bonds. Raw strain energy remained available for those records.

PPS showed higher central strain than PR even after normalization by
rotatable-bond count. Both populations contained extreme upper-tail strain
outliers.

### Relationship with docking score

Spearman correlations with docking score:

| PoseCheck measurement | PPS | PR |
|---|---:|---:|
| Clash count | -0.364 | -0.510 |
| Clashes per heavy atom | -0.149 | -0.406 |
| Strain energy | -0.336 | -0.237 |
| Strain per rotatable bond | -0.124 | +0.030 |

Because a more negative docking score is treated as more favourable, a
negative coefficient means that better-scoring poses tended to have more
clashes or higher strain. This effect was particularly clear for PR clash
measurements. Docking score alone should therefore not be used as a structural
plausibility filter.

### Current pose-quality pass flag

All 40,108 poses currently have `pose_quality_pass=True`. This must not be
interpreted as evidence that all poses are structurally acceptable.

The multiprocessing configuration did not specify:

```text
max_clashes
max_clashes_per_heavy_atom
max_strain_energy
```

Consequently, the pass column reflects the absence of active thresholds
rather than a discriminating quality-control decision. Threshold selection
and sensitivity analysis remain required.

## ProLIF interaction coverage

| Statistic | PPS | PR |
|---|---:|---:|
| Interaction occurrences | 109,622 | 81,794 |
| Poses with at least one interaction | 19,601 | 19,438 |
| Pose coverage | 96.0% | 98.7% |
| Poses without a detected interaction | 819 | 250 |
| Mean occurrences per all poses | 5.37 | 4.15 |
| Mean occurrences per interacting pose | 5.59 | 4.21 |

PPS produced more interaction occurrences per pose. This may partly reflect
its larger average ligand size and should be re-evaluated after size and
pose-quality stratification.

## Interaction-type distributions

| Interaction type | PPS occurrences | PR occurrences |
|---|---:|---:|
| HBDonor | 42,519 | 44,157 |
| HBAcceptor | 39,442 | 27,793 |
| Cationic | 11,012 | 7,133 |
| PiStacking | 8,749 | 349 |
| PiCation | 5,801 | 1,653 |
| CationPi | 1,700 | 0 |
| Anionic | 311 | 635 |
| XBDonor | 88 | 74 |

The strongest unadjusted population contrast was the much higher PPS
PiStacking count. This is a candidate state-associated signal, but it may also
reflect ligand aromaticity, size, binding-site occupancy, or pose-quality
composition. It should not be interpreted as enrichment until molecule-level
frequencies and covariate sensitivity are calculated.

## Dominant PPS residue-interaction patterns

Pose prevalence was calculated as the number of PPS poses containing a given
residue-interaction type divided by all 20,420 PPS poses.

| Residue | Interaction | Poses | Prevalence |
|---|---|---:|---:|
| ASP168.A | HBDonor | 14,787 | 72.4% |
| ASN711.A | HBAcceptor | 14,011 | 68.6% |
| ARG712.A | HBAcceptor | 13,689 | 67.0% |
| HIS666.A | PiStacking | 7,394 | 36.2% |
| GLY768.A | HBAcceptor | 5,292 | 25.9% |
| ARG712.A | PiCation | 4,044 | 19.8% |
| ASP168.A | Cationic | 2,648 | 13.0% |
| ASP89.A | HBDonor | 2,552 | 12.5% |
| ARG712.A | HBDonor | 2,265 | 11.1% |
| GLU170.A | HBDonor | 2,253 | 11.0% |

The unfiltered PPS pattern is dominated by the ASP168 hydrogen-donor contact,
ASN711/ARG712 acceptor contacts, and HIS666 aromatic stacking.

## Dominant PR residue-interaction patterns

| Residue | Interaction | Poses | Prevalence |
|---|---|---:|---:|
| ALA91.A | HBDonor | 8,195 | 41.6% |
| SER118.A | HBAcceptor | 7,519 | 38.2% |
| ARG712.A | HBAcceptor | 6,890 | 35.0% |
| PRO710.A | HBDonor | 6,543 | 33.2% |
| ASN711.A | HBAcceptor | 5,605 | 28.5% |
| GLU497.A | HBDonor | 5,020 | 25.5% |
| GLU497.A | Cationic | 4,436 | 22.5% |
| ARG712.A | HBDonor | 4,207 | 21.4% |
| SER118.A | HBDonor | 3,409 | 17.3% |
| TYR722.A | HBAcceptor | 3,323 | 16.9% |

The unfiltered PR pattern is more distributed across ALA91, SER118, PRO710,
GLU497, ARG712, and ASN711.

## Preliminary cross-state interpretation

ARG712 and ASN711 were prominent in both populations, although their
prevalences differed substantially. Candidate state-associated features
include:

- PPS: ASP168 HBDonor/Cationic, HIS666 PiStacking, GLY768 HBAcceptor, and
  elevated PiStacking generally;
- PR: ALA91 HBDonor, SER118 HBAcceptor/HBDonor, PRO710 HBDonor, and GLU497
  HBDonor/Cationic.

These are hypothesis-generating observations. Direct frequency differences
are potentially confounded by:

- independent PPS and PR ligand populations;
- different molecular-size and aromaticity distributions;
- different docking-score distributions;
- pose clash and strain;
- chemical-series and scaffold dependence;
- receptor-state pocket geometry.

The next analysis should calculate molecule-level prevalence, effect sizes,
confidence intervals, and quality-threshold sensitivity rather than compare
raw occurrence totals alone.

## Residue-numbering and mapping audit

No verified mapping CSV was supplied to this run. Accordingly:

```text
common_residue_id:       empty
residue_mapping_status: not_requested
```

However, the two prepared receptor PDB files were audited directly using
fixed-column PDB residue and atom identifiers.

| Audit item | PPS | PR |
|---|---:|---:|
| Protein chain | A | A |
| Residue range | 1–783 | 1–783 |
| ATOM residues | 783 | 783 |
| Parsed named atoms | 12,524 | 12,522 |

Cross-receptor comparison:

```text
Common chain/residue/insertion-code keys: 783
Matching residue names at common keys:    783
Residue-name mismatches:                  0
PPS-only residue keys:                    0
PR-only residue keys:                     0
Common CA atoms:                          783
```

The residue sequence, chain identifier, residue numbering, and insertion-code
scheme are therefore identical across PPS and PR. An identity mapping is
appropriate:

```text
PPS A:<residue number>:<residue name>
    -> common A:<residue number>:<residue name>
PR  A:<residue number>:<residue name>
    -> common A:<residue number>:<residue name>
```

The receptors are not coordinate-identical. The fitted CA RMSD across all 783
residues was approximately 9.05 Å, consistent with substantially different
receptor conformations despite identical residue identity and numbering.
This is expected for the PPS and PR structural states and does not invalidate
the identity residue mapping.

For reproducibility, a formal identity-mapping CSV should still be generated
and supplied in future runs. The present interaction table can also be
annotated post hoc because every parsed residue key has an unambiguous common
identity.

## Recommended pose-quality sensitivity analysis

Thresholds should not be selected solely to maximize the apparent
state-specific interaction contrast. A staged sensitivity analysis is
recommended.

Possible empirical starting subsets include:

```text
All valid poses
Moderate QC: clash/heavy atom <= state-specific 75th percentile
Strict QC:   clash/heavy atom <= state-specific median
Additional strain filters at the 75th, 50th, and 25th percentiles
```

Direct physical cutoffs can be added after reviewing PoseCheck definitions and
representative structures. At each threshold, record:

- retained pose count and fraction by state;
- heavy-atom, aromaticity, formal-charge, and rotatable-bond distributions;
- docking-score distributions;
- residue-interaction prevalence and effect-size changes;
- scaffold representation;
- whether the direction of each candidate interaction contrast is preserved.

## Recommended interaction statistics

The next primary table should have one row per molecule and one Boolean or
count column per common-residue/interaction pair. Recommended outputs are:

- PPS and PR molecule-level prevalence;
- absolute prevalence difference;
- odds ratio with confidence interval;
- Fisher exact or an appropriate regression-based test;
- false-discovery-rate correction;
- bootstrap confidence intervals;
- sensitivity to pose-quality filters;
- sensitivity to molecular size, aromatic-ring count, charge, and scaffold.

Because the original PPS and PR populations contain different molecules, the
initial comparison is unpaired. Paired retention/loss statistics should be
reserved for the later double-arm cross-docking dataset.

## Suggested report-ready Methods text

> Full-population docked poses from the PPS- and PR-directed AHC runs were
> analysed using PoseCheck and ProLIF. PoseCheck clash count and ligand strain
> energy were retained as separate measurements together with normalization
> by ligand heavy-atom and rotatable-bond counts. ProLIF interactions were
> exported at occurrence level with receptor state, molecule, pose, residue,
> interaction type, atom indices, and available geometry. The computation used
> 30 process workers for chunked PoseCheck evaluation and 30 ProLIF workers.
> All valid population records were retained for the initial exploratory
> summary.

## Suggested report-ready Results text

> The analysis completed for 20,420 PPS and 19,688 PR poses without missing
> source structures, duplicated pose indices, or interaction-extraction
> failures. Median clash counts were 7 for PPS and 6 for PR, whereas median
> clashes per heavy atom were similar at 0.195 and 0.200, respectively. Median
> ligand strain energy was higher for PPS than PR (30.72 versus 20.67).
> ProLIF identified at least one interaction for 96.0% of PPS poses and 98.7%
> of PR poses. The most frequent PPS residue-interaction patterns included
> ASP168 HBDonor, ASN711 HBAcceptor, ARG712 HBAcceptor, and HIS666 PiStacking.
> PR was characterized by ALA91 HBDonor, SER118 HBAcceptor, ARG712 HBAcceptor,
> PRO710 HBDonor, and GLU497 HBDonor/Cationic. Residue-level comparison of the
> two prepared receptors confirmed identical chain-A residue identity and
> numbering across all 783 positions, permitting a direct identity mapping
> between receptor states.

## Interpretation limits

The PoseCheck measurements are computational plausibility indicators, not
probabilities that the predicted poses are correct. No calibrated or
pre-specified pose-quality thresholds were applied in the current run.

ProLIF interaction calls depend on geometric definitions and the prepared
protein and ligand protonation states. Interaction occurrence is not proof of
an energetic or causal contribution to binding.

The PPS and PR populations differ in molecular size and chemical composition.
Raw interaction counts and prevalences therefore cannot be attributed solely
to receptor conformation.

The residue-numbering audit supports identity mapping between these specific
prepared PDB files. It does not imply that residue coordinates are equivalent;
the receptor conformations differ substantially.

## Reproducibility record

Primary outputs:

```text
Q:\coding_dir\701_Project\Resultsbin\MP_Interrun\
    interaction_run_summary.json
    manifests\interaction_manifest.csv
    aggregate\pose_metadata.csv
    aggregate\pose_quality.csv
    aggregate\interactions_long.csv
    aggregate\run_summary.csv
```

Implementation:

```text
Q:\coding_dir\701_Project\workspace\Modules\Crossdocking\analysis\MP_interaction
```

## Current conclusion

The multiprocessing workflow successfully produced a complete
full-population structural-interaction dataset. The initial results show
distinct PPS and PR residue-interaction patterns, with a particularly strong
unadjusted PPS PiStacking signal and different dominant contact regions in the
two receptor states.

The two prepared receptors already share an exact residue identity and
numbering system across chain A residues 1–783. A formal identity mapping can
therefore be generated without sequence alignment and applied post hoc to the
current interaction table.

Before state-associated interaction hypotheses are promoted, the analysis
must apply explicit pose-quality sensitivity filters and molecule-level
effect-size statistics while controlling for ligand size and chemical-series
composition.
