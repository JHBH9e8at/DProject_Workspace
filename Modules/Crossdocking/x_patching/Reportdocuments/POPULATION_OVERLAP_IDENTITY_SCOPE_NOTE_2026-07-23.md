# Population-overlap molecular identity: problem definition and project scope

Date: 2026-07-23  
Project: Myosin PPS/PR state-selectivity analysis

## Problem definition

Before applying the final filtering workflow, the raw PPS and PR AHC
populations can be searched for molecules occurring in both `scores.csv`
files. Because the two populations were evaluated with equivalent LigPrep and
Glide protocols but against different myosin-state receptor grids, an
overlapping molecule provides a naturally paired PPS-versus-PR docking-score
comparison.

The main methodological problem is defining what counts as the same molecule.
Raw SMILES strings cannot be used directly because the same chemical structure
can be written with different atom orderings. Conversely, overly broad
standardization could incorrectly merge stereoisomers, protonation states, or
other chemically distinct states that may bind differently.

Three-dimensional conformational differences, including different bond,
torsion, and dihedral angles, should not create different molecular identities.
These conformational differences are not encoded by ordinary SMILES and are
therefore naturally ignored during SMILES-based identity matching.

## Solution within the current project scope

The primary overlap analysis will use **canonical isomeric SMILES** as the
molecular identity key.

Each input SMILES will be parsed and rewritten in canonical isomeric form
before matching the PPS and PR populations. This approach:

- ignores the original SMILES atom ordering and writing direction;
- ignores three-dimensional conformation and docking pose;
- preserves atom connectivity and bond order;
- preserves formal charge;
- preserves tetrahedral and double-bond stereochemistry;
- does not neutralize charged structures;
- does not collapse different tautomeric representations.

Accordingly, differently written SMILES that describe the same defined
chemical structure will match, while neutral and formally charged states will
remain separate. Enantiomers and other defined stereoisomers will also remain
separate because they can have different receptor-binding behaviour.

Exact raw-SMILES matching will not be used because it adds representation-level
strictness without improving the chemical interpretation of the selectivity
analysis.

## Planned overlap-score analysis

For each canonical isomeric SMILES present in both populations, repeated
instances will first be aggregated independently within the PPS and PR runs.
The analysis should retain occurrence count, best score, median score, mean
score, and score variability for each state.

The primary paired metric will be:

```text
state_score_difference = PR_score - PPS_score
```

Because more negative docking scores are more favourable:

- a positive difference indicates PPS preference;
- a negative difference indicates PR preference;
- a difference close to zero indicates weak or inconclusive state preference.

Best-score and median-score differences should both be retained. Agreement
between their signs will provide a simple robustness indicator, while repeated
instances will be aggregated before statistical analysis to avoid treating
duplicate generations as independent observations.

## Analyses outside the current implementation scope

Two broader identity analyses are scientifically interesting but will not be
implemented in the initial population-overlap module.

### Tautomer-normalized matching

Canonical tautomer matching could identify input SMILES that represent the
same parent compound through different tautomeric forms. This may be useful
because Epik/LigPrep can generate related tautomer ensembles downstream.
However, it changes the identity definition and is better treated as a
secondary sensitivity or structure-based analysis.

### Canonical non-isomeric matching

Non-isomeric canonical SMILES would merge structures that differ only in
specified stereochemistry, including enantiomeric or diastereomeric cases.
Comparing their docking behaviour could reveal stereochemistry-dependent
binding, but those structures should not be pooled in the primary
state-selectivity analysis. They may instead be studied later as related
structure pairs.

## Scope decision

The initial implementation will therefore:

1. read the unfiltered PPS and PR `scores.csv` populations;
2. remove only technically unusable records, such as unparsable SMILES and
   missing, zero, or non-negative docking scores;
3. construct a canonical isomeric SMILES key without charge neutralization or
   tautomer normalization;
4. aggregate repeated instances within each population;
5. identify canonical molecules shared by PPS and PR;
6. calculate paired PPS/PR score differences and robustness summaries.

Tautomer-family and non-isomeric stereochemical analyses are retained as
planned future structure-based analyses and are intentionally separated from
the primary selectivity workflow.
