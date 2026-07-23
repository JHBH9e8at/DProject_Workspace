# Stage 4: pose quality and residue interactions

This module is the planned structural-interpretation layer after the
score-based analyses in `s3_analysis`.

It separates two measurements:

```text
PoseCheck -> pose-quality metrics (clashes and ligand strain)
ProLIF    -> residue-level interaction fingerprints
```

The first implementation deliberately does not convert the raw PoseCheck
metrics into one weighted confidence probability. Raw measurements are retained
and explicit thresholds are recorded so that sensitivity to the QC definition
can be evaluated later.

## Current implementation scope

The current code supports one prepared receptor and one Glide pose file:

```text
prepared receptor (.pdb or .mol2)
docked poses (.sdf or .sdfgz)
```

ProLIF supports the documented receptor formats above. The current PoseCheck
adapter requires its receptor input as `.pdb`; therefore `--analysis both` also
requires a PDB receptor.

It provides:

- validated pose metadata and stable zero-based pose indices;
- transparent PoseCheck clash and strain metrics;
- optional pose-quality thresholds;
- ProLIF occurrence-level interaction fingerprints;
- optional PPS/PR common residue mapping;
- separate primary and broad secondary interaction sets.

Population-scale statistics and double-arm batch runners are intentionally
deferred until a real PPS and PR receptor/pose pair passes this single-complex
test.

## Files

```text
s4_interaction_analysis/
├── input_models.py
├── receptor_mapping.py
├── pose_quality.py
├── interaction_fingerprint.py
├── run_single_complex.py
├── residue_mapping_template.csv
└── README.md
```

## Dependencies

The base schema and mapping code requires:

```text
pandas
RDKit
```

Pose-quality analysis additionally requires PoseCheck. Interaction analysis
requires ProLIF and MDAnalysis. These are optional imports: the package can be
inspected and its non-tool-dependent functions tested before the server
analysis environment is prepared.

The receptor must contain the hydrogens and preparation state intended for the
analysis. A Glide grid ZIP is not a replacement for the receptor coordinate
file.

## Single-complex run

Run both analysis components:

```bash
python s4_interaction_analysis/run_single_complex.py \
  --protein /path/to/prepared_PPS_receptor.pdb \
  --poses /path/to/docked_poses.sdfgz \
  --receptor-state PPS \
  --output-dir /path/to/interaction_test \
  --analysis both
```

Run ProLIF only:

```bash
python s4_interaction_analysis/run_single_complex.py \
  --protein /path/to/prepared_PR_receptor.pdb \
  --poses /path/to/docked_poses.sdf \
  --receptor-state PR \
  --output-dir /path/to/prolif_test \
  --analysis prolif
```

Apply explicit PoseCheck thresholds:

```bash
python s4_interaction_analysis/run_single_complex.py \
  --protein /path/to/prepared_PPS_receptor.pdb \
  --poses /path/to/docked_poses.sdf \
  --receptor-state PPS \
  --output-dir /path/to/posecheck_test \
  --analysis posecheck \
  --max-clashes 2 \
  --max-clashes-per-heavy-atom 0.10 \
  --max-strain-energy 20
```

Threshold values above are usage examples, not validated project cutoffs.

## Primary ProLIF interactions

The main fingerprint excludes broad hydrophobic and van der Waals contacts by
default:

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

Use `--include-secondary-interactions` to add `Hydrophobic` and `VdWContact`.

## Residue mapping

`residue_mapping_template.csv` demonstrates the required schema:

```text
receptor_state
original_chain
original_residue_number
residue_name
common_residue_id
```

The example rows are placeholders and must not be treated as the real myosin
mapping. After constructing a verified mapping, pass it using:

```bash
--residue-map /path/to/verified_residue_mapping.csv
```

## Outputs

Depending on the selected analysis:

```text
pose_metadata.csv
pose_quality.csv
interactions_long.csv
run_summary.csv
intermediate/*.sdf
```

`interactions_long.csv` contains detected interaction occurrences only.
`pose_metadata.csv` remains the authoritative table for all input poses,
including poses with no detected primary interaction.

## Required validation before batch implementation

1. Confirm that the receptor is the prepared structure used for the matching
   Glide grid.
2. Confirm explicit hydrogen and protonation handling.
3. Run one PPS and one PR example.
4. Compare selected interactions manually in Maestro or PyMOL.
5. Cross-check representative interactions with PLIP.
6. Verify the PPS/PR common residue mapping.
7. Freeze the intermediate CSV schemas before adding population runners.
