# Stage 4: pose quality and residue interactions

This module is the planned structural-interpretation layer after the
score-based analyses in `analysis.score`.

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

The primary entry point now runs PPS and PR original-docking interaction
analysis in batch from one flat `.in` file. It accepts the two AHC `scores.csv`
files and resolves each selected best variant using:

```text
<run_dir>/<run_name>_GlideDock/<step>/<best_variant>_lib.sdfgz
```

The existing single-complex function remains the worker used for each resolved
pose file.

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

Population-scale interaction enrichment and cross-docked-pose interaction
comparison remain separate future stages. This runner performs auditable batch
extraction over the original PPS and PR AHC docking results.

## Files

```text
analysis/interaction/
├── input_models.py
├── receptor_mapping.py
├── pose_quality.py
├── fingerprint.py
├── single_complex.py
├── batch_runner.py
├── interaction_runner.py
├── test_batch_runner.py
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

### PPS/PR receptor PDB requirement

For this project, use the prepared hydrogen-containing PPS and PR receptor PDBs
with all `CONECT` records removed:

```text
PPS_Protein_no_conect.pdb
PR_Protein_no_conect.pdb
```

Do not pass the original Maestro-exported PDBs containing only a partial set of
`CONECT` records to ProLIF. MDAnalysis may treat that partial topology as the
available bond topology, which prevents ProLIF from recognizing protein
donor/acceptor patterns correctly. In the validated five-pose test, the partial
`CONECT` receptors produced zero interactions, whereas the corresponding
no-`CONECT` receptors produced 44 PPS and 39 PR primary interaction
occurrences.

The no-`CONECT` files must otherwise preserve the original prepared
coordinates, atom records, explicit hydrogens, protonation state, and residue
identifiers. They can be produced without modifying the source files:

```bash
grep -v '^CONECT' PPS_Protein.pdb > PPS_Protein_no_conect.pdb
grep -v '^CONECT' PR_Protein.pdb > PR_Protein_no_conect.pdb
```

With these files, the following MDAnalysis warning is expected during ProLIF
conversion and indicates the intended bond-inference path:

```text
No `bonds` attribute in this AtomGroup. Guessing bonds based on atoms coordinates
```

PoseCheck also accepts these no-`CONECT` PDBs because it uses the preserved
prepared coordinates and explicit hydrogens.

## Batch run from `.in`

Start from the repository root:

```bash
python -m analysis.interaction.interaction_runner \
  --config interaction.in
```

Required configuration:

```ini
pps_results=/path/to/PPS/scores.csv
pr_results=/path/to/PR/scores.csv
pps_run_dir=/path/to/PPS/AHC_run
pr_run_dir=/path/to/PR/AHC_run
pps_receptor=/path/to/PPS_Protein_no_conect.pdb
pr_receptor=/path/to/PR_Protein_no_conect.pdb
output_dir=/path/to/interaction_results
```

These format keys have built-in defaults and normally do not need to appear:

```ini
pps_run_name=PPS
pr_run_name=PR
step_col=step
pps_variant_col=PPS_best_variant
pr_variant_col=PR_best_variant
smiles_col=smiles
```

`testmode_top_n=` selects the complete usable population after validity,
uniqueness, score, and variant checks. A positive value selects the best N
own-state docking scores independently from PPS and PR.

Prepared receptor coordinates are mandatory. A `run_dir` resolves ligand poses;
it does not replace the receptor required by ProLIF or PoseCheck.

The previous one-file workflow remains available for validation:

```bash
python -m analysis.interaction.single_complex \
  --protein /path/to/prepared_PPS_receptor.pdb \
  --poses /path/to/docked_poses.sdfgz \
  --receptor-state PPS \
  --output-dir /path/to/single_test \
  --analysis both
```

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

Set `include_secondary_interactions=on` in the `.in` file to add `Hydrophobic`
and `VdWContact`.

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

Each state is materialized and analyzed as one population batch under
`populations/PPS/` and `populations/PR/`, avoiding repeated receptor loading for
every molecule. Consolidated tables are written under `aggregate/`:

```text
aggregate/pose_metadata.csv
aggregate/pose_quality.csv
aggregate/interactions_long.csv
aggregate/run_summary.csv
manifests/interaction_manifest.csv
interaction_run_summary.json
```

`aggregate/interactions_long.csv` contains detected interaction occurrences
only. `aggregate/pose_metadata.csv` remains the authoritative table for all poses,
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
