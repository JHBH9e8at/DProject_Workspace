# LigPrep protocol consistency and validation notes

Date: 2026-07-23  
Project: Myosin PPS/PR state cross-docking  
Purpose: Record the evidence required to justify direct comparison between the
original AHC own-state docking scores and newly generated opposite-state
cross-docking scores.

## 1. Experimental context

Two independent SMILES-RNN + MolScore AHC populations are used:

```text
PPS-generated population -> PPS score from original AHC run
                         -> cross-docking to PR receptor grid

PR-generated population  -> PR score from original AHC run
                         -> cross-docking to PPS receptor grid
```

For a meaningful paired state-selectivity comparison, ligand preparation and
Glide docking conditions must be consistent between the original AHC scoring
and the cross-docking calculation. Glide settings were separately confirmed to
be equivalent. This note documents the LigPrep consistency check.

## 2. Original AHC configuration

The PPS AHC run configuration contains the following GlideDock preparation
setting:

```json
{
  "name": "GlideDock",
  "parameters": {
    "prefix": "PPS",
    "glide_template": "/home/andy/proj/701/input/PPSGlide.in",
    "cluster": "tcp://138.37.52.153:8786",
    "timeout": 600.0,
    "ligand_preparation": "ligprep"
  }
}
```

Relevant local copy:

```text
Q:\coding_dir\Workspace1\Docking\2026_06_02_SMILES-RNN_PPS_tester_run_2\PPS_tester_run_2_config.json
```

The AHC JSON requests MolScore's `ligprep` preparation route but does not
override its individual LigPrep parameters. Therefore the defaults in the
installed MolScore implementation were used.

## 3. MolScore implementation used by the AHC run

Server environment:

```text
Conda environment: molscore
Python: 3.9
MolScore package:
/home/andy/Project/conda_envs/molscore/lib/python3.9/site-packages/molscore
```

Relevant implementation file:

```text
/home/andy/Project/conda_envs/molscore/lib/python3.9/site-packages/molscore/scoring_functions/_ligand_preparation.py
```

The installed preparation class defines:

```python
pH: float = 7.0
pHt: float = 1.0
bff: int = 16
max_stereo: int = 8
use_epik: bool = True
host: Optional[str] = None
nojobid: bool = True
```

It constructs the following command:

```python
cmd = [
    ligprep_executable,
    "-ismi", smi_input,
    "-osd", sdf_output,
    "-ph", "7.0",
    "-pht", "1.0",
    "-bff", "16",
    "-s", "8",
    "-WAIT",
    "-epik",
    "-NOJOBID",
]
```

The installed file contains a local modification changing `fail_fast` from
`True` to `False`. This changes error handling only: a failed molecule can be
skipped and assigned a fallback score rather than terminating the full AHC
batch. It does not change the successful molecule preparation parameters
listed above.

## 4. Meaning of the LigPrep arguments

| Argument | Value | Meaning in this workflow |
|---|---:|---|
| `-ph` | `7.0` | Target/effective pH for ionization-state generation |
| `-pht` | `1.0` | pH tolerance around the target pH |
| `-bff` | `16` | S-OPLS for final LigPrep geometry optimization |
| `-s` | `8` | Maximum stereoisomer generation setting used by MolScore |
| `-epik` | enabled | Use Epik Classic for ionization and tautomer generation |
| `-WAIT` | enabled | Wait for LigPrep completion before continuing |
| `-NOJOBID` | enabled | Suppress standard Job Control ID output behavior |

`-WAIT` and `-NOJOBID` are execution-control options and do not change the
prepared chemical structures.

The total number of final SDF records can exceed eight because the stereoisomer
setting is not a global cap on every combination of stereochemistry,
ionization, and tautomeric state generated through the complete workflow.

## 5. Server-side confirmation of the force-field mapping

The exact Schrodinger installation used on the server is:

```text
SCHRODINGER=/home/apps/schrodinger2024-2
```

The installed LigPrep help was queried with:

```bash
"$SCHRODINGER/ligprep" -long_help \
  | grep -A 5 -B 2 -- '-bff'
```

The server returned:

```text
Force-field based geometry optimization:
  -bff {14,16}          Force-field to be used for the final geometry
                        optimization. Default: 14 (OPLS_2005). For S-OPLS
                        specify 16.
```

This establishes the mapping for the exact server installation:

```text
-bff 14 -> OPLS_2005
-bff 16 -> S-OPLS
```

Therefore the MolScore default `bff=16` explicitly requested S-OPLS.

## 6. Empirical validation from actual AHC output structures

The AHC configuration and source code were not used as the sole evidence. SDF
properties embedded in actual AHC `*_lib.sdfgz` output files were audited.

Thirty readable output files were sampled from distinct PPS AHC iteration
directories. All 30 contained:

```text
i_lp_mmshare_version = 66139
s_lp_Force_Field     = S-OPLS
s_epik_cmdline       = present
```

Results:

```text
Files audited:                 30
S-OPLS records:               30 / 30
Records containing Epik cmd:  30 / 30
mmshare version 66139:        30 / 30
Exceptions observed:           0
```

The Base64-encoded `s_epik_cmdline` property was decoded. Every sampled file
contained the same internal Epik call:

```text
'epik_python',
'-pht', '1.0',
'-ph', '7.0',
'-tn', '8',
'-ma', '200',
'-imae', '<infile.mae>',
'-omae', '<outfile.mae>'
```

This confirms that the requested settings were not merely present in source
code; they were used to generate the actual AHC ligand structures.

## 7. Cross-docking LigPrep implementation

Cross-docking preparation module:

```text
Crossdocking/s1_Prep/prepare_ligand.py
```

The cross-docking command was updated to reproduce the MolScore AHC settings:

```bash
ligprep \
  -ismi <input.smi> \
  -osd <prepared.sdf> \
  -ph 7.0 \
  -pht 1.0 \
  -bff 16 \
  -s 8 \
  -WAIT \
  -epik \
  -NOJOBID
```

For reliable file placement, the implementation also runs LigPrep with the
job-specific LigPrep output directory as both the process working directory
and `PWD`, while using an explicit output path. This execution-path handling
does not alter the preparation chemistry.

## 8. Empirical validation of the updated cross-docking output

A five-molecule PPS-to-PR test was completed after matching the AHC protocol.

Quality-control outcome:

```text
Selected instances:       5
Completed instances:      5
Failed instances:         0
Unavailable instances:    0
Prepared variants:       34
Docked poses:            34
Best-score rows:          5
NaN/zero/positive score:  0
```

The new cross-docking SDF-derived table contained:

```text
i_lp_mmshare_version = 66139
s_lp_Force_Field     = S-OPLS
s_epik_cmdline       = present
```

Its decoded Epik command was identical to the original AHC output:

```text
-pht 1.0
-ph 7.0
-tn 8
-ma 200
```

The cross-docking output additionally contained expected Epik properties such
as total charge and state-penalty fields. These properties had been absent from
the earlier unmatched default-LigPrep test.

Relevant downloaded test files:

```text
/home/andy/proj/Results/temptesting/PPStoPRtest_2/best_docking_scores.csv
/home/andy/proj/Results/temptesting/PPStoPRtest_2/combined_docking_scores2.csv
/home/andy/proj/Results/temptesting/PPStoPRtest_2/best_docking_scores2.csv
```

## 9. Why matching the protocol mattered

The first pipeline test used a minimal LigPrep command without explicit Epik,
pH, force-field, or stereoisomer options. Its output recorded OPLS_2005 and
generated 16 variants across the five molecules later retested.

The AHC-matched protocol recorded S-OPLS and generated 34 variants across the
same five molecules. Individual best PR-grid scores changed by approximately
`-3.12` to `+1.12` Glide score units relative to the unmatched test.

Examples:

| Molecule | Unmatched variants | Matched variants | Unmatched PR score | Matched PR score |
|---|---:|---:|---:|---:|
| `244_14-6` | 2 | 12 | -8.56047 | -9.67572 |
| `237_40-1` | 1 | 4 | -5.85758 | -8.97312 |
| `301_44-1` | 4 | 8 | -10.24850 | -9.12449 |
| `354_49-1` | 1 | 2 | -8.54508 | -8.13662 |

This demonstrates that ligand preparation was a material experimental
variable and validates the decision to align the cross-docking preparation
with the original AHC workflow before performing selectivity analysis.

## 10. Distinction between LigPrep and Glide force fields

The SDF property `s_lp_Force_Field = S-OPLS` describes the force field used for
the final LigPrep geometry optimization. It does not state the force field used
inside the subsequent Glide docking calculation.

In this workflow:

```text
LigPrep geometry optimization: S-OPLS (-bff 16)
Glide docking force field:     defined separately in the Glide input
```

Using S-OPLS in LigPrep and OPLS4 in Glide is therefore not contradictory;
they belong to separate computational stages. The Glide settings should be
reported independently in the docking-method section.

## 11. Interpretation for state-selectivity analysis

Following protocol alignment, the following paired comparison is justified:

```text
Original own-state AHC best score
vs
AHC-matched LigPrep opposite-state cross-docking best score
```

For a PPS-origin molecule:

```text
selectivity_margin = PR_grid_score - PPS_grid_score
```

For a PR-origin molecule:

```text
selectivity_margin = PPS_grid_score - PR_grid_score
```

Because lower docking scores are better, a positive `selectivity_margin`
indicates preference for the population's original/own receptor state.

The initial five-molecule PPS-to-PR matched test produced positive margins for
all five molecules (`+6.01` to `+7.94`), providing a preliminary PPS-preference
signal. This small test verifies the workflow but is not sufficient by itself
for a population-level biological conclusion.

## 12. Suggested report-ready Methods text

> Ligands were prepared using the same LigPrep protocol employed during the
> original MolScore AHC runs. Ionization and tautomeric states were generated
> with Epik at a target pH of 7.0 with a tolerance of 1.0 pH unit. Final LigPrep
> geometry optimization used S-OPLS (`-bff 16`), and stereoisomer generation was
> limited using `-s 8`. Prepared variants were docked independently, and the
> lowest (most favourable) Glide docking score across all variants returned for
> each input molecule was retained as the molecule-level score.

Suggested validation sentence:

> Protocol consistency was confirmed both from the installed MolScore command
> construction and empirically from SDF metadata. Thirty sampled AHC output
> structures and the matched cross-docking test outputs consistently reported
> mmshare version 66139, S-OPLS LigPrep optimization, and identical embedded
> Epik settings (`pH 7.0`, `pH tolerance 1.0`).

Suggested selectivity definition:

> State selectivity was quantified as the opposite-state docking score minus
> the own-state docking score. Since more negative Glide scores indicate more
> favourable docking, positive values indicate preferential docking to the
> molecule's original target state.

## 13. Information to preserve for reproducibility

The following should be archived with the final analysis:

- AHC run JSON configuration;
- exact MolScore `_ligand_preparation.py` version used;
- `$SCHRODINGER` installation version/path;
- LigPrep `-long_help` excerpt defining `-bff 16`;
- LigPrep command and preparation mode recorded by the cross-docking runner;
- original and cross-docking manifests;
- all-pose and best-score CSV files;
- receptor grid filenames and Glide input templates;
- `s_lp_Force_Field`, `s_epik_cmdline`, and `i_lp_mmshare_version` properties;
- failure/unavailable counts and variant-to-pose reconciliation results.

## 14. Final conclusion

Three independent evidence sources establish LigPrep consistency:

1. the MolScore source constructs LigPrep with pH 7.0, tolerance 1.0,
   `-bff 16`, `-s 8`, and Epik;
2. the exact server LigPrep help defines `-bff 16` as S-OPLS;
3. both original AHC outputs and updated cross-docking outputs empirically
   report the same mmshare version, S-OPLS force-field metadata, and embedded
   Epik settings.

The updated AHC-matched cross-docking results are therefore suitable for
paired own-state versus opposite-state docking-score analysis, subject to the
usual limitations of docking scores as computational ranking metrics rather
than direct experimental binding free energies.
