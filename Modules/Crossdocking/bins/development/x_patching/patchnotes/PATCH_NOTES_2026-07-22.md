# Cross-docking pipeline patch notes

Date: 2026-07-22

## Objective

Build a reusable pipeline for identifying myosin state-selective molecules by
cross-docking two independently generated AHC populations:

```text
PPS AHC population -> filtering -> docking to PR grid
PR AHC population  -> filtering -> docking to PPS grid
                                -> comparative selectivity analysis
```

The pipeline was designed to preserve the connection to the original AHC run,
avoid unnecessary copying of ligand files, support small terminal tests, and
later scale to both populations.

## Development history

### 1. Existing project and AHC output inspection

- Reviewed the existing project context, development code, and previous AHC
  result-processing package.
- Inspected the PPS test run directory and confirmed that `scores.csv` is the
  selection/filtering source.
- Confirmed the relationship between a selected row and its corresponding AHC
  Glide output:

  ```text
  <run_dir>/<run_name>_GlideDock/<step>/<best_variant>_lib.sdfgz
  ```

- Confirmed that an AHC `*_lib.sdfgz` generally contains a prepared ligand
  pose and can be decompressed directly to SDF.
- Incomplete or unusually small files near the end of the locally downloaded
  run were treated as transfer artifacts rather than a change in directory
  format.

### 2. Filtering decisions

- Kept the existing filtering implementation because its iteration-level
  reporting was useful and easier to read.
- Retained duplicate and validity filtering.
- Treated docking score `0` as docking failure because valid docking scores are
  negative and lower values are better.
- Kept PAINS and BRENK as masks rather than mandatory removal filters.
- Added reporting of the number of PAINS- and BRENK-masked instances.

Filtering remains an upstream stage and produces the selected CSV consumed by
the cross-docking package.

### 3. Glide CLI environment confirmation

- Confirmed the server Schrodinger installation:

  ```text
  SCHRODINGER=/home/apps/schrodinger2024-2
  ```

- Confirmed that both executables are available:

  ```text
  $SCHRODINGER/glide
  $SCHRODINGER/ligprep
  ```

- Decided to reuse the existing Glide CLI workflow rather than introduce a new
  docking implementation.
- Confirmed that a single docking job does not require Dask. Dask remains a
  possible later orchestration layer for larger parallel runs.

### 4. Package reorganization

The Crossdocking package was reorganized into preparation and docking stages:

```text
Crossdocking/
├── single_case_docking.py
├── single_arm_docking.py
├── double_arm_docking.py
├── s1_Prep/
│   ├── fetch_target.py
│   ├── unpack_p201.py
│   └── prepare_ligand.py
├── s2_Dock/
│   ├── make_glide_input.py
│   ├── run_glide.py
│   └── extract_scores.py
└── x_patching/
```

Development and test artifacts were separated from the main reusable modules.

### 5. SDFGZ unpacking module

`s1_Prep/unpack_p201.py` was established as a content-preserving SDFGZ
decompression utility.

It now:

- accepts one `.sdfgz` file;
- writes an `.sdf` without changing its contents;
- records console and file logs;
- refuses to overwrite an existing output;
- removes partial output after decompression failure;
- counts `$$$$` records and rejects an output containing no complete SDF
  record;
- closes and replaces old logging handlers safely when called repeatedly.

### 6. Target manifest instead of ligand copying

`s1_Prep/fetch_target.py` was introduced to connect a filtered CSV to the
corresponding AHC ligand files.

The selected source files are not copied into a temporary directory. Instead,
the script writes a manifest containing paths to the original SDFGZ files.
This avoids redundant disk usage while preserving traceability.

The manifest records:

- source run and molecule ID;
- expected `*_lib.sdfgz` path;
- whether the file exists;
- file size;
- whether it is ready for downstream processing.

Duplicate best-variant IDs are rejected because they would make job identity
ambiguous. `--allow-unavailable` permits manifest inspection or partial runs
when some downloaded source files are missing.

### 7. Standalone single-case docking

`single_case_docking.py` was created for small manual tests outside the batch
pipeline.

It:

- accepts one `.sdf` or `.sdfgz` ligand and one Glide grid ZIP;
- automatically unpacks SDFGZ input;
- does not use Dask;
- does not run LigPrep;
- creates the Glide input;
- runs Glide directly;
- extracts docking scores to CSV;
- rejects output directories containing stale `*_lib.sdfgz` results.

### 8. Single-arm sequential pipeline

The original runner was renamed `single_arm_docking.py` and established as the
terminal runner for one cross-docking direction.

Example arm:

```text
filtered PPS population -> PR receptor grid
```

Each selected instance receives an isolated directory:

```text
jobs/0001_<molecule_id>/
```

The runner performs:

```text
fetch target -> prepare ligand -> create Glide input -> run Glide
             -> extract all poses -> update manifest and combined scores
```

Implemented operational features include:

- sequential execution without Dask;
- `--limit` for small tests;
- `--allow-unavailable` for incomplete inputs;
- `--resume` for interrupted runs;
- `--fail-fast` for stopping on the first job failure;
- continued processing by default when an individual docking job fails;
- per-ligand status and error reporting in `target_manifest.csv`;
- progressive manifest and score updates after each job.

### 9. Double-arm pipeline

`double_arm_docking.py` was added to run the complete bidirectional experiment:

```text
L_PPS_to_R_PR
L_PR_to_R_PPS
```

The original molecule ID is retained rather than renaming the source SDFGZ
file. Collision-safe metadata is added to combined outputs:

```text
molecule_id    = 2_3-2
ligand_state   = PPS
receptor_state = PR
crossdock_arm  = L_PPS_to_R_PR
crossdock_id   = L_PPS_to_R_PR__2_3-2
```

This prevents collisions when both independent AHC populations contain the
same local variant identifier.

Double-arm outputs include:

```text
L_PPS_to_R_PR/
L_PR_to_R_PPS/
double_arm_manifest.csv
double_arm_summary.csv
double_arm_docking_scores.csv
double_arm_best_docking_scores.csv
```

### 10. Optional LigPrep route

Both arm runners now accept:

```text
--ligand-prep-mode off
--ligand-prep-mode on
```

`off` is the default and preserves the original cross-docking route:

```text
existing AHC *_lib.sdfgz -> unpack -> opposite-state docking
```

`on` starts again from the filtered CSV SMILES:

```text
filtered SMILES -> per-instance SMI -> LigPrep variants
                -> dock every prepared variant
```

`--smiles-col` can be used when the filtered CSV does not name its SMILES
column `smiles`.

For LigPrep mode, one prepared SDF may contain multiple protonation,
tautomeric, or stereochemical variants. Glide docks all prepared records.
Every returned pose is retained, and the most negative docking score across
all variants and poses is selected as the representative score for that
filtered instance.

Outputs are separated into:

```text
combined_docking_scores.csv  # all variants and returned poses
best_docking_scores.csv      # one lowest-scoring row per selected instance
```

Each job also contains `docking_scores.csv` and `best_docking_score.csv`.
The manifest records the preparation mode, SMILES, prepared file, variant
count, LigPrep status, and best docking score.

### 11. Testing and server observations

- Python syntax and CLI import checks passed in the local DL environment.
- Manifest-only smoke tests confirmed unavailable-target handling.
- Mock docking tests confirmed both LigPrep modes and best-score selection.
- A mock instance with scores `-7.0`, `-9.0`, and `-8.0` correctly selected
  `-9.0` as its representative result.
- Double-arm combination tests confirmed arm labels and collision-safe IDs.
- The first server LigPrep test exposed Schrodinger's requirement that the
  relative output belongs to the launch working directory. The LigPrep launch
  was adjusted to use a relative output filename with the job's LigPrep output
  directory as its working directory.
- After the launch-path adjustment, LigPrep and Glide proceeded successfully.
- Initial live output confirmed that multiple LigPrep variants were retained
  and docked. For example, one selected instance produced four prepared
  variants and four returned docking poses.
- Schrodinger reported that the configured shared scratch directory was not
  writable and automatically fell back to:

  ```text
  /home/andy/.schrodinger/tmp
  ```

  This is non-fatal, but available disk space should be checked before a large
  production run.

## Current status

The functional pipeline is now available for:

1. standalone single-ligand testing;
2. one-direction population cross-docking;
3. complete PPS-to-PR and PR-to-PPS cross-docking;
4. optional reuse of AHC prepared poses or fresh LigPrep from SMILES;
5. preservation of all pose-level results and generation of one best score per
   selected instance.

The current server test should be allowed to finish before final validation.
Final validation should compare:

- terminal completion counts;
- `target_manifest.csv` statuses;
- row counts in `combined_docking_scores.csv`;
- one row per completed instance in `best_docking_scores.csv`;
- LigPrep variant counts versus returned docking pose counts;
- available space under `/home/andy/.schrodinger/tmp`.

The next development stage is comparative state-selectivity analysis using the
two arm-level best-score tables.



# working 

Attempted to change the log output by
` run_env = os.environ.copy() ` \
` run_env["PWD"] = str(output_sdf.parent) `

Attempted to match the lig prep conditions to the AHC run variables.
