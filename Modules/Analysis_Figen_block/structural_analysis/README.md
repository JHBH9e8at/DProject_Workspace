# Structural analysis

Run the enabled structural blocks through one entry point:

```bash
python Modules/Analysis_Figen_block/structural_analysis/run_structural_analysis.py \
  --config Modules/Analysis_Figen_block/structural_analysis/structural_analysis.in
```

`pocket_volume` and `redocking` are independently selectable with ON/OFF, but
at least one must be enabled. Implementation files under `blocks/` are internal.

## Configuration reference

Shared options:

| Key | Required | Default | Meaning |
|---|---:|---:|---|
| `results_root` | Yes | — | Output root; its directory name must be `work_results`. |
| `run_id` | Yes | — | Identifier appended to this module's analysis and figure directories. |
| `resume` | No | `OFF` | `OFF` rejects an existing run; `ON` creates collision-safe output names. |
| `pocket_volume` | No | `OFF` | Enable ligand-local Fpocket alpha-sphere volume estimation. |
| `redocking` | No | `OFF` | Enable PR/PPS docking-score and RMSD validation. |

Pocket-volume options (required only when `pocket_volume=ON`):

| Key | Required | Default | Meaning |
|---|---:|---:|---|
| `fpocket_dir` | Yes | — | One Fpocket `<protein>_out` directory containing pocket files. |
| `targets` | Yes | — | Comma-separated Fpocket pocket IDs, for example `2,38,51`. |
| `ligand` | Yes | — | Reference-ligand PDB used to select nearby alpha spheres. |
| `ligand_resname` | No | blank | Restrict ligand atoms to this PDB residue name; blank uses all parsed ligand atoms. |
| `distance_threshold` | No | `3.0` | Maximum centre-to-heavy-atom distance in Å for selecting a sphere. |
| `all_spheres` | No | `OFF` | `ON` integrates all spheres in each target pocket; `OFF` keeps only ligand-local spheres. |
| `iterations` | No | `1000000` | Number of Monte Carlo sample points. More points reduce sampling noise but increase runtime. |
| `seed` | No | `42` | Random seed making the Monte Carlo estimate reproducible. |
| `chunk_size` | No | `100000` | Points evaluated per memory batch; it changes memory/runtime behaviour, not the intended estimator. |
| `radius_offset` | No | `-1.6` | Value in Å added to every raw Fpocket alpha-sphere radius before integration. |
| `output_prefix` | No | `ligand_local_pocket_volume` | Base name for the generated CSV and JSON tables. |

Redocking options (required only when `redocking=ON`):

| Key | Required | Default | Meaning |
|---|---:|---:|---|
| `pps_docking_csv`, `pr_docking_csv` | Yes | — | PPS/PR redocking-score tables. |
| `pps_rmsd_csv`, `pr_rmsd_csv` | Yes | — | Matching PPS/PR RMSD tables in the same pose order. |
| `score_ylim` | No | `-12,-1` | Lower and upper y-axis bounds for docking-score plots. |
| `rmsd_ylim` | No | `0,2.5` | Lower and upper y-axis bounds for RMSD plots. |

Relative input paths are resolved from `structural_analysis.in`. Unknown keys,
invalid ON/OFF values, non-positive `iterations`/`chunk_size`, and reversed plot
bounds are rejected.

## Pocket volume

`blocks/pocket_volume.py` contains pure parsing, sphere selection, and Monte Carlo
volume calculations. It does not create output files.

For alpha-sphere centre `c_i`, ligand heavy-atom coordinate `l_j`, and distance
threshold `d`, a sphere is selected when:

```text
min[j] ||c_i − l_j|| ≤ d
```

The current default is `d = 3.0 Å`. This is a centre-to-atom criterion. The
effective radius used for integration is:

```text
r_eff = r_raw + radius_offset
```

The current `radius_offset=-1.6 Å` reproduces the relevant Fpocket volume
convention. The union volume is estimated by:

```text
V_union(est.) = V_box × N_inside / N
SE(V_est.)    = V_box × sqrt[p(1 − p) / N]
```

A point inside multiple spheres is counted once. The combined result is one
union over all selected pockets and is not generally the arithmetic sum of the
individual pocket volumes.

`blocks/pocket_volume_runner.py` is the internal output boundary. It writes to:

```text
work_results/analysis/pocket_volume/<run_id>/
|-- inputs/
|-- intermediate/
|-- tables/
|   |-- <output_prefix>.csv
|   `-- <output_prefix>.json
|-- logs/
`-- run_manifest.json
```

Example from the `work_space` root:

```text
python -m Modules.Analysis_Figen_block.structural_analysis.pocket_volume_runner \
  --fpocket-dir <fpocket-result-directory> \
  --target 2,38,51 \
  --ligand <reference-ligand.pdb> \
  --target-label PPS \
  --output-prefix PPS_reference_3A
```

If `--run-id` is omitted, a timestamped run ID is generated. Existing runs
are rejected unless `--resume` is explicit. A resumed run allocates matching
CSV/JSON suffixes rather than overwriting existing scientific results.

The configured `targets=2,38,51` are specific to the current
`PPS_Protein_out` Fpocket result. Pocket IDs must be rechecked after rerunning
Fpocket or changing the receptor structure.

## Redocking

When `redocking=ON`, all four PR/PPS docking-score and RMSD CSV paths are
required. Tables are written below `work_results/analysis/redocking/`; PNG/SVG
pairs are written below `work_results/Figs/redocking_validation/`. `score_ylim`
and `rmsd_ylim` each accept two comma-separated bounds.

Docking-score and RMSD rows are paired by position, so the input files must use
the same pose ordering. A successful redocking pose is defined here as
`RMSD < 2.0 Å`; the reported fraction is the number passing this threshold
divided by the total pose count.

The docking CSV must contain a numeric column named `docking score`; the RMSD
CSV must contain a numeric column named `RMS`. For each state, the two files
must contain the same nonzero number of rows, and the PR and PPS pose counts
must also match. The success fraction is:

```text
f_success = N(RMSD < 2.0 Å) / N(all paired poses)
```

Relative paths are resolved from `structural_analysis.in`. `resume=OFF` rejects
an existing run; `resume=ON` versions related tables and figures together.
