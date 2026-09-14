# Completed cross-docking analysis

Run completed score and optional interaction analysis with:

```bash
python -m Modules.Crossdocking.CompletedAnalysis.crossdocking_analysis \
  --config Modules/Crossdocking/crossdocking.in
```

Run this command from the repository root. It consumes an existing completed
cross-docking calculation; it does not launch Glide.

Enabled operations are controlled by `crossdocking.in`:

- `score_analysis`: validate and pair own/opposite-state scores, classify
  selectivity, and analyse raw-population overlap;
- `interaction_analysis`: run PoseCheck and/or ProLIF on cross-docked poses;
- `reference_interaction_analysis`: run ProLIF for the PPS and PR reference
  complexes and create `reference_interactions.csv`.

```text
selectivity_margin = opposite_state_score − own_state_score
```

Because more-negative Glide scores are more favourable, a positive margin
supports preference for the original AHC state.

`multiprocessing` affects only PoseCheck and ProLIF. Glide cross-docking remains
sequential by design because concurrent external Glide jobs may exhaust licence
slots or destabilize the execution host.

With `posecheck_chunk_size=auto`, each receptor-state population uses:

```text
ceil(valid_pose_count / posecheck_workers)
```

The minimum chunk size is one. A positive integer can be supplied instead.

Principal outputs are written below:

```text
work_results/analysis/
├── crossdocking_scores/<run_id>/
├── crossdocking_interactions/<run_id>/
├── reference_interactions/<run_id>/reference_interactions.csv
└── completed_crossdocking/<run_id>/crossdocking_analysis.json
```

The current production defaults use 30 PoseCheck workers, automatic chunk
sizing, and 30 ProLIF workers. Reduce worker counts if the host reports memory
pressure or `BrokenProcessPool`.

## Configuration reference

This runner reuses the completed-analysis section of `crossdocking.in`.

| Key | Required | Default | Meaning |
|---|---:|---:|---|
| `results_root` | Yes | — | Root `work_results` directory for analysis outputs. |
| `run_id` | Yes | — | Identifier shared across score, interaction, reference, and completion outputs. |
| `output_dir` | Normally | — | Completed calculation directory. When `calculation_manifest` is omitted, the manifest is read from this directory. |
| `calculation_manifest` | Conditional | derived | Explicit calculation-manifest path; required if it cannot be derived from `output_dir`. |
| `score_analysis` | No | `ON` | Pair own/opposite scores and calculate selectivity. |
| `interaction_analysis` | No | `OFF` | Run cross-docked-pose analysis; requires both receptor PDBs. |
| `reference_interaction_analysis` | No | `OFF` | Run reference ProLIF; requires both receptors and both reference ligands. |
| `pps_receptor`, `pr_receptor` | Conditional | — | Receptor PDBs for PoseCheck/ProLIF. |
| `pps_reference_ligand`, `pr_reference_ligand` | Conditional | — | Reference poses for reference interaction analysis. |
| `interaction_mode` | No | `both` | `posecheck`, `prolif`, or `both`. |
| `include_secondary_interactions` | No | `OFF` | Include supported secondary interaction categories. |
| `multiprocessing` | No | `OFF` | Parallelise completed interaction analysis only. |
| `posecheck_workers`, `prolif_workers` | No | `30` | Positive worker-process counts. |
| `posecheck_chunk_size` | No | `auto` | `auto` or a positive number of poses per task. |
| `testmode_top_n` | No | blank | Optional positive cap for a small validation run. |
| `selectivity_threshold` | No | `2.0` | Minimum positive score difference for the configured selective class. |
| `strong_score_threshold` | No | `-8.0` | Score boundary used to annotate strong-scoring candidates. |
| `resume` | No | `OFF` | Permit eligible completed-output reuse/versioning. |

Relative paths resolve from `crossdocking.in`. The calculation manifest is the
provenance boundary; do not combine it with loose CSV files from another run.

## Selectivity classification

Let `S_own` be the AHC score against the generation state and `S_opp` the best
cross-docking score against the opposite receptor:

```text
margin = S_opp - S_own
```

For selectivity threshold `delta` and strong-score threshold `S_strong`, the
classes are assigned as follows:

```text
margin >= delta and S_own <= S_strong  -> selective_own_state
margin >= delta                         -> own_state_preferred_weak
margin <= -delta and S_opp <= S_strong -> reverse_selective
margin <= -delta                        -> reverse_preferred_weak
both scores <= S_strong                 -> nonselective_strong
otherwise                               -> inconclusive
```

Ranking first follows this class priority, then decreasing margin, then the
more-negative own-state score. This ordering is a prioritisation rule, not a
statistical significance test.
