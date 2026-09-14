# Completed AHC analysis

Run the completed PR/PPS analysis through one entry point:

```powershell
python run_analysis.py --config analysis.in
```

Only `run_analysis.py` is user-facing. Files under `blocks/` are internal
orchestration units built on the validated analysis and figure workflows.

## Required configuration

| Key | Meaning |
|---|---|
| `pr_input` | Completed PR score CSV |
| `pps_input` | Completed PPS score CSV |
| `results_root` | Must be a directory named `work_results` |
| `run_id` | Shared identifier used by all enabled blocks |

## Blocks and order

| Key | Default | Work performed |
|---|---:|---|
| `physchem` | `on` | Three-property PR/PPS distributions, scatter plots, and trajectories |
| `last_window` | `on` | Combined last-N PR/PPS physicochemical comparison |
| `scaffold` | `off` | Consolidated scaffold and fingerprint-diversity tables; membership is prepared automatically |
| `chemical_space` | `off` | Joint UMAP/t-SNE calculation followed by matching figures |

Enabled blocks run in dependency-safe order: physicochemical, last-window,
chemical-space, then scaffold. The runner stops at the first
error and does not silently skip a requested block.

## Options and override rules

- `reference` defaults to the bundled reference table. Physicochemical figures
  use OM by default.
- `last_n` defaults to `50` and affects last-window and scaffold analyses.
- `resume=off` rejects an existing run directory. `resume=on` explicitly adds
  collision-safe suffixes; it never silently overwrites figures.
- `usecache=on` (default) reuses the newest chemical-space cache only when the
  PR, PPS, and reference SHA-256 values and all calculation settings match.
  `usecache=off` always recalculates the requested embeddings.
- `membership` is not a user input. With `scaffold=on`, the runner automatically
  obtains it from the joint UMAP/fingerprint calculation.
- If both `scaffold` and `chemical_space` are on, the chemical-space result is
  reused and `umap` plus `fp` are automatically included if omitted.
- If only `scaffold` is on, the runner performs the minimum UMAP/fingerprint
  calculation internally without generating chemical-space figures.
- `methods` and `features` otherwise affect only `chemical_space=on`.
- `random_state` affects chemical-space and defaults to `42`.
- The completion summary records the decision as `chemical_space_cache_hit`.
- Physicochemical comparison uses PR blue `#2878B5`, PPS red `#D9534F`, removes
  score-zero placeholder rows, uses OM as the default reference, shares the PR/
  PPS molecular-weight axis, and uses integer rotatable-bond ticks from 0–17.

All outputs are written below the configured `work_results` directory.
