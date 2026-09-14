# Completed AHC analysis

Run the completed PR/PPS analysis through one entry point:

```bash
python Modules/ACHresultspackage/CompletedAnalysis/run_analysis.py \
  --config Modules/ACHresultspackage/CompletedAnalysis/analysis.in
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
- `usecache=on` (default) automatically reuses an existing chemical-space
  result only when it was generated from the same PR, PPS, and reference file
  contents with the same calculation settings. File identity is checked with
  SHA-256 checksums; the user does not enter these values manually.
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

Fingerprint diversity is defined as:

```text
D_FP = 1 − mean pairwise Tanimoto similarity
```

Scaffold summaries use Bemis-Murcko scaffolds and report unique and singleton
fractions, top-scaffold shares, Shannon entropy, and normalized entropy. UMAP
fits the generated PR/PPS population and transforms reference compounds into
that space. t-SNE has no transform operation and therefore uses a joint fit.

## Complete configuration reference

| Key | Required | Default | Meaning |
|---|---:|---:|---|
| `pr_input`, `pps_input` | Yes | — | Completed PR and PPS AHC score CSVs. |
| `results_root` | Yes | — | Output root named `work_results`. |
| `run_id` | Yes | — | Shared output identifier. |
| `physchem` | No | `ON` | Physicochemical distributions, score scatter, and trajectories. |
| `last_window` | No | `ON` | Comparison over the final `last_n` iterations. |
| `scaffold` | No | `OFF` | Bemis-Murcko scaffold and fingerprint-diversity analysis. |
| `chemical_space` | No | `OFF` | Requested UMAP/t-SNE calculations and figures. |
| `reference` | No | bundled CSV | Reference-compound descriptor table. |
| `last_n` | No | `50` | Positive final-window size. |
| `resume` | No | `OFF` | Existing-run handling. |
| `usecache` | No | `ON` | Reuse a cache after automatic provenance validation. |
| `methods` | No | `umap,tsne` | Comma-separated subset of `umap` and `tsne`. |
| `features` | No | `fp,all_desc,select_desc` | Fingerprint, all-descriptor, and/or selected-descriptor representations. |
| `random_state` | No | `42` | Seed for reproducible stochastic embeddings. |

## Mathematical definitions

For binary fingerprints A and B, Tanimoto similarity is:

```text
T(A,B) = |A intersection B| / |A union B|
       = c / (a + b - c)
```

Fingerprint diversity over all unordered molecule pairs is:

```text
D_FP = 1 - [sum over i<j of T(F_i,F_j)] / [n(n-1)/2]
```

For scaffold proportions `p_i`, Shannon entropy and normalized entropy are:

```text
H      = -sum_i p_i ln(p_i)
H_norm = H / ln(K)
```

Here `K` is the number of observed scaffold classes. Higher values indicate a
more even scaffold distribution, not necessarily biological diversity.

## Expected input tables

Both AHC CSVs must contain `smiles`, `step`, and the state-specific docking
score column (`PR_r_i_docking_score` or `PPS_r_i_docking_score`). Population
cleaning additionally uses `valid` and `unique`. Physicochemical blocks require
their `desc_*` columns; the standard production plots use `desc_MolWt`,
`desc_NumRotatableBonds`, and `desc_CLogP`. Chemical-space descriptor modes
require the descriptor set selected by the implementation.

The reference CSV must contain `name` and `smiles` plus the descriptor columns
needed by enabled comparisons. Score zero is treated as an unavailable docking
result and excluded from score-based numerical summaries.
