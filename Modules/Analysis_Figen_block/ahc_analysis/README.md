# Completed AHC result analysis

Run this module through its single user-facing entry point:

```bash
python Modules/Analysis_Figen_block/ahc_analysis/run_ahc_analysis.py \
  --config Modules/Analysis_Figen_block/ahc_analysis/ahc_analysis.in
```

Implementation files under `blocks/` are internal. The same analysis remains
available through `ACHresultspackage/CompletedAnalysis/run_analysis.py` for the
package-level workflow; both entry points call the same implementation.

## Required configuration

| Key | Meaning |
|---|---|
| `pr_input` | Completed PR score CSV |
| `pps_input` | Completed PPS score CSV |
| `results_root` | Directory whose name is exactly `work_results` |
| `run_id` | Shared output identifier |

## Blocks

The fixed execution order is physicochemical, last-window, chemical-space, and
scaffold. `physchem` and `last_window` default to ON; the two heavier blocks
default to OFF.

- `physchem`: PR/PPS distributions, docking scatter, and trajectory figures.
- `last_window`: combined last-N PR/PPS physicochemical comparison.
- `chemical_space`: joint UMAP/t-SNE calculations and figures.
- `scaffold`: scaffold and fingerprint-diversity tables.

## Important option interactions

- `membership` is never supplied by the user.
- `scaffold=ON` automatically prepares or reuses the required UMAP/fingerprint
  membership.
- `chemical_space=ON` and `scaffold=ON` share one calculation. UMAP and
  fingerprint features are automatically added if scaffold needs them.
- `usecache=ON` reuses only a cache whose input hashes and calculation settings
  match. `usecache=OFF` forces recalculation.
- `resume=OFF` rejects an existing run. `resume=ON` uses collision-safe output
  names and does not overwrite figures.
- Config-relative paths are resolved from `ahc_analysis.in`, not the shell's
  current directory.

Figures preserve the agreed PR blue `#2878B5`, PPS red `#D9534F`, OM reference,
shared molecular-weight axis, and integer rotatable-bond ticks from 0 to 17.

## Calculation notes

- Molecular identity is based on canonical RDKit SMILES.
- For duplicate canonical molecules, the most negative own-state Glide score is
  retained.
- Zero docking scores are treated as unavailable/failed scores.
- Descriptor features are standardized using the generated population before
  dimensionality reduction.
- Fingerprint diversity is `1 − mean pairwise Tanimoto similarity`.
- Scaffold analysis uses Bemis-Murcko scaffolds and reports frequency and
  entropy-based diversity metrics.

Outputs are grouped under `work_results/analysis/` and `work_results/Figs/` by
feature and run ID. The completion record is written below
`analysis/completed_ahc/<run_id>/`.

## Configuration reference

| Key | Required | Default | Meaning |
|---|---:|---:|---|
| `pr_input`, `pps_input` | Yes | — | Completed PR/PPS score CSV files. |
| `results_root` | Yes | — | Output root named `work_results`. |
| `run_id` | Yes | — | Identifier for this completed-analysis run. |
| `physchem` | No | `ON` | Physicochemical, docking-scatter, and trajectory results. |
| `last_window` | No | `ON` | Compare the final `last_n` PR/PPS iterations. |
| `chemical_space` | No | `OFF` | Requested UMAP/t-SNE representations. |
| `scaffold` | No | `OFF` | Bemis-Murcko and fingerprint-diversity tables. |
| `reference` | No | bundled CSV | Reference descriptor table. |
| `last_n` | No | `50` | Final-window size. |
| `resume` | No | `OFF` | Existing-run handling. |
| `usecache` | No | `ON` | Reuse a provenance-compatible chemical-space cache. |
| `methods` | No | `umap,tsne` | Requested embedding algorithms. |
| `features` | No | `fp,all_desc,select_desc` | Requested molecular representations. |
| `random_state` | No | `42` | Stochastic embedding seed. |

For the full Tanimoto, fingerprint-diversity, and scaffold-entropy definitions,
see `Modules/ACHresultspackage/CompletedAnalysis/analysis.md`; both entry points
call the same implementation.

## Expected input tables

The PR/PPS files must contain `smiles`, `step`, `valid`, `unique`, their
state-specific `*_r_i_docking_score`, and the descriptor columns required by
enabled blocks. The reference table requires `name`, `smiles`, and matching
descriptor fields. A docking score of zero is handled as a failed/unavailable
score rather than a physical zero-valued Glide result.
