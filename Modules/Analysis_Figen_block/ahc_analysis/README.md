# Completed AHC result analysis

Run this module through its single user-facing entry point:

```powershell
python run_ahc_analysis.py --config ahc_analysis.in
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
