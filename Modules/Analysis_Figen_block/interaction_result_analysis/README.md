# Interaction result analysis

This module only post-processes existing interaction result files.

```powershell
python run_interaction_result_analysis.py --config interaction_result_analysis.in
```

It does **not** run ProLIF, PoseCheck, receptor preparation, or pose extraction.
Those calculations belong to `Modules/Crossdocking/CompletedAnalysis/interaction`.

## Blocks

- `residue_hotspot`: residue, interaction-type, and dominant-type prevalence.
- `full_last`: full-population versus last-window hotspot tables.
- `top_candidates`: own/opposite prevalence for selected candidates.
- `own_opposite`: matched own-state and opposite-state hotspot tables.

At least one block must be ON. Each enabled block validates only its required
inputs. Named counts use comma-separated `NAME=positive_integer` values, for
example `PPS=200,PR=200`.

Paths are resolved relative to the `.in` file. Outputs remain below
`work_results/analysis/<feature>/<run_id>/`; `resume=ON` versions related output
tables together instead of overwriting them.

The retired adapter that launched interaction calculations is preserved only at
`bins/legacy_interaction_adapter/` for historical comparison and is not exposed
by the canonical workflow registry.
