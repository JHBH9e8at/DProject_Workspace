# Interaction result analysis

This module only post-processes existing interaction result files.

```bash
python Modules/Analysis_Figen_block/interaction_result_analysis/run_interaction_result_analysis.py \
  --config Modules/Analysis_Figen_block/interaction_result_analysis/interaction_result_analysis.in
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

`residue_hotspot` deduplicates repeated ProLIF occurrences at the molecule and
residue level before counting. Its principal calculation is:

```text
residue prevalence
    = unique molecules interacting with the residue
      / analysed molecules in the population
```

Type-specific prevalence additionally groups by ProLIF interaction type. The
dominant type is the type with the largest prevalence for a residue and state.
The configured denominators must match the populations actually analysed; they
must be updated when cross-docking selection or pose availability changes.

Paths are resolved relative to the `.in` file. Outputs remain below
`work_results/analysis/<feature>/<run_id>/`; `resume=ON` versions related output
tables together instead of overwriting them.

Historical implementations under `bins/` are not exposed by the canonical
workflow.

## Configuration reference

Shared options:

| Key | Required | Meaning |
|---|---:|---|
| `results_root` | Yes | Output root whose directory name is `work_results`. |
| `run_id` | Yes | Output identifier; the top-level runner may override it. |
| `resume` | No | `OFF` rejects existing outputs; `ON` allocates versioned names. |
| `residue_hotspot`, `full_last`, `top_candidates`, `own_opposite` | No | Independent ON/OFF block switches; at least one must be ON. |

Inputs are conditionally required as follows:

| Enabled block | Required keys | What the inputs represent |
|---|---|---|
| `residue_hotspot` | `interactions_csv`, `denominators` | Long-form ProLIF occurrences and the analysed molecule count for each population. |
| `full_last` | `interactions_csv`, `metadata_csv`, `total_poses`, `last_start` | Full interaction table, pose ordering/identity, total population sizes, and the first one-based position included in each last window. |
| `top_candidates` | `summary_csv`, `transitions_csv` | Candidate summary and state-transition tables from upstream interaction analysis. |
| `own_opposite` | `own_metadata`, `own_interactions`, `opposite_metadata`, `opposite_interactions` | Matched pose metadata and interaction occurrences for own-state and opposite-state docking. |

Named integer mappings use `NAME=value` pairs separated by commas:

```ini
denominators=PPS=205,PR=197
total_poses=PPS=205,PR=197
last_start=PPS=156,PR=148
```

Counts must be positive integers. Blank paths are allowed only when the block
that requires them is OFF.

## Calculation details

For state `s` and residue `r`, repeated occurrences are first collapsed to one
molecule–residue observation. The reported prevalence is:

```text
P(r | s) = N_unique_molecules_interacting(r, s) / N_analysed_molecules(s)
```

For an interaction type `t`, the same calculation is performed after grouping
by `(state, residue, type)`:

```text
P(r, t | s) = N_unique_molecules_interacting(r, t, s) / N_analysed_molecules(s)
```

Therefore the denominator is a molecule count, not an interaction-occurrence
count. A molecule producing several hydrogen-bond rows for the same residue is
counted once for residue prevalence. Supplying a denominator larger or smaller
than the actual analysed population directly biases every reported percentage.
