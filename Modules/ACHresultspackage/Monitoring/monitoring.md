# AHC in-run monitoring

Use one Python entry point for a normal monitoring run:

```bash
python Modules/ACHresultspackage/Monitoring/run_monitor.py \
  --config Modules/ACHresultspackage/Monitoring/t_mp_pr.in
python Modules/ACHresultspackage/Monitoring/run_monitor.py \
  --config Modules/ACHresultspackage/Monitoring/t_mp_pps.in
```

Scripts under `blocks/` are internal implementation units and normally should
not be run individually.

## Execution order

1. Merge completed `*_scores.csv` iteration files, excluding the newest file.
2. Clean and canonicalise the merged population.
3. Generate physicochemical distribution and docking-score scatter figures.
4. Generate the docking-score trajectory figure.
5. Generate UMAP outputs when enabled.
6. Generate t-SNE outputs when enabled.

The newest iteration file is intentionally excluded because the active AHC
process may still be writing it.

At least two `*_scores.csv` files are required. With files `000001` through
`000354`, the runner merges `000001` through `000353` and excludes `000354`.

## Configuration

The `.in` format is one `key=value` per line. Blank lines and `#` comments are
ignored. Relative paths are resolved from the `.in` file directory. Duplicate
keys are rejected.

### Required keys

| Key | Meaning |
|---|---|
| `analysismode` | `mp` for active iterations; `fa` for finished-CSV compatibility |
| `basedir` | Required for `mp`; directory containing `*_scores.csv` files |
| `indir` | Required for `fa`; completed score CSV |
| `odir` | Output directory |
| `jn` | `PR` or `PPS` |
| `threshold` | Docking threshold, normally PR `-6.0` or PPS `-9.0` |

### Optional keys

| Key | Default | Values / meaning |
|---|---:|---|
| `ref` | `OM` | Comma-separated `OM`, `MAV`, `AFI`; empty disables references |
| `ref_csv` | bundled table | Override reference-descriptor CSV |
| `umap` | `off` | `on` or `off` |
| `tsne` | `off` | `on` or `off` |
| `mode` | `both` | `all`, `select`, or `both` |
| `feat` | `both` | `fp`, `desc`, or `both` |

## Override and dependency rules

- `analysismode=mp` uses `basedir`; `indir` is ignored.
- `analysismode=fa` uses `indir`; `basedir` is ignored.
- `ref_csv` overrides the bundled reference table.
- `mode` and `feat` affect only UMAP/t-SNE.
- When `umap=off`, UMAP-specific settings are ignored.
- When `tsne=off`, t-SNE-specific settings are ignored.
- Disabling UMAP and t-SNE does not disable merge, cleaning, physicochemical
  figures, or the docking trajectory.

The `fa` route remains temporarily for compatibility. Finished-run analysis
will be exposed through its own top-level runner.

## Outputs

```text
<odir>/
├── <job>_<merged_iteration_count>_<timestamp>.csv
├── <merged_name><job>_cleaned.csv
└── figures/
```

The cleaning stage retains valid and unique rows, canonicalises SMILES, and
keeps the best docking-score instance of each canonical molecule. A score of
zero represents a failed or unavailable docking result rather than a physical
score. UMAP and t-SNE are optional; disabling them does not disable the core
physicochemical and trajectory figures.

## Expected score input

In `mp` mode, `basedir` should contain iteration files named like:

```text
000001_scores.csv
000002_scores.csv
000003_scores.csv
...
```

The files must retain the AHC score-table fields used by the original analysis:
SMILES, validity/uniqueness indicators, iteration or step identity, and receptor
docking score. `jn` selects the PR or PPS score convention.

## Calculation interpretation

Canonical SMILES define molecule identity. When retained rows map to the same
canonical molecule, the most favourable (most negative) score is kept:

```text
representative(molecule) = arg min[row] docking_score(row)
```

The threshold defines the configured favourable-score region. UMAP and t-SNE
coordinates are similarity embeddings, not physical coordinates or binding
energies. When descriptor scaling is used, each feature is standardized as:

```text
z = (x - mean(x)) / standard_deviation(x)
```
