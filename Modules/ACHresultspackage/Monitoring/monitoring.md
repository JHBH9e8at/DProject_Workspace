# AHC in-run monitoring

Use one Python entry point for a normal monitoring run:

```powershell
python run_monitor.py --config t_mp.in
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
