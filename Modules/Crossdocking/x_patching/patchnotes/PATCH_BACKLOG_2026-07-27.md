# Cross-docking patch backlog

Date: 2026-07-27

## Context

This backlog was recorded after reviewing the production cross-docking run
configured with the PPS and PR `scores.csv` files as both:

1. raw inputs for filtering and selected cross-docking; and
2. raw inputs for population-overlap analysis.

The running production job is not to be modified. These items are follow-up
changes for a later patch and validation cycle.

## Patch items

### 1. Remove redundant raw score-path configuration

Current configuration requires the same raw files to be repeated:

```text
pps_input=.../PPSrun/scores.csv
pr_input=.../PRrun/scores.csv
pps_scores=.../PPSrun/scores.csv
pr_scores=.../PRrun/scores.csv
```

`pps_input` and `pr_input` are used for population filtering, while
`pps_scores` and `pr_scores` are passed to raw population-overlap analysis.
In the current full workflow these are normally the same source files.

Planned behavior:

- use `pps_input` and `pr_input` as the default population-overlap inputs;
- retain optional `pps_scores` and `pr_scores` only if an explicit override is
  useful for backward compatibility;
- if override keys remain supported, validate them as a complete pair;
- update configuration examples and documentation;
- record the resolved overlap input paths in `pipeline_run_summary.json`.

### 2. Reconcile the documented `limit` option

The README lists `limit` as a full-run configuration option, but
`crossdocking_runner.py` does not parse or forward it. The top-level runner
always passes `limit=None` to the double-arm pipeline.

Planned behavior:

- treat `testmode_top_n` as the supported full-workflow selection override;
- remove the misleading top-level `limit` entry from the README, unless a
  distinct post-filter docking limit is intentionally required;
- keep the lower-level single-arm `--limit` option documented only in its
  appropriate CLI scope.

### 3. Add regression coverage for configuration resolution

Add tests covering:

- omitted `pps_scores` and `pr_scores` defaulting to `pps_input` and
  `pr_input`;
- explicit paired overlap-path overrides;
- rejection of a one-sided override;
- blank `testmode_top_n` selecting percentage mode;
- documentation/config examples matching the actual top-level parser.

## Validation required before deployment

- run the filtering stage against small PPS and PR fixtures;
- confirm overlap analysis reads the resolved raw inputs rather than filtered
  CSV files;
- run a small `testmode_top_n` double-arm test;
- compare output schemas and overlap counts with the pre-patch behavior;
- do not replace or resume the currently running production output directory
  during validation.
