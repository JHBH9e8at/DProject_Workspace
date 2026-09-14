# Analysis and figure workflows

## Redocking validation

The redocking workflow consumes four explicit CSV files:

- PPS docking scores;
- PPS RMSD values;
- PR docking scores;
- PR RMSD values.

It creates two linked runs with the same run ID:

```text
work_results/analysis/redocking/<run_id>/
`-- tables/
    |-- redocking_poses.csv
    `-- redocking_summary.csv

work_results/Figs/redocking_validation/<run_id>/
|-- redocking_docking_score_comparison.png/.svg
|-- redocking_rmsd_comparison.png/.svg
|-- redocking_PR_mirrored.png/.svg
`-- redocking_PPS_mirrored.png/.svg
```

Example:

```text
python -m Modules.Analysis_Figen_block.workflows.redocking_validation \
  --pps-docking <PPS_docking.csv> \
  --pps-rmsd <PPS_RMSD.csv> \
  --pr-docking <PR_docking.csv> \
  --pr-rmsd <PR_RMSD.csv>
```

The figure manifest records the pose table as its direct input and the shared
run ID as its upstream analysis reference. Existing runs require explicit
`--resume`; resumed tables and all eight figure outputs receive matching
numeric suffixes without overwriting prior outputs.

