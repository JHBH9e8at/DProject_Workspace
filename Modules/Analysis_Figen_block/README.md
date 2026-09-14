# Analysis and figure block

Run all selected modules through one entry point:

```powershell
python run_analysis.py --config analysis.in
```

The top-level `.in` enables modules and points to each module's own `.in` file.
The fixed order is AHC analysis, structural analysis, interaction-result
analysis, and standalone figures. A failure stops the sequence immediately.

`run_id` and `results_root` in the top-level file override the corresponding
values in enabled child configs. Other child options are never overridden.
This makes the same child config usable both independently and as part of a
combined run.

Use `dry_run=ON` to validate the top-level configuration and print the resolved
execution plan without starting a calculation. The provided template defaults
to dry-run for safety; change it to OFF for a real run.

Subdirectories expose only their runner, `.in`, README, and internal `blocks/`.
`common/` is shared infrastructure and is not directly executable. Historical
or duplicated code is retained below `bins/` and is not part of the canonical
workflow.
