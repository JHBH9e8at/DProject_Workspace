# Common workflow infrastructure

This package is the shared runtime foundation used by Analysis_Figen_block,
ACHresultspackage, and Crossdocking. It is not a user-facing runnable module.
It centralizes strict configuration parsing, result paths, isolated run
directories, collision-safe figure names, and provenance manifests.

## Components

- `config.py`: strict UTF-8 `key=value` parsing, required-key validation, and
  typed ON/OFF, integer, and numeric conversion.
- `paths.py`: project discovery, config-relative path resolution, and the
  `work_results` safety boundary.
- `output_naming.py`: non-overwriting filename allocation.
- `run_context.py`: standard analysis run directories.
- `provenance.py`: input/output hashes and run manifests.

## Safety rules

- Production outputs belong below the project-level `work_results` directory.
- Callers must use explicit feature components and a validated run ID.
- Existing run directories fail by default; resuming requires `resume=True`.
- Figure/data output groups receive one shared numeric suffix when any requested
  filename already exists.
- Manifest outputs must resolve inside their associated run directory.
- Configuration paths resolve relative to the `.in` file, not the process
  working directory.
- Inline comments are not supported; comments must occupy their own line.

## Example

```python
from Modules.Analysis_Figen_block.common import (
    allocate_versioned_group,
    create_run_context,
    update_manifest,
)

context = create_run_context(
    ("analysis", "interaction"),
    "PPS_PR_hotspot_20260822_201530",
)

png_path, csv_path = allocate_versioned_group(
    context.run_dir,
    ("residue_hotspot_full_last50.png", "residue_hotspot_full_last50_data.csv"),
)

# Generate png_path and csv_path, then record them.
update_manifest(
    context,
    script_id="PY-012",
    status="completed",
    outputs=(png_path, csv_path),
)
```

Run the unit tests from `work_space`:

```text
python -m unittest tests.test_common_infrastructure -v
```
