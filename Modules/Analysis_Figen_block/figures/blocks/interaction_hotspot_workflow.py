"""Versioned workflow for the base PR/PPS interaction hotspot figure."""

from __future__ import annotations

from pathlib import Path

from ...common.output_naming import allocate_versioned_group
from ...common.provenance import update_manifest
from ...common.run_context import create_run_context
from .interaction_hotspot import SCRIPT_ID, render_interaction_hotspot


def run_interaction_hotspot(*, residue_csv, typed_csv, pocket_files,
                            reference_csv, run_id, results_root=None,
                            resume=False, upstream_run_ids=()):
    context = create_run_context(
        ("Figs", "interaction_hotspot"), run_id,
        results_root=results_root, resume=resume,
    )
    png, data = allocate_versioned_group(
        context.run_dir,
        ("interaction_hotspot.png", "interaction_hotspot_data.csv"),
    )
    inputs = (
        Path(residue_csv), Path(typed_csv), Path(reference_csv),
        *(Path(path) for path in pocket_files.values()),
    )
    try:
        _, generated_data, table = render_interaction_hotspot(
            residue_csv=residue_csv, typed_csv=typed_csv,
            pocket_files=pocket_files, reference_csv=reference_csv,
            output_png=png,
        )
        if generated_data != data:
            generated_data.replace(data)
        update_manifest(
            context, script_id=SCRIPT_ID, status="completed", config={},
            inputs=inputs, outputs=(png, data),
            upstream_run_ids=upstream_run_ids,
        )
    except Exception:
        update_manifest(
            context, script_id=SCRIPT_ID, status="failed", config={},
            upstream_run_ids=upstream_run_ids,
        )
        raise
    return context, png, data, table
