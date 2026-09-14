"""Figure adapters for full/last-window and top-candidate hotspot tables."""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import pandas as pd

from ...interaction_result_analysis.blocks.residue_hotspot import read_pocket_residues


SCRIPT_ID = "PY-133"


def _legacy_module(module_name: str, filename: str):
    directory = Path(__file__).resolve().parents[3] / "ACHresultspackage" / "CompletedAnalysis" / "legacy" / "Fig_gen" / "reshotstopt"
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
    spec = importlib.util.spec_from_file_location(module_name, directory / filename)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load retained renderer: {filename}")
    module = importlib.util.module_from_spec(spec)
    # The retained full-style renderer imports the old calculation CLI even
    # though draw_figure does not use it. Keep that unrelated matplotlib/RDKit
    # dependency outside this Figure-only boundary.
    dependency_name = "generate_crossdock_own_opposite_hotspots"
    previous = sys.modules.get(dependency_name)
    if filename == "generate_crossdock_own_opposite_hotspots_fullstyle.py" and previous is None:
        sys.modules[dependency_name] = types.ModuleType(dependency_name)
    try:
        spec.loader.exec_module(module)
    finally:
        if filename == "generate_crossdock_own_opposite_hotspots_fullstyle.py" and previous is None:
            sys.modules.pop(dependency_name, None)
    return module


def _inputs(full_csv, last_csv, typed_csv, pocket_files, reference_csv):
    full = pd.read_csv(full_csv)
    last = pd.read_csv(last_csv)
    typed = pd.read_csv(typed_csv)
    pockets = {name: read_pocket_residues(path) for name, path in pocket_files.items()}
    reference = pd.read_csv(reference_csv, index_col="residue")
    return full, last, typed, pockets, reference


def render_full_last(*, full_csv, last_csv, typed_csv, pocket_files,
                     reference_csv, denominators, output_png):
    """Render PY-015-compatible output from explicit, stable inputs."""
    full, last, typed, pockets, reference = _inputs(
        full_csv, last_csv, typed_csv, pocket_files, reference_csv
    )
    output_png = Path(output_png)
    renderer = _legacy_module("retained_full_last_hotspot", "generate_residue_hotspot_figure_full_vs_last50.py")
    renderer.OUTPUT_DIR = output_png.parent
    renderer.OUTPUT_STEM = output_png.stem
    table = renderer.draw_figure(full, last, typed, pockets, reference, denominators)
    return output_png, output_png.with_name(output_png.stem + "_data.csv"), table


def render_top_candidates(*, own_csv, opposite_csv, typed_csv, pocket_files,
                          reference_csv, denominators, output_png,
                          title="Top-5 state-selective candidates: own vs opposite receptor state"):
    """Render PY-016-compatible output from explicit, stable inputs."""
    own, opposite, typed, pockets, reference = _inputs(
        own_csv, opposite_csv, typed_csv, pocket_files, reference_csv
    )
    output_png = Path(output_png)
    renderer = _legacy_module("retained_top_candidate_hotspot", "generate_crossdock_own_opposite_hotspots_fullstyle.py")
    renderer.OUTPUT_STEM = output_png.stem
    table = renderer.draw_figure(
        own, opposite, typed, pockets, reference, denominators, output_png.parent, title
    )
    return output_png, output_png.with_name(output_png.stem + "_data.csv"), table
