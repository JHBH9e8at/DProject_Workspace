"""Adapter from stable own/opposite tables to the retained full-style renderer."""

from __future__ import annotations
import importlib.util, sys
from pathlib import Path
import pandas as pd
from ...interaction_result_analysis.blocks.own_opposite_hotspot import ARM_SPECS
from ...interaction_result_analysis.blocks.residue_hotspot import read_pocket_residues

SCRIPT_ID = "PY-128"

def prepare_renderer_tables(prevalence_csv, typed_csv):
    summary = pd.read_csv(prevalence_csv); typed = pd.read_csv(typed_csv)
    mapping = {arm: spec["ligand"] for arm, spec in ARM_SPECS.items()}
    summary["population_state"] = summary.crossdock_arm.map(mapping)
    summary["poses"] = summary.molecule_count
    return summary[summary.condition.eq("own")].copy(), summary[summary.condition.eq("opposite")].copy(), typed

def _legacy_renderer():
    directory = Path(__file__).resolve().parents[3] / "ACHresultspackage" / "CompletedAnalysis" / "legacy" / "Fig_gen" / "reshotstopt"
    sys.path.insert(0, str(directory)) if str(directory) not in sys.path else None
    path = directory / "generate_crossdock_own_opposite_hotspots_fullstyle.py"
    spec = importlib.util.spec_from_file_location("retained_fullstyle_hotspot", path)
    module = importlib.util.module_from_spec(spec); assert spec and spec.loader; spec.loader.exec_module(module)
    return module

def render_fullstyle(*, prevalence_csv, typed_csv, pocket_files: dict[str, str | Path],
                     reference_csv, denominators: dict[str, int], output_png: str | Path,
                     title="PR/PPS residue interaction hotspots: own vs opposite receptor state"):
    own, opposite, typed = prepare_renderer_tables(prevalence_csv, typed_csv)
    pockets = {name: read_pocket_residues(path) for name, path in pocket_files.items()}
    reference = pd.read_csv(reference_csv, index_col="residue")
    renderer = _legacy_renderer(); output_png = Path(output_png)
    renderer.OUTPUT_STEM = output_png.stem
    table = renderer.draw_figure(own, opposite, typed, pockets, reference,
                                 denominators, output_png.parent, title)
    return output_png, output_png.with_name(output_png.stem + "_data.csv"), table
