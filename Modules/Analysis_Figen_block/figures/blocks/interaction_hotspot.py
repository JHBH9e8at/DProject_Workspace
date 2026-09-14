"""Render the base PR/PPS interaction hotspot figure from stable summary tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from PIL import ImageFont

from ...interaction_result_analysis.blocks.residue_hotspot import read_pocket_residues
from .hotspot_variants import _legacy_module


SCRIPT_ID = "PY-141"


def _portable_truetype(original):
    regular = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
    bold = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf")

    def load(font, size, *args, **kwargs):
        requested = Path(str(font))
        if requested.is_file():
            return original(str(requested), size, *args, **kwargs)
        fallback = bold if "bd" in requested.stem.lower() else regular
        if fallback.is_file():
            return original(str(fallback), size, *args, **kwargs)
        return ImageFont.load_default(size=size)

    return load


def render_interaction_hotspot(*, residue_csv, typed_csv, pocket_files,
                               reference_csv, output_png):
    residue = pd.read_csv(residue_csv)
    typed = pd.read_csv(typed_csv)
    pockets = {
        name: read_pocket_residues(path) for name, path in pocket_files.items()
    }
    reference = pd.read_csv(reference_csv, index_col="residue")
    output_png = Path(output_png)
    renderer = _legacy_module(
        "retained_base_interaction_hotspot", "generate_residue_hotspot_figure.py"
    )
    renderer.OUTPUT_DIR = output_png.parent
    original_truetype = renderer.ImageFont.truetype
    renderer.ImageFont.truetype = _portable_truetype(original_truetype)
    try:
        table = renderer.build_figure(residue, typed, pockets, reference)
    finally:
        renderer.ImageFont.truetype = original_truetype

    generated_png = output_png.parent / "pr_pps_residue_hotspot_pocket_map.png"
    generated_data = output_png.parent / "pr_pps_residue_hotspot_pocket_map_data.csv"
    target_data = output_png.with_name(output_png.stem + "_data.csv")
    if generated_png != output_png:
        generated_png.replace(output_png)
    if generated_data != target_data:
        generated_data.replace(target_data)
    return output_png, target_data, table
