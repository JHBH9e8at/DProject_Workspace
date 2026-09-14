"""Adapter from stable last-window tables to the retained PY-009 layout."""

import importlib.util
from pathlib import Path
import pandas as pd


SCRIPT_ID = "PY-148"


def _legacy_renderer():
    path = Path(__file__).resolve().parents[3] / "ACHresultspackage" / "CompletedAnalysis" / "legacy" / "Fig_gen" / "physchem_scatter_last_iterations.py"
    spec = importlib.util.spec_from_file_location("retained_physchem_last_window", path)
    module = importlib.util.module_from_spec(spec); assert spec and spec.loader
    spec.loader.exec_module(module); return module


def render_last_window(*, rows_csv, summary_csv, output_png, last_n,
                       pr_threshold=None, pps_threshold=None):
    rows = pd.read_csv(rows_csv); summaries = pd.read_csv(summary_csv).set_index("population_state")
    renderer = _legacy_renderer()

    def prepared(_input, job_name, requested_last_n):
        if int(requested_last_n) != int(last_n):
            raise ValueError("Renderer last_n does not match stable table")
        part = rows[rows.population_state.eq(job_name)].copy()
        summary = summaries.loc[job_name]
        score_col = str(summary.score_column)
        return {"jn": job_name, "score_col": score_col,
                "earlier": part[part.window.eq("earlier")],
                "recent": part[part.window.eq("recent")],
                "first_last_step": summary.first_last_step,
                "max_step": summary.max_step}

    renderer.prepare_data = prepared
    output_png = Path(output_png); output_png.parent.mkdir(parents=True, exist_ok=True)
    generated = Path(renderer.plot_pr_pps_comparison(
        "stable-table", "stable-table", str(output_png.parent), last_n=int(last_n),
        pr_threshold=pr_threshold, pps_threshold=pps_threshold))
    if generated != output_png:
        generated.replace(output_png)
    return output_png
