"""Render joint chemical-space figures from stable embedding cache tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
import pandas as pd


SCRIPT_ID = "PY-112"
COLORS = {"PR": "#2878B5", "PPS": "#D9534F"}
REF_COLOR = "#7F3C8D"
STEP_CMAPS = {key: LinearSegmentedColormap.from_list(f"{key}_steps", ["#F2F2F2", color])
              for key, color in COLORS.items()}
REQUIRED = {"source", "source_index", "smiles", "canon_smiles", "step", "docking_score",
            "label", "embedding_x", "embedding_y", "method", "feature"}


def load_embedding_cache(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    missing = sorted(REQUIRED.difference(frame.columns))
    if missing: raise ValueError(f"Embedding cache missing columns: {missing}")
    if frame.empty: raise ValueError("Embedding cache is empty")
    if frame["source"].eq("REF").sum() == 0: raise ValueError("Embedding cache contains no reference")
    for column in ("step", "docking_score", "embedding_x", "embedding_y"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if frame[["embedding_x", "embedding_y"]].isna().any().any():
        raise ValueError("Embedding coordinates contain non-numeric values")
    return frame


def _base(frame: pd.DataFrame):
    fig, axes = plt.subplots(1, 3, figsize=(19, 7), sharex=True, sharey=True,
                             constrained_layout=True)
    population = frame[frame["source"].isin(["PR", "PPS"])]
    references = frame[frame["source"] == "REF"]
    for source in ("PR", "PPS"):
        part = population[population["source"] == source]
        axes[0].scatter(part["embedding_x"], part["embedding_y"], s=5, alpha=.35,
                        color=COLORS[source], linewidths=0, rasterized=True, label=source)
    for _, ref in references.iterrows():
        label = ref.get("label", "Reference")
        for ax in axes:
            ax.scatter(ref["embedding_x"], ref["embedding_y"], marker="*", s=190,
                       color=REF_COLOR, edgecolor="white", linewidth=.9, zorder=10,
                       label=label if ax is axes[0] else "_nolegend_")
            ax.annotate(label, (ref["embedding_x"], ref["embedding_y"]), xytext=(7, 6),
                        textcoords="offset points", fontsize=10, fontweight="bold", zorder=11)
    axes[0].set_title("PR + PPS", fontsize=16, fontweight="normal")
    axes[0].legend(loc="upper center", bbox_to_anchor=(.5, -.15), ncol=3,
                   frameon=False, fontsize=12, markerscale=2.5)
    return fig, axes, population


def render_iteration(cache_path: str | Path, output_path: str | Path) -> Path:
    frame = load_embedding_cache(cache_path); fig, axes, population = _base(frame)
    method = str(frame["method"].iloc[0]).upper(); feature = frame["feature"].iloc[0]
    population = population.dropna(subset=["step"]); norm = Normalize(population.step.min(), population.step.max())
    for ax, source in zip(axes[1:], ("PR", "PPS")):
        part = population[population.source == source].sort_values("step")
        ax.scatter(part.embedding_x, part.embedding_y, c=part.step, cmap=STEP_CMAPS[source],
                   norm=norm, s=6, alpha=.55, linewidths=0, rasterized=True)
        ax.set_title(source, fontsize=16, fontweight="normal")
        mappable = ScalarMappable(norm=norm, cmap=STEP_CMAPS[source]); mappable.set_array([])
        bar = fig.colorbar(mappable, ax=ax, orientation="horizontal", location="bottom", shrink=.82, pad=-.1)
        bar.set_label("Iteration step", labelpad=8)
    axes[0].set_xlabel(f"{method} 1"); axes[0].set_ylabel(f"{method} 2")
    for ax in axes: ax.grid(alpha=.1)
    fig.suptitle(f"Joint PR/PPS Chemical Space by Iteration — {method} ({feature})", fontsize=15)
    output = Path(output_path); output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight"); plt.close(fig); return output


def render_docking(cache_path: str | Path, output_path: str | Path) -> Path:
    frame = load_embedding_cache(cache_path); fig, axes, population = _base(frame)
    method = str(frame["method"].iloc[0]).upper(); feature = frame["feature"].iloc[0]
    population = population.dropna(subset=["docking_score"])
    norm = Normalize(population.docking_score.min(), population.docking_score.max())
    for ax, source in zip(axes[1:], ("PR", "PPS")):
        part = population[population.source == source].sort_values("docking_score", ascending=False)
        ax.scatter(part.embedding_x, part.embedding_y, c=part.docking_score, cmap="viridis",
                   norm=norm, s=6, alpha=.55, linewidths=0, rasterized=True)
        ax.set_title(source, fontsize=16, fontweight="normal")
        mappable = ScalarMappable(norm=norm, cmap="viridis"); mappable.set_array([])
        bar = fig.colorbar(mappable, ax=ax, orientation="horizontal", location="bottom", shrink=.82, pad=-.1)
        bar.set_label("Docking score", labelpad=8)
    axes[0].set_xlabel(f"{method} 1"); axes[0].set_ylabel(f"{method} 2")
    for ax in axes: ax.grid(alpha=.1)
    fig.suptitle(f"Joint PR/PPS Chemical Space by Docking Score — {method} ({feature})", fontsize=15)
    output = Path(output_path); output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight"); plt.close(fig); return output
