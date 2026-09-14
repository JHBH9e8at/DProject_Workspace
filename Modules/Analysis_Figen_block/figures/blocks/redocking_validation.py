"""Render redocking validation figures from a stable pose-level table."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


COLORS = {"PR": "#2878B5", "PPS": "#D9534F"}


def split_states(poses: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"state", "pose_index", "docking_score", "rmsd"}
    missing = required - set(poses.columns)
    if missing:
        raise ValueError(f"Pose table is missing columns: {sorted(missing)}")
    pps = poses.loc[poses["state"] == "PPS"].sort_values("pose_index")
    pr = poses.loc[poses["state"] == "PR"].sort_values("pose_index")
    if pps.empty or pr.empty or len(pps) != len(pr):
        raise ValueError(
            f"Pose table requires equal non-empty PPS/PR rows; "
            f"received PPS={len(pps)}, PR={len(pr)}"
        )
    return pps.reset_index(drop=True), pr.reset_index(drop=True)


def style_axis(ax, ylabel: str) -> None:
    ax.set_xlabel("Docking pose")
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)
    ax.set_axisbelow(True)


def plot_docking_scores(
    pps: pd.DataFrame,
    pr: pd.DataFrame,
    score_ylim: tuple[float, float],
) -> plt.Figure:
    x = np.arange(1, len(pps) + 1)
    fig, ax = plt.subplots(figsize=(11, 5.8))
    ax.bar(
        x,
        pps["docking_score"],
        width=0.82,
        color=COLORS["PPS"],
        edgecolor="none",
        label="PPS",
        zorder=2,
    )
    ax.bar(
        x,
        pr["docking_score"],
        width=0.50,
        color=COLORS["PR"],
        edgecolor="none",
        label="PR",
        zorder=3,
    )
    style_axis(ax, "Docking score (kcal/mol)")
    ax.set_title("Redocking validation: Docking score")
    ax.set_xlim(0.3, len(x) + 0.7)
    ax.set_ylim(*score_ylim)
    ax.set_xticks(np.arange(1, len(x) + 1, 2))
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    return fig


def plot_rmsd(
    pps: pd.DataFrame,
    pr: pd.DataFrame,
    rmsd_ylim: tuple[float, float],
) -> plt.Figure:
    x = np.arange(1, len(pps) + 1)
    fig, ax = plt.subplots(figsize=(11, 5.8))
    ax.plot(
        x,
        pps["rmsd"],
        color=COLORS["PPS"],
        marker="o",
        markersize=4.5,
        linewidth=2,
        label="PPS",
        zorder=2,
    )
    ax.plot(
        x,
        pr["rmsd"],
        color=COLORS["PR"],
        marker="s",
        markersize=4.5,
        linewidth=2,
        label="PR",
        zorder=3,
    )
    ax.axhline(
        2.0,
        color="#4D4D4D",
        linestyle="--",
        linewidth=1.2,
        label="RMSD cutoff (2 Å)",
        zorder=1,
    )
    style_axis(ax, "RMSD (Å)")
    ax.set_title("Redocking validation: RMSD")
    ax.set_xlim(0.3, len(x) + 0.7)
    ax.set_ylim(*rmsd_ylim)
    ax.set_xticks(np.arange(1, len(x) + 1, 2))
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    return fig


def plot_mirrored_dataset(
    data: pd.DataFrame,
    state: str,
    score_ylim: tuple[float, float],
    rmsd_ylim: tuple[float, float],
) -> plt.Figure:
    x = np.arange(1, len(data) + 1)
    color = COLORS[state]
    fig, (ax_rmsd, ax_score) = plt.subplots(
        2,
        1,
        figsize=(11, 7.2),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1.35], "hspace": 0.04},
    )
    ax_rmsd.plot(
        x,
        data["rmsd"],
        color=color,
        marker="o",
        markersize=4.5,
        linewidth=2,
        label=f"{state} RMSD",
        zorder=3,
    )
    ax_rmsd.axhline(
        2.0,
        color="#4D4D4D",
        linestyle="--",
        linewidth=1.2,
        label="RMSD cutoff (2 Å)",
        zorder=2,
    )
    ax_rmsd.set_ylabel("RMSD (Å)")
    ax_rmsd.set_ylim(*rmsd_ylim)
    ax_rmsd.spines[["top", "right"]].set_visible(False)
    ax_rmsd.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)
    ax_rmsd.set_axisbelow(True)
    ax_rmsd.legend(frameon=False, ncol=2, loc="upper right")
    ax_score.bar(
        x,
        data["docking_score"],
        width=0.72,
        color=color,
        edgecolor="none",
        label=f"{state} docking score",
        zorder=3,
    )
    ax_score.set_ylabel("Docking score (kcal/mol)")
    ax_score.set_xlabel("Docking pose")
    ax_score.set_ylim(*score_ylim)
    ax_score.spines[["right", "bottom"]].set_visible(False)
    ax_score.spines["top"].set_color("#4D4D4D")
    ax_score.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)
    ax_score.set_axisbelow(True)
    ax_score.legend(frameon=False, loc="lower right")
    for axis in (ax_rmsd, ax_score):
        axis.set_xlim(0.3, len(x) + 0.7)
    ax_score.set_xticks(np.arange(1, len(x) + 1, 2))
    fig.suptitle(f"{state} redocking validation", y=0.98)
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.10, top=0.92)
    return fig


def render_figure_set(
    poses: pd.DataFrame,
    output_paths: dict[str, tuple[Path, Path]],
    *,
    score_ylim: tuple[float, float] = (-12, -1),
    rmsd_ylim: tuple[float, float] = (0, 2.5),
) -> list[Path]:
    """Render four PNG/SVG pairs to preallocated non-overwriting paths."""

    pps, pr = split_states(poses)
    figures = {
        "redocking_docking_score_comparison": plot_docking_scores(
            pps, pr, score_ylim
        ),
        "redocking_rmsd_comparison": plot_rmsd(pps, pr, rmsd_ylim),
        "redocking_PR_mirrored": plot_mirrored_dataset(
            pr, "PR", score_ylim, rmsd_ylim
        ),
        "redocking_PPS_mirrored": plot_mirrored_dataset(
            pps, "PPS", score_ylim, rmsd_ylim
        ),
    }
    written: list[Path] = []
    try:
        for stem, figure in figures.items():
            png_path, svg_path = output_paths[stem]
            figure.savefig(png_path, dpi=300, bbox_inches="tight")
            figure.savefig(svg_path, bbox_inches="tight")
            written.extend((png_path, svg_path))
    finally:
        for figure in figures.values():
            plt.close(figure)
    return written
