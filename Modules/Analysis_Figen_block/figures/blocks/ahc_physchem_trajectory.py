"""Render legacy-equivalent AHC physicochemical and trajectory figures from tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from ...ahc_analysis.blocks.physchem_trajectory import DESCRIPTORS


SCRIPT_ID = "PY-106"
REF_MAP = {"OM": "OmecamtivMecarbil", "MAV": "Mavacamten", "AFI": "Aficamten"}
REF_COLORS = {"OM": "#7F3C8D", "MAV": "#E6A700", "AFI": "#2CA02C"}
DISPLAY_NAMES = {
    "desc_CLogP": "CLogP",
    "desc_MolWt": "Molecular weight",
    "desc_NumRotatableBonds": "Rotatable bonds",
}
SIGNATURE_COLORS = {
    "PPS": {"full": "#F2AAA6", "last": "#D9534F"},
    "PR": {"full": "#9CC4E4", "last": "#2878B5"},
}


def _make_property_axes(descriptors: list[str]):
    if descriptors == DESCRIPTORS:
        fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        flat = list(axes.flatten())
        flat[-1].set_visible(False)
        return fig, flat
    fig, axes = plt.subplots(1, len(descriptors), figsize=(6 * len(descriptors), 5), squeeze=False)
    return fig, list(axes.flatten())


def _apply_property_x_axis(ax, column: str, mw_upper_bound: float) -> None:
    if column == "desc_MolWt":
        ax.set_xlim(0, mw_upper_bound)
    elif column == "desc_NumRotatableBonds":
        ax.set_xlim(0, 17)
        ax.set_xticks(range(18))


def load_references(ref_csv: str | Path | None, ref_keys: list[str]) -> dict[str, pd.Series]:
    if not ref_keys:
        return {}
    if ref_csv is None:
        raise ValueError("ref_csv is required when reference keys are selected")
    frame = pd.read_csv(ref_csv)
    references = {}
    for key in ref_keys:
        rows = frame.loc[frame["name"] == REF_MAP[key]]
        if rows.empty:
            raise ValueError(f"Reference not found: {REF_MAP[key]}")
        references[key] = rows.iloc[0]
    return references


def render_figures(*, rows_csv: str | Path, trajectory_csv: str | Path,
                   output_paths: dict[str, Path], job_name: str, threshold: float,
                   ref_csv: str | Path | None = None, ref_keys: list[str] | None = None,
                   descriptors: list[str] | None = None,
                   mw_upper_bound: float | None = None) -> None:
    rows = pd.read_csv(rows_csv)
    trajectory = pd.read_csv(trajectory_csv)
    refs = load_references(ref_csv, ref_keys or [])
    score_col = f"{job_name}_r_i_docking_score"
    valid_rows = rows.loc[rows[score_col] != 0].copy()
    if valid_rows.empty:
        raise ValueError(f"{job_name}: no non-zero docking scores are available")
    mw_upper = (
        float(mw_upper_bound)
        if mw_upper_bound is not None
        else float(valid_rows["desc_MolWt"].max()) + 100.0
    )
    if mw_upper <= 0:
        raise ValueError("mw_upper_bound must be greater than zero")
    selected = descriptors or DESCRIPTORS
    legacy_layout = selected == DESCRIPTORS
    signature = SIGNATURE_COLORS.get(job_name.upper(),
                                     {"full": "steelblue", "last": "tomato"})

    fig, axes = _make_property_axes(selected)
    for index, column in enumerate(selected):
        ax = axes[index]
        ax.hist(valid_rows[column].dropna(), bins=50, color=signature["last"], edgecolor="none", alpha=0.7)
        for key, ref in refs.items():
            if column in ref and pd.notna(ref[column]):
                ax.axvline(ref[column], color=REF_COLORS[key], linewidth=1.5,
                           linestyle="--", label=REF_MAP[key])
        label = column if legacy_layout else DISPLAY_NAMES.get(column, column)
        ax.set_title(label, fontsize=12); ax.set_xlabel(""); ax.set_ylabel("Count")
        _apply_property_x_axis(ax, column, mw_upper)
    if refs:
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.01),
                   ncol=max(1, len(labels)), fontsize=11, frameon=False)
    plt.suptitle(f"[{job_name}] Physicochemical Properties Distribution", fontsize=16, y=1.02)
    plt.tight_layout(rect=(0, 0.10 if refs else 0, 1, 1)); fig.savefig(output_paths["distribution_png"], dpi=150, bbox_inches="tight")
    fig.savefig(output_paths["distribution_svg"], bbox_inches="tight"); plt.close(fig)

    fig, axes = _make_property_axes(selected)
    max_step = valid_rows["step"].max(); min_last50 = max_step - 50
    for index, column in enumerate(selected):
        ax = axes[index]; sample = valid_rows[[column, score_col, "step"]].dropna()
        recent = sample["step"] >= min_last50; good = sample[score_col] <= threshold; weak = ~good
        for data, color, label in ((sample[~recent], signature["full"], f"steps 0 ~ {min_last50 - 1}"),
                                   (sample[recent], signature["last"], f"steps {min_last50} ~ {max_step}")):
            ax.scatter(data.loc[good.reindex(data.index, fill_value=False), column],
                       data.loc[good.reindex(data.index, fill_value=False), score_col],
                       alpha=0.5, s=5, color=color, label=label)
            ax.scatter(data.loc[weak.reindex(data.index, fill_value=False), column],
                       data.loc[weak.reindex(data.index, fill_value=False), score_col],
                       alpha=0.05, s=5, color=color)
        ax.axhline(threshold, color="black", linewidth=1, linestyle="--", alpha=0.7,
                   label=f"docking = {threshold}")
        label = column if legacy_layout else DISPLAY_NAMES.get(column, column)
        ax.set_title(label, fontsize=11); ax.set_xlabel(label); ax.set_ylabel(score_col)
        _apply_property_x_axis(ax, column, mw_upper)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.01),
               ncol=len(labels), fontsize=11, frameon=False)
    plt.suptitle(f"[{job_name}] Physicochemical Properties vs Docking Score", fontsize=16, y=1.02)
    plt.tight_layout(rect=(0, 0.12, 1, 1)); fig.savefig(output_paths["scatter_png"], dpi=150, bbox_inches="tight")
    fig.savefig(output_paths["scatter_svg"], bbox_inches="tight"); plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10),
        gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    plotted = trajectory.dropna(subset=["mean"])
    ax1.plot(plotted["step"], plotted["mean"], color=signature["last"], linewidth=1.5, label="mean")
    ax1.fill_between(plotted["step"], plotted["mean"] + plotted["std"],
                     plotted["mean"] - plotted["std"], color=signature["full"],
                     alpha=0.5, label="±1 std")
    ax1.set_xlabel("Step"); ax1.tick_params(labelbottom=True); ax1.set_ylabel("Mean Doccking score")
    failed = trajectory.loc[trajectory["failure_count"] > 0]
    ax2.bar(failed["step"], failed["failure_count"], color=signature["last"],
            alpha=0.6, label="zero count")
    ax2.set_ylabel("Zero count"); ax2.set_xlabel("Step"); ax2.legend()
    lines1, labels1 = ax1.get_legend_handles_labels(); lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2); fig.suptitle("Docking Score Trajectory", fontsize=14)
    plt.tight_layout(); fig.savefig(output_paths["trajectory_png"], dpi=150, bbox_inches="tight")
    fig.savefig(output_paths["trajectory_svg"], bbox_inches="tight"); plt.close(fig)
