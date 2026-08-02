import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent / ".matplotlib"),
)

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DATA_DIR = Path(r"Q:\coding_dir\701_Project\Resultsbin\redockingvalidation")
OUTPUT_DIR = Path(r"Q:\coding_dir\701_Project\Writeup_Figures")

COLORS = {
    "PR": "#2878B5",
    "PPS": "#D9534F",
}


def load_dataset(name: str) -> pd.DataFrame:
    docking = pd.read_csv(DATA_DIR / f"{name}_docking.csv")
    rmsd = pd.read_csv(DATA_DIR / f"{name}_RMSD.csv")

    rmsd.columns = (
        rmsd.columns.str.strip().str.replace('"', "", regex=False).str.strip()
    )
    rmsd_values = pd.to_numeric(
        rmsd["RMS"].astype(str).str.replace('"', "", regex=False).str.strip(),
        errors="raise",
    )

    if len(docking) != len(rmsd_values):
        raise ValueError(
            f"{name}: docking ({len(docking)}) and RMSD "
            f"({len(rmsd_values)}) row counts differ."
        )

    return pd.DataFrame(
        {
            "docking_score": pd.to_numeric(docking["docking score"]),
            "rmsd": rmsd_values,
        }
    )


def style_axis(ax, ylabel: str) -> None:
    ax.set_xlabel("Docking pose")
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)
    ax.set_axisbelow(True)


def plot_docking_scores(
    pps: pd.DataFrame, pr: pd.DataFrame, score_ylim: tuple[float, float]
) -> plt.Figure:
    x = np.arange(1, len(pps) + 1)
    fig, ax = plt.subplots(figsize=(11, 5.8))

    # PPS is drawn first and wider, so it stays behind PR.
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
    pps: pd.DataFrame, pr: pd.DataFrame, rmsd_ylim: tuple[float, float]
) -> plt.Figure:
    x = np.arange(1, len(pps) + 1)
    fig, ax = plt.subplots(figsize=(11, 5.8))

    # PPS is drawn first so the PR series remains in the foreground.
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
    name: str,
    score_ylim: tuple[float, float],
    rmsd_ylim: tuple[float, float],
) -> plt.Figure:
    x = np.arange(1, len(data) + 1)
    color = COLORS[name]
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
        label=f"{name} RMSD",
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
        label=f"{name} docking score",
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

    for ax in (ax_rmsd, ax_score):
        ax.set_xlim(0.3, len(x) + 0.7)
    ax_score.set_xticks(np.arange(1, len(x) + 1, 2))

    fig.suptitle(f"{name} redocking validation", y=0.98)
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.10, top=0.92)
    return fig


def save_figure(fig: plt.Figure, stem: str) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    pps = load_dataset("PPS")
    pr = load_dataset("PR")

    if len(pps) != len(pr):
        raise ValueError(f"PPS ({len(pps)}) and PR ({len(pr)}) row counts differ.")

    score_ylim = (-12, -1)
    rmsd_ylim = (0, 2.5)

    save_figure(
        plot_docking_scores(pps, pr, score_ylim),
        "redocking_docking_score_comparison",
    )
    save_figure(
        plot_rmsd(pps, pr, rmsd_ylim),
        "redocking_rmsd_comparison",
    )
    save_figure(
        plot_mirrored_dataset(pr, "PR", score_ylim, rmsd_ylim),
        "redocking_PR_mirrored",
    )
    save_figure(
        plot_mirrored_dataset(pps, "PPS", score_ylim, rmsd_ylim),
        "redocking_PPS_mirrored",
    )


if __name__ == "__main__":
    main()