from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


COLORS = {"PR": "#2878B5", "PPS": "#D9534F"}
DEFAULT_RESULTS = Path(r"Q:\coding_dir\701_Project\Resultsbin\CDRes_LP_On")
SCORE_COLUMN = "r_i_docking_score"


def gaussian_kde(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Evaluate a Gaussian KDE using Silverman's bandwidth rule."""
    count = values.size
    standard_deviation = values.std(ddof=1)
    bandwidth = 1.06 * standard_deviation * count ** (-1 / 5)
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError("KDE requires scores with non-zero variance.")
    scaled = (grid[:, None] - values[None, :]) / bandwidth
    kernels = np.exp(-0.5 * scaled**2) / np.sqrt(2 * np.pi)
    return kernels.mean(axis=1) / bandwidth


def parse_config(path: Path) -> dict[str, str]:
    config: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, value = line.split("=", 1)
                config[key.strip().lower().replace("-", "_")] = value.strip()
    for key in ("results_dir", "output_dir"):
        if not config.get(key):
            raise ValueError(f"Missing required config key: {key}")
    config.setdefault("input_mode", "own-state")
    config.setdefault("pps_color", COLORS["PPS"])
    config.setdefault("pr_color", COLORS["PR"])
    config.setdefault("opposite_color", "#9A9A9A")
    config.setdefault("figure_width", "7.2")
    config.setdefault("figure_height", "4.8")
    config.setdefault("dpi", "300")
    config.setdefault("fill_alpha", "0.22")
    config.setdefault("line_width", "2.2")
    config.setdefault("show_median", "on")
    config.setdefault("show_median_label", "on")
    config.setdefault("median_label_x_offset", "8")
    config.setdefault("median_label_y_offset", "22")
    config.setdefault("bottom_median_label_x_offset", "8")
    config.setdefault("bottom_median_label_y_offset", "22")
    config.setdefault("bottom_left_median_label_x_offset",
                      config["bottom_median_label_x_offset"])
    config.setdefault("bottom_left_median_label_y_offset",
                      config["bottom_median_label_y_offset"])
    config.setdefault("bottom_right_median_label_x_offset",
                      config["bottom_median_label_x_offset"])
    config.setdefault("bottom_right_median_label_y_offset",
                      config["bottom_median_label_y_offset"])
    config.setdefault("bottom_left_own_median_label_x_offset",
                      config["bottom_left_median_label_x_offset"])
    config.setdefault("bottom_left_own_median_label_y_offset",
                      config["bottom_left_median_label_y_offset"])
    config.setdefault("bottom_right_own_median_label_x_offset",
                      config["bottom_right_median_label_x_offset"])
    config.setdefault("bottom_right_own_median_label_y_offset",
                      config["bottom_right_median_label_y_offset"])
    config.setdefault("composite_figure_width", "10.0")
    config.setdefault("composite_figure_height", "8.0")
    config.setdefault("panel_wspace", "0.12")
    config.setdefault("panel_hspace", "0.12")
    config.setdefault("bottom_legend_loc", "upper right")
    config.setdefault("bottom_left_legend_loc", config["bottom_legend_loc"])
    config.setdefault("bottom_right_legend_loc", "upper left")
    config.setdefault("legend_font_size", "9")
    return config


def resolve_path(value: str, config_path: Path) -> Path:
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    if not path.is_absolute():
        path = config_path.resolve().parent / path
    return path


def load_scores(results_dir: Path, arm: str) -> np.ndarray:
    path = results_dir / "crossdocking" / arm / "best_docking_scores.csv"
    table = pd.read_csv(path, usecols=[SCORE_COLUMN])
    scores = pd.to_numeric(table[SCORE_COLUMN], errors="coerce").dropna().to_numpy()
    if scores.size < 2:
        raise ValueError(f"At least two valid scores are required: {path}")
    return scores


def load_double_arm_scores(results_dir: Path, arm: str) -> np.ndarray:
    path = results_dir / "crossdocking" / "double_arm_docking_scores.csv"
    table = pd.read_csv(path, usecols=["crossdock_arm", SCORE_COLUMN])
    scores = pd.to_numeric(
        table.loc[table["crossdock_arm"] == arm, SCORE_COLUMN], errors="coerce"
    ).dropna().to_numpy()
    if scores.size < 2:
        raise ValueError(f"At least two valid scores are required for {arm}: {path}")
    return scores


def load_own_state_scores(results_dir: Path, arm: str) -> np.ndarray:
    path = results_dir / "analysis" / "selectivity" / "paired_state_scores.csv"
    table = pd.read_csv(path, usecols=["crossdock_arm", "own_state_score"])
    scores = pd.to_numeric(
        table.loc[table["crossdock_arm"] == arm, "own_state_score"], errors="coerce"
    ).dropna().to_numpy()
    if scores.size < 2:
        raise ValueError(f"At least two valid own-state scores are required for {arm}: {path}")
    return scores


def load_paired_state_scores(results_dir: Path) -> pd.DataFrame:
    path = results_dir / "analysis" / "selectivity" / "paired_state_scores.csv"
    table = pd.read_csv(
        path,
        usecols=["crossdock_arm", "own_state_score", "opposite_state_score"],
    )
    for column in ("own_state_score", "opposite_state_score"):
        table[column] = pd.to_numeric(table[column], errors="coerce")
    return table.dropna(subset=["own_state_score", "opposite_state_score"])


def style_axis(ax: plt.Axes, x_label: str) -> None:
    ax.set_xlabel(x_label)
    ax.set_ylabel("Density")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")


def draw_density_axis(ax: plt.Axes, series: list[tuple[str, np.ndarray, str]],
                      config: dict[str, str], title: str | None = None,
                      legend_loc: str = "best",
                      median_x_offset: float | None = None,
                      median_y_offset: float | None = None,
                      median_offsets_by_label: dict[str, tuple[float, float]] | None = None) -> None:
    if median_x_offset is None:
        median_x_offset = float(config["median_label_x_offset"])
    if median_y_offset is None:
        median_y_offset = float(config["median_label_y_offset"])
    median_offsets_by_label = median_offsets_by_label or {}
    all_scores = np.concatenate([scores for _, scores, _ in series])
    padding = max(0.4, 0.06 * np.ptp(all_scores))
    x_grid = np.linspace(all_scores.min() - padding, all_scores.max() + padding, 600)
    maximum_density = 0.0
    legend_handles: list[Patch] = []
    for label, scores, color in series:
        density = gaussian_kde(scores, x_grid)
        maximum_density = max(maximum_density, float(density.max()))
        ax.fill_between(x_grid, density, color=color,
                        alpha=float(config["fill_alpha"]), linewidth=0)
        ax.plot(x_grid, density, color=color,
                linewidth=float(config["line_width"]),
                label="_nolegend_")
        legend_handles.append(
            Patch(
                facecolor=color,
                edgecolor=color,
                alpha=float(config["fill_alpha"]),
                label=f"{label} (n={scores.size})",
            )
        )
        if config["show_median"].lower() == "on":
            median = float(np.median(scores))
            ax.axvline(median, color=color, linewidth=1.4,
                       linestyle="--", alpha=0.9)
            if config["show_median_label"].lower() == "on":
                median_density = float(np.interp(median, x_grid, density))
                label_x_offset, label_y_offset = median_offsets_by_label.get(
                    label, (median_x_offset, median_y_offset)
                )
                ax.annotate(
                    f"Median = {median:.2f}", xy=(median, median_density),
                    xytext=(label_x_offset, label_y_offset),
                    textcoords="offset points", ha="left", va="bottom",
                    color=color, fontsize=8, fontweight="bold",
                )
    if title:
        ax.set_title(title, pad=10)
    ax.set_ylim(0, maximum_density * 1.22)
    ax.legend(handles=legend_handles, frameon=False, loc=legend_loc,
              fontsize=float(config["legend_font_size"]))


def plot_density(results_dir: Path, output_dir: Path, input_mode: str,
                 config: dict[str, str]) -> None:
    colors = {"PPS": config["pps_color"], "PR": config["pr_color"]}
    datasets = [
        ("L_PPS_to_R_PR", "PPS ligand to PR receptor", "PPS"),
        ("L_PR_to_R_PPS", "PR ligand to PPS receptor", "PR"),
    ]
    if input_mode == "own-state":
        datasets = [
            ("L_PPS_to_R_PR", "L_PPS to R_PR", "PPS"),
            ("L_PR_to_R_PPS", "L_PR to R_PPS", "PR"),
        ]
        loader = load_own_state_scores
    elif input_mode == "double-arm":
        loader = load_double_arm_scores
    else:
        loader = load_scores
    loaded = [(arm, label, source, loader(results_dir, arm))
              for arm, label, source in datasets]

    x_label = ("Docking score (kcal/mol)" if input_mode == "double-arm"
               else "Best docking score (kcal/mol)")

    if input_mode == "own-state":
        fig = plt.figure(figsize=(float(config["composite_figure_width"]),
                                  float(config["composite_figure_height"])),
                         layout="constrained")
        grid = fig.add_gridspec(
            2,
            2,
            height_ratios=(1.15, 1),
            hspace=float(config["panel_hspace"]),
            wspace=float(config["panel_wspace"]),
        )
        top_ax = fig.add_subplot(grid[0, :])
        draw_density_axis(
            top_ax,
            [(label, scores, colors[source]) for _, label, source, scores in loaded],
            config,
            "Own-state docking-score distributions",
        )
        style_axis(top_ax, x_label)

        paired = load_paired_state_scores(results_dir)
        bottom_specs = [
            ("L_PPS_to_R_PR", "PPS ligands", "PPS"),
            ("L_PR_to_R_PPS", "PR ligands", "PR"),
        ]
        bottom_axes = (fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1]))
        bottom_legend_locs = (
            config["bottom_left_legend_loc"],
            config["bottom_right_legend_loc"],
        )
        bottom_median_offsets = (
            (float(config["bottom_left_median_label_x_offset"]),
             float(config["bottom_left_median_label_y_offset"])),
            (float(config["bottom_right_median_label_x_offset"]),
             float(config["bottom_right_median_label_y_offset"])),
        )
        bottom_own_median_offsets = (
            (float(config["bottom_left_own_median_label_x_offset"]),
             float(config["bottom_left_own_median_label_y_offset"])),
            (float(config["bottom_right_own_median_label_x_offset"]),
             float(config["bottom_right_own_median_label_y_offset"])),
        )
        for ax, (arm, title, own_color), legend_loc, median_offsets, own_offsets in zip(
                bottom_axes, bottom_specs, bottom_legend_locs,
                bottom_median_offsets, bottom_own_median_offsets):
            arm_table = paired.loc[paired["crossdock_arm"] == arm]
            draw_density_axis(
                ax,
                [("Own state", arm_table["own_state_score"].to_numpy(), colors[own_color]),
                 ("Opposite state", arm_table["opposite_state_score"].to_numpy(),
                  config["opposite_color"])],
                config,
                title,
                legend_loc,
                median_offsets[0],
                median_offsets[1],
                {"Own state": own_offsets},
            )
            style_axis(ax, x_label)
    else:
        fig, ax = plt.subplots(figsize=(float(config["figure_width"]),
                                        float(config["figure_height"])))
        draw_density_axis(
            ax,
            [(label, scores, colors[source]) for _, label, source, scores in loaded],
            config,
        )
        style_axis(ax, x_label)
        fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    stems = {
        "best": "crossdocking_score_density",
        "double-arm": "crossdocking_score_density_double_arm_all_poses",
        "own-state": "own_state_best_docking_score_density",
    }
    stem = stems[input_mode]
    for suffix in ("png", "svg"):
        fig.savefig(output_dir / f"{stem}.{suffix}",
                    dpi=int(config["dpi"]), bbox_inches="tight")
    plt.close(fig)

    for arm, _, _, scores in loaded:
        print(
            f"{arm}: n={scores.size}, mean={scores.mean():.3f}, "
            f"median={np.median(scores):.3f}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot density distributions for the two cross-docking arms."
    )
    parser.add_argument("--config", type=Path, required=True,
                        help="Path to the key=value .in configuration file.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    cfg = parse_config(args.config)
    mode = cfg["input_mode"].lower()
    if mode not in {"best", "double-arm", "own-state"}:
        raise ValueError("input_mode must be best, double-arm, or own-state")
    plot_density(
        resolve_path(cfg["results_dir"], args.config),
        resolve_path(cfg["output_dir"], args.config),
        mode,
        cfg,
    )
