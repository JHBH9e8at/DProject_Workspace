import argparse
import os
from datetime import datetime
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent / ".matplotlib"),
)

import matplotlib
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCATTER_DESC = [
    "desc_MolWt",
    "desc_NumRotatableBonds",
    "desc_CLogP",
    "desc_TPSA",
]
XLIMS = {
    "desc_MolWt": (100, 1000),
    "desc_NumRotatableBonds": (0, 17),
    "desc_CLogP": (-3, 9),
    "desc_TPSA": (0, 300),
}

REQUIRED_CONFIG_KEYS = ["prdir", "ppsdir", "outdir"]
COLORS = {"PR": "#2878B5", "PPS": "#D9534F"}
STEP_CMAPS = {
    source: LinearSegmentedColormap.from_list(
        f"{source}_steps",
        ["#F2F2F2", color],
    )
    for source, color in COLORS.items()
}


def parse_config(path: str) -> dict:
    config = {}
    with open(path, "r", encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(
                    f"{path}:{lineno}: invalid config line; expected key=value"
                )
            key, value = line.split("=", 1)
            config[key.strip().lower().replace("-", "_")] = value.strip()

    missing = [key for key in REQUIRED_CONFIG_KEYS if not config.get(key)]
    if missing:
        raise ValueError(f"Config is missing required key(s): {missing}")

    config["last_n"] = int(config.get("last_n", "50"))
    if config["last_n"] < 1:
        raise ValueError("last_n must be at least 1")

    for key in ("threshold", "pr_threshold", "pps_threshold"):
        if config.get(key, "") != "":
            config[key] = float(config[key])
        else:
            config[key] = None

    return config


def resolve_input_path(value: str, jn: str, config_path: str) -> str:
    value = os.path.expandvars(os.path.expanduser(value))
    if not os.path.isabs(value):
        value = os.path.join(os.path.dirname(os.path.abspath(config_path)), value)
    value = os.path.normpath(value)

    if os.path.isdir(value):
        value = os.path.join(value, f"{jn}_scores.csv")
    if not os.path.isfile(value):
        raise FileNotFoundError(f"{jn} input CSV not found: {value}")
    return value


def resolve_output_path(value: str, config_path: str) -> str:
    value = os.path.expandvars(os.path.expanduser(value))
    if not os.path.isabs(value):
        value = os.path.join(os.path.dirname(os.path.abspath(config_path)), value)
    return os.path.normpath(value)


def validate_input(df: pd.DataFrame, score_col: str) -> None:
    required = ["step", score_col, *SCATTER_DESC]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Input CSV is missing required column(s): {missing}")


def prepare_data(ind: str, jn: str, last_n: int) -> dict:
    df = pd.read_csv(ind)
    score_col = f"{jn}_r_i_docking_score"
    validate_input(df, score_col)

    plot_cols = ["step", score_col, *SCATTER_DESC]
    work = df[plot_cols].copy()
    for col in plot_cols:
        work[col] = pd.to_numeric(work[col], errors="coerce")

    loaded_rows = len(work)
    work = work.dropna(subset=["step", score_col])
    work = work[work[score_col] != 0].copy()
    if work.empty:
        raise ValueError(
            f"No plottable {jn} rows remain after excluding missing values and "
            f"{score_col} == 0"
        )

    max_step = work["step"].max()
    first_last_step = max(work["step"].min(), max_step - last_n + 1)
    earlier = work[work["step"] < first_last_step]
    recent = work[work["step"] >= first_last_step]

    print(f"[{jn}] Loaded rows: {loaded_rows}")
    print(f"[{jn}] Plottable rows (docking score != 0): {len(work)}")
    print(
        f"[{jn}] Last {last_n} iteration window: step {first_last_step:g} "
        f"to {max_step:g} ({len(recent)} rows)"
    )

    return {
        "jn": jn,
        "score_col": score_col,
        "earlier": earlier,
        "recent": recent,
        "first_last_step": first_last_step,
        "max_step": max_step,
    }


def draw_row(
    axes,
    prepared: dict,
    last_n: int,
    threshold: float | None,
    step_norm: Normalize,
):
    score_col = prepared["score_col"]
    earlier = prepared["earlier"]
    recent = prepared["recent"]
    first_last_step = prepared["first_last_step"]
    job_name = prepared["jn"]
    recent_scatter = None

    for ax, col in zip(axes, SCATTER_DESC):
        earlier_xy = earlier[[col, score_col]].dropna()
        recent_xy = recent[[col, score_col, "step"]].dropna()

        if not earlier_xy.empty:
            ax.scatter(
                earlier_xy[col],
                earlier_xy[score_col],
                color="lightgray",
                alpha=0.2,
                s=8,
                linewidths=0,
                rasterized=True,
                label="_nolegend_",
            )

        if not recent_xy.empty:
            recent_scatter = ax.scatter(
                recent_xy[col],
                recent_xy[score_col],
                c=recent_xy["step"],
                cmap=STEP_CMAPS[job_name],
                norm=step_norm,
                alpha=0.65,
                s=12,
                linewidths=0,
                rasterized=True,
                label="_nolegend_",
            )

        if threshold is not None:
            ax.axhline(
                threshold,
                color="black",
                linewidth=1,
                linestyle="--",
                alpha=0.8,
                label=f"{job_name} docking threshold = {threshold:g}",
            )

        # ax.set_ylabel("Docking_Score")
        ax.grid(alpha=0.15)

    return recent_scatter


def plot_pr_pps_comparison(
    pr_ind: str,
    pps_ind: str,
    outdir: str,
    last_n: int = 50,
    pr_threshold: float | None = None,
    pps_threshold: float | None = None,) -> str:
    datasets = [
        prepare_data(pr_ind, "PR", last_n),
        prepare_data(pps_ind, "PPS", last_n),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(21, 10))
    fig.subplots_adjust(
        left=0.06,
        right=0.88,
        bottom=0.08,
        top=0.86,
        wspace=0.18,
        hspace=0.0,
    )
    thresholds = [pr_threshold, pps_threshold]
    step_min = min(data["recent"]["step"].min() for data in datasets)
    step_max = max(data["recent"]["step"].max() for data in datasets)
    step_norm = Normalize(vmin=step_min, vmax=step_max)

    for row_index, (prepared, threshold) in enumerate(zip(datasets, thresholds)):
        row_axes = axes[row_index, :]
        draw_row(row_axes, prepared, last_n, threshold, step_norm)

    fig.text(
        0.012,
        0.47,
        "Docking score (kcal/mol)",
        rotation=90,
        va="center",
        ha="center",
        fontsize=13,
    )

    for col_index, descriptor in enumerate(SCATTER_DESC):
        axes[0, col_index].set_xlabel("")
        axes[0, col_index].tick_params(
            axis="x",
            bottom=True,
            labelbottom=False,
        )
        axes[1, col_index].set_xlabel(descriptor, fontsize=12)
        for row_index in range(2):
            axes[row_index, col_index].set_xlim(*XLIMS[descriptor])

    legend_items = {}
    legend_items["Y-axis = Docking score (kcal/mol)"] = Line2D(
        [],
        [],
        linestyle="none",
        marker="",
        color="none",
    )
    for row_index in range(2):
        handles, labels = axes[row_index, 0].get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            if label and not label.startswith("_"):
                legend_items.setdefault(label, handle)
    if legend_items:
        fig.legend(
            legend_items.values(),
            legend_items.keys(),
            loc="upper right",
            bbox_to_anchor=(0.99, 0.99),
            frameon=False,
        )

    # Use one shared set of tick labels between two tightly grouped colorbars.
    # The bars share the same normalization, so duplicated scales are unnecessary.
    pps_cax = fig.add_axes([0.900, 0.15, 0.012, 0.69])
    pr_cax = fig.add_axes([0.940, 0.15, 0.012, 0.69])

    pps_mappable = ScalarMappable(norm=step_norm, cmap=STEP_CMAPS["PPS"])
    pps_mappable.set_array([])
    pps_colorbar = fig.colorbar(pps_mappable, cax=pps_cax)
    pps_colorbar.ax.set_title("PPS", pad=8)
    pps_colorbar.ax.yaxis.set_ticks_position("right")
    pps_colorbar.ax.tick_params(labelright=True, labelleft=False, pad=7)

    pr_mappable = ScalarMappable(norm=step_norm, cmap=STEP_CMAPS["PR"])
    pr_mappable.set_array([])
    pr_colorbar = fig.colorbar(pr_mappable, cax=pr_cax)
    pr_colorbar.ax.set_title("PR", pad=8)
    pr_colorbar.set_ticks([])

    pr = datasets[0]
    pps = datasets[1]
    fig.suptitle(
        f"Physicochemical Property Distribution of the Last {last_n} Iterations\n"
        f"[PR steps {pr['first_last_step']:g}–{pr['max_step']:g}; "
        f"PPS steps {pps['first_last_step']:g}–{pps['max_step']:g}]",
        fontsize=20,
    )

    os.makedirs(outdir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = os.path.join(
        outdir,
        f"PR_PPS_physchem_docking_scatter_last{last_n}_{timestamp}.png",
    )
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")
    return output_path


def run_from_config(config_path: str) -> str:
    config = parse_config(config_path)
    pr_ind = resolve_input_path(config["prdir"], "PR", config_path)
    pps_ind = resolve_input_path(config["ppsdir"], "PPS", config_path)
    outdir = resolve_output_path(config["outdir"], config_path)

    common_threshold = config["threshold"]
    pr_threshold = config["pr_threshold"]
    pps_threshold = config["pps_threshold"]
    if pr_threshold is None:
        pr_threshold = common_threshold
    if pps_threshold is None:
        pps_threshold = common_threshold

    return plot_pr_pps_comparison(
        pr_ind=pr_ind,
        pps_ind=pps_ind,
        outdir=outdir,
        last_n=config["last_n"],
        pr_threshold=pr_threshold,
        pps_threshold=pps_threshold,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Generate a 2x4 PR/PPS comparison of physicochemical descriptors "
            "against docking score using settings from a .in config file."
        )
    )
    parser.add_argument("--config", required=True, help="Path to the .in config file")
    args = parser.parse_args()
    run_from_config(args.config)


# Config example:
# prdir=/path/to/PR
# ppsdir=/path/to/PPS
# outdir=/path/to/figures
# last_n=50
# threshold=-6
#
# Run:
# python physchem_scatter_last_iterations.py --config physchem_scatter.in
