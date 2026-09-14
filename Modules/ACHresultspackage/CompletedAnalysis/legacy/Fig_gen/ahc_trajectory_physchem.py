import argparse
import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent / ".matplotlib"),
)

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_OUTPUT_DIR = Path(r"Q:\coding_dir\701_Project\Writeup_Figures")
DEFAULT_REF_CSV = Path(
    r"Q:\coding_dir\701_Project\workspace\AHC_Related"
    r"\ref_structures\ref_desc.csv"
)
DEFAULT_YLIM = (-14.0, -5.0)
DIST_DESC = [
    "desc_MolWt",
    "desc_NumRotatableBonds",
    "desc_CLogP",
    "desc_TPSA",
]
COLORS = {"PR": "#2878B5", "PPS": "#D9534F"}
REF_MAP = {
    "OM": "OmecamtivMecarbil",
    "MAV": "Mavacamten",
    "AFI": "Aficamten",
}
REF_COLORS = {"OM": "#7F3C8D", "MAV": "#E6A700", "AFI": "#2CA02C"}


def validate_columns(df: pd.DataFrame, required: list[str], input_path: Path) -> None:
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"{input_path}: missing required columns: {missing}")


def load_data(input_path: Path, job_name: str) -> tuple[pd.DataFrame, str]:
    df = pd.read_csv(input_path)
    score_col = f"{job_name}_r_i_docking_score"
    validate_columns(df, ["step", score_col, *DIST_DESC], input_path)

    numeric_columns = ["step", score_col, *DIST_DESC]
    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    if df["step"].isna().any() or df[score_col].isna().any():
        raise ValueError(f"{input_path}: step or docking score contains non-numeric values.")
    return df, score_col


def load_refs(ref_csv: Path | None, ref_keys: list[str]) -> dict[str, pd.Series]:
    if not ref_keys:
        return {}
    if ref_csv is None:
        raise ValueError("--ref-csv is required when --ref is used.")
    if not ref_csv.is_file():
        raise FileNotFoundError(f"Reference CSV not found: {ref_csv}")

    ref_df = pd.read_csv(ref_csv)
    validate_columns(ref_df, ["name", *DIST_DESC], ref_csv)
    refs = {}
    for key in ref_keys:
        name = REF_MAP[key]
        rows = ref_df.loc[ref_df["name"] == name]
        if rows.empty:
            raise ValueError(f"{ref_csv}: reference '{name}' ({key}) was not found.")
        refs[key] = rows.iloc[0]
    return refs


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_dir / stem}.[png|svg]")


def plot_trajectory(
    df: pd.DataFrame,
    score_col: str,
    job_name: str,
    output_dir: Path,
    ylim: tuple[float, float] | None,
) -> None:
    color = COLORS[job_name]
    failed = df.loc[df[score_col] == 0].groupby("step").size()
    valid = df.loc[df[score_col] != 0].copy()
    if valid.empty:
        raise ValueError(f"{job_name}: no non-zero docking scores are available.")

    grouped = valid.groupby("step")[score_col]
    mean = grouped.mean()
    std = grouped.std().fillna(0)
    all_steps = pd.Index(sorted(df["step"].unique()))
    failed = failed.reindex(all_steps, fill_value=0)

    fig, (ax_score, ax_fail) = plt.subplots(
        2,
        1,
        figsize=(11, 7.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
    )
    ax_score.plot(mean.index, mean, color=color, linewidth=1, label="Mean") #ttt
    ax_score.fill_between(
        mean.index,
        mean - std,
        mean + std,
        color=color,
        alpha=0.22,
        linewidth=0,
        label="±1 SD",
    )
    if ylim is not None:
        ax_score.set_ylim(*ylim)
    ax_score.set_ylabel("Mean docking score (kcal/mol)")
    ax_score.legend(frameon=False, ncol=2)

    ax_fail.bar(
        failed.index,
        failed.values,
        width=0.85,
        color=color,
        alpha=0.75,
        label="Docking failures",
    )
    ax_fail.set_ylabel("Failure count")
    ax_fail.set_xlabel("Generation step")
    ax_fail.legend(frameon=False)

    for ax in (ax_score, ax_fail):
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)
        ax.set_axisbelow(True)

    fig.suptitle(f"{job_name} docking score trajectory")
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.10, top=0.92)
    save_figure(fig, output_dir, f"{job_name}_docking_score_trajectory")
    print(f"{job_name} docking failures: {int(failed.sum())}")


def plot_distribution(
    df: pd.DataFrame,
    job_name: str,
    output_dir: Path,
    refs: dict[str, pd.Series],
    suffix: str,
    title_suffix: str,
) -> None:
    color = COLORS[job_name]
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    for ax, column in zip(axes.flat, DIST_DESC):
        values = df[column].dropna()
        ax.hist(values, bins=50, color=color, edgecolor="none", alpha=0.78)
        for key, row in refs.items():
            value = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
            if pd.notna(value):
                ax.axvline(
                    value,
                    color=REF_COLORS[key],
                    linewidth=1.5,
                    linestyle="--",
                    label=REF_MAP[key],
                )
        ax.set_title(column.removeprefix("desc_"))
        ax.set_ylabel("Count")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)
        ax.set_axisbelow(True)

    if refs:
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=len(refs), frameon=False)

    fig.suptitle(f"{job_name} physicochemical property distributions{title_suffix}")
    bottom = 0.12 if refs else 0.08
    fig.subplots_adjust(left=0.09, right=0.98, bottom=bottom, top=0.91, hspace=0.32)
    save_figure(fig, output_dir, f"{job_name}_physchem_distribution{suffix}")


def run_single(args: argparse.Namespace, input_path: Path, job_name: str) -> None:
    if not 0 < args.top_pct <= 100:
        raise ValueError("--top-pct must be greater than 0 and at most 100.")
    if args.ylim is not None and args.ylim[0] >= args.ylim[1]:
        raise ValueError("--ylim requires YMIN < YMAX.")

    df, score_col = load_data(input_path, job_name)
    refs = load_refs(args.ref_csv, args.ref)

    if args.which in ("trajectory", "both"):
        plot_trajectory(
            df,
            score_col,
            job_name,
            args.output_dir,
            tuple(args.ylim) if args.ylim else None,
        )

    if args.which in ("distribution", "both"):
        valid = df.loc[df[score_col] != 0].copy()
        if valid.empty:
            raise ValueError(f"{job_name}: no non-zero docking scores are available.")

        plot_distribution(valid, job_name, args.output_dir, refs, "", "")
        threshold = valid[score_col].quantile(args.top_pct / 100)
        top = valid.loc[valid[score_col] <= threshold].copy()
        pct_label = f"{args.top_pct:g}"
        plot_distribution(
            top,
            job_name,
            args.output_dir,
            refs,
            f"_top{pct_label}pct",
            f" — top {pct_label}% by docking score",
        )
        print(
            f"{job_name} top {pct_label}%: {len(top)} / {len(valid)} "
            f"(docking score <= {threshold:.4f})"
        )


def run(args: argparse.Namespace) -> None:
    paired_mode = args.pr_input is not None or args.pps_input is not None
    single_mode = args.input is not None or args.job_name is not None

    if paired_mode and single_mode:
        raise ValueError(
            "Use either --pr-input/--pps-input together or "
            "--input/--job-name, not both modes."
        )
    if paired_mode:
        if args.pr_input is None or args.pps_input is None:
            raise ValueError("Both --pr-input and --pps-input are required.")
        run_single(args, args.pr_input, "PR")
        run_single(args, args.pps_input, "PPS")
        return
    if args.input is None or args.job_name is None:
        raise ValueError(
            "Provide both --pr-input and --pps-input, or provide "
            "--input with --job-name."
        )
    run_single(args, args.input, args.job_name)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate an AHC docking-score trajectory and four-property "
            "distribution figures."
        )
    )
    parser.add_argument(
        "--pr-input",
        type=Path,
        default=None,
        help="PR cleaned scores CSV for combined PR/PPS execution.",
    )
    parser.add_argument(
        "--pps-input",
        type=Path,
        default=None,
        help="PPS cleaned scores CSV for combined PR/PPS execution.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Single-dataset cleaned scores CSV (legacy mode).",
    )
    parser.add_argument(
        "--output-dir",
        "--odir",
        dest="output_dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--job-name",
        choices=["PR", "PPS"],
        default=None,
        help="Dataset identity for --input single-dataset mode.",
    )
    parser.add_argument(
        "--which",
        choices=["trajectory", "distribution", "both"],
        default="both",
    )
    parser.add_argument(
        "--ylim",
        type=float,
        nargs=2,
        metavar=("YMIN", "YMAX"),
        default=list(DEFAULT_YLIM),
        help="Fixed trajectory y-axis (default: -14 -5).",
    )
    parser.add_argument(
        "--top-pct",
        type=float,
        default=1.0,
        help="Best low-tail docking-score percentage (default: 1).",
    )
    parser.add_argument(
        "--ref",
        nargs="*",
        choices=sorted(REF_MAP),
        default=["OM"],
        help="Reference molecules (default: OM). Pass --ref with no values to disable.",
    )
    parser.add_argument(
        "--ref-csv",
        type=Path,
        default=DEFAULT_REF_CSV,
        help=f"Reference descriptor CSV (default: {DEFAULT_REF_CSV}).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())


#python .\ahc_trajectory_physchem.py --pr-input "Q:\coding_dir\701_Project\Resultsbin\AHC_results\PR\PR_scores.csv" --pps-input "Q:\coding_dir\701_Project\Resultsbin\AHC_results\PPS\PPS_scores.csv" --odir "Q:\coding_dir\701_Project\Writeup_Figures\PC" --which both --ref OM