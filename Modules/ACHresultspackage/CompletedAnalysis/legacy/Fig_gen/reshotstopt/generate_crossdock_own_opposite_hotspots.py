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


ARM_SPECS = {
    "L_PPS_to_R_PR": {"ligand": "PPS", "own": "PPS", "opposite": "PR"},
    "L_PR_to_R_PPS": {"ligand": "PR", "own": "PR", "opposite": "PPS"},
}
COLORS = {"PPS": "#D9534F", "PR": "#2878B5", "opposite": "#9A9A9A"}
REFERENCE_RESIDUES = {"ALA91", "TYR164", "ASP168", "HIS666", "ASN711", "ARG712", "GLU774"}
POCKET_SPECS = {
    "PPS": (
        ("P2", "PPS_Protein_out/pockets/pocket2_atm.pdb", "#00AFC1"),
        ("P38", "PPS_Protein_out/pockets/pocket38_atm.pdb", "#4D4D4D"),
        ("P51", "PPS_Protein_out/pockets/pocket51_atm.pdb", "#F28E2B"),
    ),
    "PR": (
        ("P11", "PR_potein_out/pockets/pocket11_atm.pdb", "#C2185B"),
        ("P15", "PR_potein_out/pockets/pocket15_atm.pdb", "#C62828"),
        ("P60", "PR_potein_out/pockets/pocket60_atm.pdb", "#4D4D4D"),
    ),
}


def parse_config(path: Path) -> dict[str, str]:
    cfg: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(f"{path}:{lineno}: expected key=value")
            key, value = line.split("=", 1)
            cfg[key.strip().lower()] = value.strip()
    required = ("own_result_dir", "opposite_result_dir", "output_dir")
    missing = [key for key in required if not cfg.get(key)]
    if missing:
        raise ValueError(f"Missing config keys: {missing}")
    cfg.setdefault("min_prevalence", "0.10")
    cfg.setdefault("top_n_per_condition", "12")
    cfg.setdefault("top_n_delta", "12")
    cfg.setdefault("pps_color", COLORS["PPS"])
    cfg.setdefault("pr_color", COLORS["PR"])
    cfg.setdefault("opposite_color", COLORS["opposite"])
    cfg.setdefault("dpi", "300")
    cfg.setdefault("residue_gap_ratio", "0.24")
    cfg.setdefault("prevalence_panel_wspace", "0.02")
    cfg.setdefault("panel_header_y", "1.0")
    cfg.setdefault("main_title_y", "1.04")
    cfg.setdefault("show_fpocket", "on")
    cfg.setdefault("fpocket_root", r"Q:\coding_dir\701_Project\Resultsbin\fpocket_res")
    return cfg


def resolve(value: str, config_path: Path) -> Path:
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    return path if path.is_absolute() else config_path.resolve().parent / path


def residue_column(frame: pd.DataFrame) -> pd.Series:
    number = pd.to_numeric(frame["original_residue_number"], errors="coerce")
    return frame["residue_name"].astype(str) + number.astype("Int64").astype(str)


def load_cohorts(opposite_metadata: Path) -> dict[str, set[str]]:
    meta = pd.read_csv(opposite_metadata, usecols=["crossdock_arm", "molecule_id", "pose_read_status"])
    meta = meta.loc[meta["pose_read_status"].eq("valid")].copy()
    return {
        arm: set(meta.loc[meta["crossdock_arm"].eq(arm), "molecule_id"].astype(str))
        for arm in ARM_SPECS
    }


def select_best_own_poses(metadata_path: Path, cohorts: dict[str, set[str]]) -> pd.DataFrame:
    usecols = ["population_state", "molecule_id", "pose_index", "docking_score", "pose_read_status"]
    meta = pd.read_csv(metadata_path, usecols=usecols)
    meta = meta.loc[meta["pose_read_status"].eq("valid")].copy()
    meta["molecule_id"] = meta["molecule_id"].astype(str)
    meta["docking_score"] = pd.to_numeric(meta["docking_score"], errors="coerce")
    selected = []
    for arm, spec in ARM_SPECS.items():
        part = meta.loc[
            meta["population_state"].eq(spec["ligand"])
            & meta["molecule_id"].isin(cohorts[arm])
        ].dropna(subset=["docking_score"])
        part = part.sort_values("docking_score").drop_duplicates("molecule_id", keep="first")
        part = part.assign(crossdock_arm=arm)
        selected.append(part)
    return pd.concat(selected, ignore_index=True)


def summarize_own(interactions_path: Path, selected: pd.DataFrame,
                  denominators: dict[str, int]) -> pd.DataFrame:
    usecols = ["population_state", "pose_index", "molecule_id", "residue_name",
               "original_residue_number", "interaction_type"]
    interactions = pd.read_csv(interactions_path, usecols=usecols)
    interactions["molecule_id"] = interactions["molecule_id"].astype(str)
    summaries = []
    for arm, spec in ARM_SPECS.items():
        keys = selected.loc[selected["crossdock_arm"].eq(arm),
                            ["population_state", "pose_index", "molecule_id"]]
        part = interactions.loc[interactions["population_state"].eq(spec["ligand"])].merge(
            keys, on=["population_state", "pose_index", "molecule_id"], how="inner",
            validate="many_to_one",
        )
        part["residue"] = residue_column(part)
        unique = part.dropna(subset=["residue"])[["molecule_id", "residue"]].drop_duplicates()
        counts = unique.groupby("residue")["molecule_id"].nunique()
        summary = counts.rename("molecule_count").reset_index()
        summary["prevalence"] = summary["molecule_count"] / denominators[arm]
        summary["crossdock_arm"] = arm
        summary["condition"] = "own"
        summaries.append(summary)
    return pd.concat(summaries, ignore_index=True)


def summarize_opposite(interactions_path: Path, denominators: dict[str, int]) -> pd.DataFrame:
    usecols = ["crossdock_arm", "molecule_id", "residue_name", "original_residue_number",
               "interaction_type"]
    interactions = pd.read_csv(interactions_path, usecols=usecols)
    interactions["molecule_id"] = interactions["molecule_id"].astype(str)
    interactions["residue"] = residue_column(interactions)
    unique = interactions.dropna(subset=["residue"])[
        ["crossdock_arm", "molecule_id", "residue"]
    ].drop_duplicates()
    summary = unique.groupby(["crossdock_arm", "residue"])["molecule_id"].nunique()
    summary = summary.rename("molecule_count").reset_index()
    summary["prevalence"] = summary.apply(
        lambda row: row["molecule_count"] / denominators[row["crossdock_arm"]], axis=1
    )
    summary["condition"] = "opposite"
    return summary


def residue_sort_key(residue: str) -> tuple[int, str]:
    digits = "".join(char for char in residue if char.isdigit())
    name = "".join(char for char in residue if not char.isdigit())
    return (int(digits) if digits else 10**9, name)


def read_pocket_residues(path: Path) -> set[str]:
    residues: set[str] = set()
    if not path.exists():
        print(f"Warning: Fpocket file not found: {path}")
        return residues
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            name = line[17:20].strip()
            number = line[22:26].strip()
            if name and number:
                residues.add(f"{name}{number}")
    return residues


def load_pockets(cfg: dict[str, str], config_path: Path) -> dict[str, list[tuple[str, set[str], str]]]:
    if cfg["show_fpocket"].lower() not in {"1", "true", "yes", "on"}:
        return {state: [] for state in POCKET_SPECS}
    root = resolve(cfg["fpocket_root"], config_path)
    return {
        state: [(label, read_pocket_residues(root / relative), color)
                for label, relative, color in specs]
        for state, specs in POCKET_SPECS.items()
    }


def prepare_table(summary: pd.DataFrame, cfg: dict[str, str]) -> tuple[pd.DataFrame, list[str]]:
    wide = summary.pivot_table(
        index="residue", columns=["crossdock_arm", "condition"], values="prevalence",
        fill_value=0.0,
    )
    expected = pd.MultiIndex.from_product([ARM_SPECS, ("own", "opposite")])
    wide = wide.reindex(columns=expected, fill_value=0.0)
    selected = set(wide.index[wide.max(axis=1) >= float(cfg["min_prevalence"])])
    top_n = int(cfg["top_n_per_condition"])
    for column in wide.columns:
        selected.update(wide[column].nlargest(top_n).index)
    top_delta = int(cfg["top_n_delta"])
    for arm in ARM_SPECS:
        delta = (wide[(arm, "opposite")] - wide[(arm, "own")]).abs()
        selected.update(delta.nlargest(top_delta).index)
    selected.update(REFERENCE_RESIDUES.intersection(wide.index))
    residues = sorted(selected, key=residue_sort_key)
    wide = wide.reindex(residues).fillna(0.0)
    for arm in ARM_SPECS:
        wide[(arm, "delta")] = wide[(arm, "opposite")] - wide[(arm, "own")]
    return wide, residues


def plot_prevalence(wide: pd.DataFrame, residues: list[str], output: Path,
                    cfg: dict[str, str], denominators: dict[str, int],
                    pockets: dict[str, list[tuple[str, set[str], str]]]) -> None:
    colors = {"PPS": cfg["pps_color"], "PR": cfg["pr_color"]}
    prevalence_columns = [
        column for column in wide.columns if column[1] in {"own", "opposite"}
    ]
    maximum = max(0.1, float(wide[prevalence_columns].max().max()))
    xmax = np.ceil(maximum * 10) / 10
    height = max(7.0, 0.34 * len(residues) + 2.0)
    fig = plt.figure(figsize=(13.5, height), layout="constrained")
    grid = fig.add_gridspec(
        1,
        3,
        width_ratios=(1, float(cfg["residue_gap_ratio"]), 1),
        wspace=float(cfg["prevalence_panel_wspace"]),
    )
    left = fig.add_subplot(grid[0, 0])
    center = fig.add_subplot(grid[0, 1], sharey=left)
    right = fig.add_subplot(grid[0, 2], sharey=left)
    y = np.arange(len(residues))
    axes = (left, right)
    for ax, (arm, spec) in zip(axes, ARM_SPECS.items()):
        own = wide[(arm, "own")].to_numpy() * 100
        opposite = wide[(arm, "opposite")].to_numpy() * 100
        own_color = colors[spec["ligand"]]
        ax.barh(y - 0.18, own, height=0.34, color=own_color, label="Own state")
        ax.barh(y + 0.18, opposite, height=0.34, color=cfg["opposite_color"],
                label="Opposite state")
        ax.set_xlim(0, xmax * 100)
        ax.grid(axis="x", color="#E3E3E3", linewidth=0.7)
        ax.set_axisbelow(True)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.set_xlabel("Molecule prevalence (%)")
        ax.legend(
            handles=[Patch(color=own_color, label="Own state"),
                     Patch(color=cfg["opposite_color"], label="Opposite state")],
            frameon=False,
            loc="upper left" if ax is left else "upper right",
        )
        ax.tick_params(axis="y", left=False, labelleft=False)

    left.invert_xaxis()
    header_y = float(cfg["panel_header_y"])
    left.set_title(
        f"PPS ligands (n={denominators['L_PPS_to_R_PR']})\nOwn PPS / opposite PR",
        y=header_y,
    )
    right.set_title(
        f"PR ligands (n={denominators['L_PR_to_R_PPS']})\nOwn PR / opposite PPS",
        y=header_y,
    )

    center.set_xlim(0, 1)
    center.set_ylim(-0.6, len(residues) - 0.4)
    center.invert_yaxis()
    center.axis("off")
    center.set_title("Pocket          Residue          Pocket", y=header_y, fontsize=10)
    for row, residue in enumerate(residues):
        center.text(0.5, row, residue, ha="center", va="center", fontsize=9)
        for state, x in (("PPS", 0.13), ("PR", 0.87)):
            hits = [(label, color) for label, members, color in pockets[state]
                    if residue in members]
            if not hits:
                center.text(x, row, "–", ha="center", va="center", fontsize=8,
                            color="#BBBBBB")
            elif len(hits) == 1:
                label, color = hits[0]
                center.text(x, row, label, ha="center", va="center", fontsize=7.5,
                            color="white", fontweight="bold",
                            bbox={"boxstyle": "square,pad=0.18", "facecolor": color,
                                  "edgecolor": "none"})
            else:
                center.text(x, row, "/".join(label.removeprefix("P") for label, _ in hits),
                            ha="center", va="center", fontsize=7.2, color="#333333",
                            fontweight="bold")

    pocket_handles = [
        Patch(facecolor=color, edgecolor="none", label=f"{state} {label}")
        for state in ("PPS", "PR")
        for label, _, color in pockets[state]
    ]
    if pocket_handles:
        fig.legend(
            handles=pocket_handles, loc="upper center", bbox_to_anchor=(0.5, 1.012),
            ncol=6, frameon=False, fontsize=8, handlelength=1.1,
            columnspacing=1.15, handletextpad=0.4,
        )

    fig.suptitle(
        "Residue interaction prevalence: own vs opposite receptor state",
        y=float(cfg["main_title_y"]),
    )
    for suffix in ("png", "svg"):
        fig.savefig(output.with_suffix(f".{suffix}"), dpi=int(cfg["dpi"]), bbox_inches="tight")
    plt.close(fig)


def plot_delta(wide: pd.DataFrame, residues: list[str], output: Path,
               cfg: dict[str, str]) -> None:
    matrix = np.column_stack([
        wide[(arm, "delta")].to_numpy() * 100 for arm in ARM_SPECS
    ])
    limit = max(1.0, float(np.ceil(np.abs(matrix).max())))
    height = max(7.0, 0.34 * len(residues) + 2.0)
    fig, ax = plt.subplots(figsize=(7.2, height))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
    ax.set_yticks(np.arange(len(residues)), residues)
    ax.set_xticks(
        [0, 1],
        ["PPS ligands\nPPS → PR", "PR ligands\nPR → PPS"],
    )
    ax.set_title("Change in residue interaction prevalence")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix[row, column]
            color = "white" if abs(value) > limit * 0.55 else "#222222"
            ax.text(column, row, f"{value:+.1f}", ha="center", va="center",
                    fontsize=8, color=color)
    colorbar = fig.colorbar(image, ax=ax, shrink=0.75, pad=0.03)
    colorbar.set_label("Opposite − own prevalence (percentage points)")
    ax.tick_params(length=0)
    fig.tight_layout()
    for suffix in ("png", "svg"):
        fig.savefig(output.with_suffix(f".{suffix}"), dpi=int(cfg["dpi"]), bbox_inches="tight")
    plt.close(fig)


def main(config_path: Path) -> None:
    cfg = parse_config(config_path)
    own_root = resolve(cfg["own_result_dir"], config_path) / "aggregate"
    opposite_root = resolve(cfg["opposite_result_dir"], config_path) / "aggregate"
    output_dir = resolve(cfg["output_dir"], config_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    pockets = load_pockets(cfg, config_path)

    cohorts = load_cohorts(opposite_root / "pose_metadata.csv")
    selected_own = select_best_own_poses(own_root / "pose_metadata.csv", cohorts)
    denominators = {
        arm: len(cohorts[arm].intersection(set(
            selected_own.loc[selected_own["crossdock_arm"].eq(arm), "molecule_id"]
        )))
        for arm in ARM_SPECS
    }
    own = summarize_own(own_root / "interactions_long.csv", selected_own, denominators)
    opposite = summarize_opposite(opposite_root / "interactions_long.csv", denominators)
    summary = pd.concat([own, opposite], ignore_index=True)
    wide, residues = prepare_table(summary, cfg)

    flat = wide.copy()
    flat.columns = [f"{arm}__{condition}" for arm, condition in flat.columns]
    for state, assignments in pockets.items():
        for label, members, _ in assignments:
            flat[f"fpocket_{state}_{label}"] = [residue in members for residue in flat.index]
    flat.to_csv(output_dir / "crossdock_own_opposite_hotspot_data.csv")
    selected_own.to_csv(output_dir / "crossdock_matched_best_own_poses.csv", index=False)
    plot_prevalence(
        wide, residues, output_dir / "crossdock_own_opposite_residue_hotspots",
        cfg, denominators, pockets,
    )
    plot_delta(
        wide, residues, output_dir / "crossdock_own_opposite_residue_delta",
        cfg,
    )
    print(f"Matched molecules: {denominators}")
    print(f"Selected residues: {len(residues)}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    main(args.config)
