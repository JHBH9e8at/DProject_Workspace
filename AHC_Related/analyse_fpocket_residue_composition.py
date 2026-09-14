#!/usr/bin/env python3
"""Compare PR/PPS fpocket residue composition around the reference OM ligand.

The selected fpocket pockets are merged per conformational state. The complete
union is compared first, followed by the subset whose protein heavy atoms are
within 5 Angstrom of any heavy atom of the reference ligand (2OW).
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle


ROOT = Path(r"Q:\coding_dir\701_Project")
DEFAULT_FPOCKET = ROOT / "Resultsbin" / "fpocket_res"
DEFAULT_STRUCTURES = ROOT / "workspace" / "AHC_Related" / "ref_structures" / "Structures"
DEFAULT_OUTPUT = ROOT / "Writeup_Figures" / "fpocket_residue_comparison"
DEFAULT_HOTSPOTS = ROOT / "Writeup_Figures" / "resanaly" / "pr_pps_residue_hotspot_full_vs_last50_data.csv"
HOTSPOT_THRESHOLD = 0.01

STATE_CONFIG = {
    "PPS": {
        "pocket_dir": "PPS_Protein_out",
        "pockets": (2, 38, 51),
        "reference": "MYH7_PPS_OM.pdb",
    },
    "PR": {
        "pocket_dir": "PR_potein_out",  # spelling matches the existing result folder
        "pockets": (11, 15, 60),
        "reference": "MYH7_PR_OM.pdb",
    },
}

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def parse_atom_line(line: str) -> dict | None:
    if not line.startswith(("ATOM  ", "HETATM")):
        return None
    try:
        element = line[76:78].strip()
        atom_name = line[12:16].strip()
        if not element:
            element = "".join(c for c in atom_name if c.isalpha())[:1]
        return {
            "record": line[:6].strip(),
            "atom": atom_name,
            "resname": line[17:20].strip(),
            "chain": line[21].strip() or "_",
            "resnum": int(line[22:26]),
            "icode": line[26].strip(),
            "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            "element": element.upper(),
        }
    except (ValueError, IndexError):
        return None


def residue_key(atom: dict) -> tuple[str, int, str, str]:
    return atom["chain"], atom["resnum"], atom["icode"], atom["resname"]


def read_pocket_residues(pockets_dir: Path, pocket_ids: tuple[int, ...]):
    membership: dict[tuple[str, int, str, str], set[int]] = defaultdict(set)
    for pocket_id in pocket_ids:
        path = pockets_dir / f"pocket{pocket_id}_atm.pdb"
        if not path.exists():
            raise FileNotFoundError(f"Missing fpocket atom file: {path}")
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            atom = parse_atom_line(line)
            if atom and atom["record"] == "ATOM" and atom["resname"] in STANDARD_AA:
                membership[residue_key(atom)].add(pocket_id)
    return membership


def read_reference_distances(path: Path, ligand_resname: str):
    protein_atoms: dict[tuple[str, int, str, str], list[np.ndarray]] = defaultdict(list)
    ligand_atoms: list[np.ndarray] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        atom = parse_atom_line(line)
        if not atom or atom["element"] == "H":
            continue
        if atom["record"] == "HETATM" and atom["resname"] == ligand_resname:
            ligand_atoms.append(atom["xyz"])
        elif atom["record"] == "ATOM" and atom["resname"] in STANDARD_AA:
            protein_atoms[residue_key(atom)].append(atom["xyz"])
    if not ligand_atoms:
        raise ValueError(f"No heavy atoms for ligand {ligand_resname!r} in {path}")

    ligand = np.vstack(ligand_atoms)
    distances = {}
    for key, coords in protein_atoms.items():
        protein = np.vstack(coords)
        distances[key] = float(np.sqrt(((protein[:, None, :] - ligand[None, :, :]) ** 2).sum(axis=2)).min())
    return distances


def display_label(key: tuple[str, int, str, str]) -> str:
    _chain, resnum, icode, resname = key
    return f"{resname}{resnum}{icode}"


def build_rows(fpocket_root: Path, structure_root: Path, ligand_resname: str):
    state_data = {}
    all_keys = set()
    for state, cfg in STATE_CONFIG.items():
        membership = read_pocket_residues(
            fpocket_root / cfg["pocket_dir"] / "pockets", cfg["pockets"]
        )
        distances = read_reference_distances(structure_root / cfg["reference"], ligand_resname)
        state_data[state] = {"membership": membership, "distances": distances}
        all_keys.update(membership)

    rows = []
    for key in sorted(all_keys, key=lambda k: (k[1], k[2], k[0], k[3])):
        row = {
            "residue": display_label(key),
            "chain": key[0],
            "resnum": key[1],
            "insertion_code": key[2],
            "resname": key[3],
        }
        for state in ("PPS", "PR"):
            pockets = state_data[state]["membership"].get(key, set())
            distance = state_data[state]["distances"].get(key, math.nan)
            row[f"{state}_pockets"] = ";".join(map(str, sorted(pockets)))
            row[f"{state}_in_selected_pockets"] = bool(pockets)
            row[f"{state}_min_OM_distance_A"] = distance
            row[f"{state}_within_5A"] = bool(pockets) and distance <= 5.0
        rows.append(row)
    return rows


def write_csv(rows: list[dict], path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sets_for_scope(rows: list[dict], scope: str):
    suffix = {
        "All selected fpocket residues": "in_selected_pockets",
        "Within 5 Å of OM": "within_5A",
    }[scope]
    return {
        state: {row["residue"] for row in rows if row[f"{state}_{suffix}"]}
        for state in ("PPS", "PR")
    }


def read_hotspot_sets(path: Path, threshold: float = HOTSPOT_THRESHOLD):
    hotspot_sets = {"PPS": set(), "PR": set()}
    prevalence = {"PPS": {}, "PR": {}}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            for state in ("PPS", "PR"):
                full = float(row[f"full_{state}"] or 0)
                last50 = float(row[f"last50_{state}"] or 0)
                prevalence[state][row["residue"]] = last50
                if max(full, last50) >= threshold:
                    hotspot_sets[state].add(row["residue"])
    return hotspot_sets, prevalence


def plot(rows: list[dict], output: Path, hotspot_path: Path):
    scopes = ["All selected fpocket residues", "Within 5 Å of OM"]
    counts = []
    scope_sets = {}
    for scope in scopes:
        sets = sets_for_scope(rows, scope)
        scope_sets[scope] = sets
        counts.append([
            len(sets["PPS"] - sets["PR"]),
            len(sets["PPS"] & sets["PR"]),
            len(sets["PR"] - sets["PPS"]),
        ])
    hotspot_sets, prevalence = read_hotspot_sets(hotspot_path)
    shared_hotspots = hotspot_sets["PPS"] & hotspot_sets["PR"]

    fig, ax_heat = plt.subplots(figsize=(15, 6), constrained_layout=True)
    all_labels = sorted(
        set().union(*(scope_sets[s][state] for s in scopes for state in ("PPS", "PR"))),
        key=lambda label: int("".join(c for c in label if c.isdigit())),
    )
    columns = []
    matrix = []
    column_order = [
        ("All selected fpocket residues", "PPS", "Full pocket"),
        ("Within 5 Å of OM", "PPS", "5 Å to OM"),
        ("Within 5 Å of OM", "PR", "5 Å to OM"),
        ("All selected fpocket residues", "PR", "Full pocket"),
    ]
    for scope, state, column_label in column_order:
        columns.append(column_label)
        values = []
        for label in all_labels:
            if label not in scope_sets[scope][state]:
                values.append(0)  # outside this pocket/scope
            elif label not in hotspot_sets[state]:
                values.append(1)  # pocket residue without a hotspot interaction
            elif label in shared_hotspots:
                values.append(3)  # hotspot in both directed populations
            elif state == "PPS":
                values.append(2)  # PPS-only hotspot
            else:
                values.append(4)  # PR-only hotspot
        matrix.append(values)
    matrix = np.array(matrix).T
    # State colours match generate_residue_hotspot_figure.py. White denotes
    # outside-pocket residues and grey denotes pocket residues without hotspots.
    cmap = ListedColormap(["#FFFFFF", "#BDBDBD", "#D9534F", "#F1C232", "#2878B5"])
    # Transpose the matrix to provide a clockwise-rotated, landscape layout
    # while retaining upright, readable text.
    ax_heat.imshow(matrix.T, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=4)
    ax_heat.set_xticks(range(len(all_labels)), all_labels, fontsize=12, rotation=90)
    ax_heat.set_yticks(range(len(columns)), columns, fontsize=12)
    ax_heat.set_xlabel("Residue identity", labelpad=12)
    ax_heat.set_xticks(np.arange(-0.5, len(all_labels), 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=0.6)
    ax_heat.tick_params(which="minor", bottom=False, left=False)
    cell_colours = ["#FFFFFF", "#BDBDBD", "#D9534F", "#F1C232", "#2878B5"]
    merged_cells = set()
    # A residue in the 5-A subset is necessarily in the corresponding full
    # pocket. Merge those two state-specific cells into one vertical block.
    for state, row_pair, y_start in (("PPS", (0, 1), -0.5), ("PR", (2, 3), 1.5)):
        five_a_row = 1 if state == "PPS" else 2
        for col_index, label in enumerate(all_labels):
            if label not in scope_sets["Within 5 Å of OM"][state]:
                continue
            code = int(matrix[col_index, five_a_row])
            value = prevalence[state].get(label, 0.0) * 100
            ax_heat.add_patch(
                Rectangle(
                    (col_index - 0.5, y_start), 1.0, 2.0,
                    facecolor=cell_colours[code], edgecolor="white", linewidth=0.6,
                    zorder=3,
                )
            )
            merged_cells.update((col_index, row) for row in row_pair)
            if code == 3:
                ax_heat.text(
                    col_index, y_start + 1.0, f"{value:.1f}",
                    ha="center", va="center", rotation=90,
                    fontsize=6.5, color="#222222", fontweight="bold", zorder=4,
                )

    # Annotate shared hotspots that remain as single, unmerged cells.
    for row_index, (_scope, state, _column_label) in enumerate(column_order):
        for col_index, label in enumerate(all_labels):
            if (col_index, row_index) not in merged_cells and matrix[col_index, row_index] == 3:
                value = prevalence[state].get(label, 0.0) * 100
                ax_heat.text(
                    col_index, row_index, f"{value:.1f}",
                    ha="center", va="center", rotation=90,
                    fontsize=6.5, color="#222222", fontweight="bold", zorder=4,
                )
    n_residues = len(all_labels)
    for y_start, state_label in ((-0.5, "PPS"), (1.5, "PR")):
        ax_heat.add_patch(
            Rectangle(
                (-0.5, y_start), n_residues, 2.0,
                fill=False, edgecolor="#333333", linewidth=2.0,
                clip_on=False, zorder=5,
            )
        )
        ax_heat.text(
            n_residues + 1.0, y_start + 1.0, state_label,
            ha="left", va="center", fontsize=13, fontweight="bold", clip_on=False,
        )
    ax_heat.legend(
        handles=[
            Patch(facecolor="#D9534F", label="PPS hotspot"),
            Patch(facecolor="#F1C232", label="Shared hotspot"),
            Patch(facecolor="#2878B5", label="PR hotspot"),
            Patch(facecolor="#BDBDBD", label="Pocket, no hotspot"),
            Patch(facecolor="#FFFFFF", edgecolor="#BBBBBB", label="Outside pocket"),
        ],
        frameon=False,
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.015),
    )

    fig.suptitle(
        "OM binding site composition of PR and PPS state cardiac myosin",
        fontsize=15,
        fontweight="bold",
    )
    fig.savefig(output, dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    return {scope: counts[i] for i, scope in enumerate(scopes)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fpocket-root", type=Path, default=DEFAULT_FPOCKET)
    parser.add_argument("--structure-root", type=Path, default=DEFAULT_STRUCTURES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--hotspot-data", type=Path, default=DEFAULT_HOTSPOTS)
    parser.add_argument("--ligand-resname", default="2OW")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = build_rows(args.fpocket_root, args.structure_root, args.ligand_resname)
    csv_path = args.output_dir / "pr_pps_fpocket_residue_comparison.csv"
    figure_path = args.output_dir / "pr_pps_fpocket_residue_comparison.png"
    write_csv(rows, csv_path)
    summary = plot(rows, figure_path, args.hotspot_data)

    print(f"Wrote: {csv_path}")
    print(f"Wrote: {figure_path}")
    print(f"Wrote: {figure_path.with_suffix('.pdf')}")
    for scope, (pps_only, shared, pr_only) in summary.items():
        ascii_scope = scope.replace("Å", "A")
        print(f"{ascii_scope}: PPS-only={pps_only}, shared={shared}, PR-only={pr_only}")


if __name__ == "__main__":
    main()
