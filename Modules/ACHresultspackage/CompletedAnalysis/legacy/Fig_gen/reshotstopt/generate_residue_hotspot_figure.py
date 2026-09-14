from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


INTERACTIONS = Path(
    r"Q:\coding_dir\701_Project\Resultsbin\MP_Interrun\aggregate\interactions_long.csv"
)
FPOCKET_ROOT = Path(r"Q:\coding_dir\701_Project\Resultsbin\fpocket_res")
REFERENCE_ROOT = Path(r"Q:\coding_dir\701_Project\Resultsbin\reference_interaction")
OUTPUT_DIR = Path(r"Q:\coding_dir\701_Project\Writeup_Figures\resanaly")

TOTAL_POSES = {"PR": 19688, "PPS": 20420}
TOP_N_PER_STATE = 12
MIN_PREVALENCE = 0.10

POCKET_FILES = {
    "PPS P2": FPOCKET_ROOT / "PPS_Protein_out" / "pockets" / "pocket2_atm.pdb",
    "PPS P38": FPOCKET_ROOT / "PPS_Protein_out" / "pockets" / "pocket38_atm.pdb",
    "PPS P51": FPOCKET_ROOT / "PPS_Protein_out" / "pockets" / "pocket51_atm.pdb",
    "PR P11": FPOCKET_ROOT / "PR_potein_out" / "pockets" / "pocket11_atm.pdb",
    "PR P15": FPOCKET_ROOT / "PR_potein_out" / "pockets" / "pocket15_atm.pdb",
    "PR P60": FPOCKET_ROOT / "PR_potein_out" / "pockets" / "pocket60_atm.pdb",
}

POCKET_COLORS = {
    "PPS P2": "#00AFC1",   # cyan
    "PPS P38": "#4D4D4D",  # dark grey
    "PPS P51": "#F28E2B",  # orange
    "PR P11": "#C2185B",   # magenta
    "PR P15": "#C62828",   # red
    "PR P60": "#4D4D4D",   # dark grey
}

# Preserve experimentally/reference-relevant anchors even if they fall below the
# population-level display threshold in one state.
REFERENCE_RESIDUES = {
    "ALA91",
    "TYR164",
    "ASP168",
    "HIS666",
    "ASN711",
    "ARG712",
    "GLU774",
}


def read_pocket_residues(path: Path) -> set[str]:
    residues: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            name = line[17:20].strip()
            number = line[22:26].strip()
            if name and number:
                residues.add(f"{name}{number}")
    return residues


def summarize_interactions() -> tuple[pd.DataFrame, pd.DataFrame]:
    usecols = [
        "population_state",
        "molecule_id",
        "residue_name",
        "original_residue_number",
        "interaction_type",
    ]
    residue_frames: list[pd.DataFrame] = []
    typed_frames: list[pd.DataFrame] = []

    for chunk in pd.read_csv(INTERACTIONS, usecols=usecols, chunksize=250_000):
        chunk = chunk.dropna(
            subset=["population_state", "molecule_id", "residue_name", "original_residue_number"]
        ).copy()
        chunk["residue_number"] = chunk["original_residue_number"].astype(int)
        chunk["residue"] = chunk["residue_name"].astype(str) + chunk["residue_number"].astype(str)
        residue_frames.append(
            chunk[["population_state", "molecule_id", "residue"]].drop_duplicates()
        )
        typed_frames.append(
            chunk[
                ["population_state", "molecule_id", "residue", "interaction_type"]
            ].drop_duplicates()
        )

    residue_unique = pd.concat(residue_frames, ignore_index=True).drop_duplicates()
    typed_unique = pd.concat(typed_frames, ignore_index=True).drop_duplicates()

    residue_counts = (
        residue_unique.groupby(["residue", "population_state"])["molecule_id"]
        .nunique()
        .rename("poses")
        .reset_index()
    )
    residue_counts["prevalence"] = residue_counts.apply(
        lambda row: row["poses"] / TOTAL_POSES[row["population_state"]], axis=1
    )

    typed_counts = (
        typed_unique.groupby(["residue", "population_state", "interaction_type"])[
            "molecule_id"
        ]
        .nunique()
        .rename("poses")
        .reset_index()
    )
    typed_counts["prevalence"] = typed_counts.apply(
        lambda row: row["poses"] / TOTAL_POSES[row["population_state"]], axis=1
    )
    return residue_counts, typed_counts


def dominant_interaction_labels(typed_counts: pd.DataFrame) -> pd.DataFrame:
    dominant = typed_counts.loc[
        typed_counts.groupby(["residue", "population_state"])["prevalence"].idxmax()
    ].copy()
    abbreviations = {
        "HBDonor": "HBD",
        "HBAcceptor": "HBA",
        "PiStacking": "pi",
        "PiCation": "pi-cat",
        "CationPi": "cat-pi",
        "Cationic": "cat",
        "Anionic": "ani",
        "XBDonor": "XB",
    }
    dominant["dominant_type"] = dominant["interaction_type"].map(abbreviations).fillna(
        dominant["interaction_type"]
    )
    return dominant[["residue", "population_state", "dominant_type"]]


def read_reference_interactions() -> pd.DataFrame:
    abbreviations = {
        "HBDonor": "HBD",
        "HBAcceptor": "HBA",
        "PiStacking": "pi",
        "PiCation": "pi-cat",
        "CationPi": "cat-pi",
        "Cationic": "cat",
        "Anionic": "ani",
        "XBDonor": "XB",
    }
    records: list[dict[str, str]] = []
    for state in ("PPS", "PR"):
        path = REFERENCE_ROOT / state / "interactions_long.csv"
        frame = pd.read_csv(
            path,
            usecols=["residue_name", "original_residue_number", "interaction_type"],
        ).dropna()
        frame["residue"] = (
            frame["residue_name"].astype(str)
            + frame["original_residue_number"].astype(int).astype(str)
        )
        frame["label"] = frame["interaction_type"].map(abbreviations).fillna(
            frame["interaction_type"]
        )
        for residue, group in frame.groupby("residue", sort=False):
            labels = list(dict.fromkeys(group["label"].astype(str)))
            records.append({"residue": residue, "population_state": state, "reference": "/".join(labels)})
    return pd.DataFrame(records).pivot(
        index="residue", columns="population_state", values="reference"
    )


def select_residues(summary: pd.DataFrame) -> list[str]:
    wide = summary.pivot(index="residue", columns="population_state", values="prevalence").fillna(0)
    for state in ("PR", "PPS"):
        if state not in wide:
            wide[state] = 0.0

    selected = set(wide.index[wide.max(axis=1) >= MIN_PREVALENCE])
    for state in ("PR", "PPS"):
        selected.update(wide[state].nlargest(TOP_N_PER_STATE).index)
    selected.update(REFERENCE_RESIDUES.intersection(wide.index))

    def residue_sort_key(residue: str) -> tuple[int, str]:
        digits = "".join(character for character in residue if character.isdigit())
        name = "".join(character for character in residue if not character.isdigit())
        return (int(digits), name)

    return sorted(selected, key=residue_sort_key)


def build_figure(
    residue_counts: pd.DataFrame,
    typed_counts: pd.DataFrame,
    pocket_membership: dict[str, set[str]],
    reference_interactions: pd.DataFrame,
) -> pd.DataFrame:
    residues = select_residues(residue_counts)
    prevalence = (
        residue_counts.pivot(index="residue", columns="population_state", values="prevalence")
        .fillna(0)
        .reindex(residues)
    )
    for state in ("PR", "PPS"):
        if state not in prevalence:
            prevalence[state] = 0.0

    dominant = dominant_interaction_labels(typed_counts).pivot(
        index="residue", columns="population_state", values="dominant_type"
    )
    dominant = dominant.reindex(residues).fillna("-")
    reference = reference_interactions.reindex(residues).fillna("-")
    for state in ("PPS", "PR"):
        if state not in reference:
            reference[state] = "-"

    pocket_names = list(pocket_membership)
    membership = pd.DataFrame(
        {
            pocket: [residue in pocket_membership[pocket] for residue in residues]
            for pocket in pocket_names
        },
        index=residues,
    )

    n = len(residues)
    scale = 2
    width = 1800 * scale
    left_outer = 130 * scale
    left_inner = 650 * scale
    backbone_left = 650 * scale
    backbone_right = 1150 * scale
    right_inner = 1150 * scale
    right_outer = 1670 * scale
    top = 150 * scale
    row_h = 32 * scale
    bottom = 95 * scale
    height = top + n * row_h + bottom

    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font_dir = Path(r"C:\Windows\Fonts")
    regular = ImageFont.truetype(str(font_dir / "arial.ttf"), 13 * scale)
    small = ImageFont.truetype(str(font_dir / "arial.ttf"), 10 * scale)
    bold = ImageFont.truetype(str(font_dir / "arialbd.ttf"), 13 * scale)
    title_font = ImageFont.truetype(str(font_dir / "arialbd.ttf"), 19 * scale)

    colors = {
        "PR": "#2878B5",
        "PPS": "#D9534F",
        "connector": "#B7B7B7",
        "grid": "#E7E7E7",
        "text": "#222222",
        "muted": "#666666",
        "pocket": "#6B5B95",
        "empty": "#F1F1F1",
    }

    draw.text((width / 2, 30 * scale), "PR/PPS residue interaction hotspots", fill=colors["text"], font=title_font, anchor="ma")
    panel_heading_y = 82 * scale
    draw.text(((left_outer + left_inner) / 2, panel_heading_y), "PPS prevalence", fill=colors["PPS"], font=bold, anchor="ma")
    draw.text(((backbone_left + backbone_right) / 2, panel_heading_y), "Residue", fill=colors["text"], font=bold, anchor="ma")
    draw.text(((right_inner + right_outer) / 2, panel_heading_y), "PR prevalence", fill=colors["PR"], font=bold, anchor="ma")
    pps_interaction_x = backbone_left + 35 * scale
    pps_pocket_x = backbone_left + 105 * scale
    pps_reference_x = backbone_left + 185 * scale
    residue_x = (backbone_left + backbone_right) / 2
    pr_reference_x = backbone_right - 185 * scale
    pr_pocket_x = backbone_right - 105 * scale
    pr_interaction_x = backbone_right - 35 * scale
    draw.text((pps_interaction_x, 118 * scale), "Interaction", fill=colors["muted"], font=small, anchor="mm")
    draw.text((pps_pocket_x, 118 * scale), "Pocket", fill=colors["muted"], font=small, anchor="mm")
    draw.text((pps_reference_x, 118 * scale), "Reference", fill=colors["muted"], font=small, anchor="mm")
    draw.text((pr_reference_x, 118 * scale), "Reference", fill=colors["muted"], font=small, anchor="mm")
    draw.text((pr_pocket_x, 118 * scale), "Pocket", fill=colors["muted"], font=small, anchor="mm")
    draw.text((pr_interaction_x, 118 * scale), "Interaction", fill=colors["muted"], font=small, anchor="mm")

    legend_y = top + n * row_h + 22 * scale
    legend_items = [
        ("PPS P2", pps_pocket_x - 75 * scale),
        ("PPS P38", pps_pocket_x - 25 * scale),
        ("PPS P51", pps_pocket_x + 30 * scale),
        ("PR P11", pr_pocket_x - 75 * scale),
        ("PR P15", pr_pocket_x - 25 * scale),
        ("PR P60", pr_pocket_x + 30 * scale),
    ]
    for pocket, x in legend_items:
        draw.rectangle((x, legend_y, x + 11 * scale, legend_y + 11 * scale), fill=POCKET_COLORS[pocket])
        draw.text(
            (x + 15 * scale, legend_y + 5 * scale),
            pocket.split()[-1],
            fill=colors["muted"],
            font=small,
            anchor="lm",
        )

    xmax = max(80, int(np.ceil(prevalence.to_numpy().max() * 100 / 10) * 10))
    for tick in range(0, xmax + 1, 10):
        x_left = left_inner - int((left_inner - left_outer) * tick / xmax)
        x_right = right_inner + int((right_outer - right_inner) * tick / xmax)
        draw.line((x_left, top - 8 * scale, x_left, top + n * row_h), fill=colors["grid"], width=1 * scale)
        draw.line((x_right, top - 8 * scale, x_right, top + n * row_h), fill=colors["grid"], width=1 * scale)
        label = str(tick)
        box = draw.textbbox((0, 0), label, font=small)
        draw.text((x_left - (box[2] - box[0]) / 2, top - 30 * scale), label, fill=colors["muted"], font=small)
        draw.text((x_right - (box[2] - box[0]) / 2, top - 30 * scale), label, fill=colors["muted"], font=small)
    draw.line((left_inner, top - 8 * scale, left_inner, top + n * row_h), fill="#999999", width=2 * scale)
    draw.line((right_inner, top - 8 * scale, right_inner, top + n * row_h), fill="#999999", width=2 * scale)

    # Frame the state-specific annotation columns so that interaction, pocket,
    # and reference information reads as one block on either side of Residue.
    annotation_top = 100 * scale
    annotation_bottom = top + n * row_h
    annotation_outline = "#B8B8B8"
    # Keep the frame clear of the outermost graph tick labels and the
    # reference columns while leaving the central residue labels unboxed.
    residue_clearance = 30 * scale
    outer_padding = 45 * scale
    draw.rectangle(
        (left_outer - outer_padding, annotation_top, residue_x - residue_clearance, annotation_bottom),
        outline=annotation_outline,
        width=2 * scale,
    )
    draw.rectangle(
        (residue_x + residue_clearance, annotation_top, right_outer + outer_padding, annotation_bottom),
        outline=annotation_outline,
        width=2 * scale,
    )
    draw.rectangle(
        (residue_x - residue_clearance, annotation_top,
         residue_x + residue_clearance, annotation_bottom),
        outline=annotation_outline,
        width=2 * scale,
    )

    for row, residue in enumerate(residues):
        yc = top + row * row_h + row_h / 2
        pps_reference = str(reference.loc[residue, "PPS"])
        pr_reference = str(reference.loc[residue, "PR"])
        pps_interaction = str(dominant.loc[residue, "PPS"])
        pr_interaction = str(dominant.loc[residue, "PR"])
        pps_reference = "" if pps_reference == "-" else pps_reference
        pr_reference = "" if pr_reference == "-" else pr_reference
        pps_interaction = "" if pps_interaction == "-" else pps_interaction
        pr_interaction = "" if pr_interaction == "-" else pr_interaction
        draw.text((residue_x, yc), residue, fill=colors["text"], font=regular, anchor="mm")
        draw.text((pps_interaction_x, yc), pps_interaction, fill=colors["muted"], font=small, anchor="mm")
        draw.text((pps_reference_x, yc), pps_reference, fill=colors["text"], font=small, anchor="mm")
        draw.text((pr_reference_x, yc), pr_reference, fill=colors["text"], font=small, anchor="mm")
        draw.text((pr_interaction_x, yc), pr_interaction, fill=colors["muted"], font=small, anchor="mm")
        for state, pocket_x in (("PPS", pps_pocket_x), ("PR", pr_pocket_x)):
            state_pockets = [name for name in pocket_names if name.startswith(state + " ")]
            hits = [name for name in state_pockets if membership.loc[residue, name]]
            if not hits:
                continue
            elif len(hits) == 1:
                pocket = hits[0]
                label = pocket.split()[-1]
                draw.rectangle(
                    (pocket_x - 24 * scale, yc - 10 * scale, pocket_x + 24 * scale, yc + 10 * scale),
                    fill=POCKET_COLORS[pocket],
                )
                draw.text((pocket_x, yc), label, fill="white", font=small, anchor="mm")
            else:
                # Overlapping pocket membership is retained as text but deliberately
                # left uncoloured to avoid assigning the residue to one pocket.
                label = "/".join(pocket.split()[-1].removeprefix("P") for pocket in hits)
                draw.text((pocket_x, yc), label, fill=colors["muted"], font=small, anchor="mm")
        pr = float(prevalence.loc[residue, "PR"] * 100)
        pps = float(prevalence.loc[residue, "PPS"] * 100)
        x_pr = right_inner + (right_outer - right_inner) * pr / xmax
        x_pps = left_inner - (left_inner - left_outer) * pps / xmax
        half_bar = 8 * scale
        if pps > 0:
            draw.rectangle((x_pps, yc - half_bar, left_inner, yc + half_bar), fill=colors["PPS"])
            label = f"{pps:.1f}%"
            bar_width = left_inner - x_pps
            if bar_width >= 55 * scale:
                draw.text((x_pps + 7 * scale, yc), label, fill="white", font=small, anchor="lm")
            else:
                draw.text((x_pps - 6 * scale, yc), label, fill=colors["text"], font=small, anchor="rm")
        if pr > 0:
            draw.rectangle((right_inner, yc - half_bar, x_pr, yc + half_bar), fill=colors["PR"])
            label = f"{pr:.1f}%"
            bar_width = x_pr - right_inner
            if bar_width >= 55 * scale:
                draw.text((x_pr - 7 * scale, yc), label, fill="white", font=small, anchor="rm")
            else:
                draw.text((x_pr + 6 * scale, yc), label, fill=colors["text"], font=small, anchor="lm")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    png_path = OUTPUT_DIR / "pr_pps_residue_hotspot_pocket_map.png"
    image.save(png_path, dpi=(1200, 1200))

    output = prevalence[["PR", "PPS"]].copy()
    output.columns = ["PR_prevalence", "PPS_prevalence"]
    output = output.join(dominant.add_prefix("dominant_"))
    output = output.join(reference.add_prefix("reference_"))
    output = output.join(membership.add_prefix("member_"))
    output.to_csv(OUTPUT_DIR / "pr_pps_residue_hotspot_pocket_map_data.csv")
    return output


def main() -> None:
    residue_counts, typed_counts = summarize_interactions()
    pocket_membership = {name: read_pocket_residues(path) for name, path in POCKET_FILES.items()}
    reference_interactions = read_reference_interactions()
    output = build_figure(
        residue_counts, typed_counts, pocket_membership, reference_interactions
    )
    print(output.to_string(float_format=lambda value: f"{value:.4f}"))


if __name__ == "__main__":
    main()
