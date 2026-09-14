from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

try:
    import generate_residue_hotspot_figure as base
except ModuleNotFoundError:
    sys.path.insert(0, r"Q:\coding_dir\701_Project\workspace\Modules\ACHresultspackage\Fig_gen\reshotstopt")
    import generate_residue_hotspot_figure as base

import generate_crossdock_own_opposite_hotspots as crossdock


OUTPUT_STEM = "crossdock_own_opposite_residue_hotspots_fullstyle"


def summarize_own_and_opposite(config_path: Path):
    cfg = crossdock.parse_config(config_path)
    own_root = crossdock.resolve(cfg["own_result_dir"], config_path) / "aggregate"
    opposite_root = crossdock.resolve(cfg["opposite_result_dir"], config_path) / "aggregate"
    cohorts = crossdock.load_cohorts(opposite_root / "pose_metadata.csv")
    selected = crossdock.select_best_own_poses(own_root / "pose_metadata.csv", cohorts)
    denominators_by_arm = {
        arm: len(cohorts[arm].intersection(set(selected.loc[selected["crossdock_arm"].eq(arm), "molecule_id"])))
        for arm in crossdock.ARM_SPECS
    }
    own = crossdock.summarize_own(own_root / "interactions_long.csv", selected, denominators_by_arm)
    opposite = crossdock.summarize_opposite(opposite_root / "interactions_long.csv", denominators_by_arm)
    arm_to_state = {arm: spec["ligand"] for arm, spec in crossdock.ARM_SPECS.items()}
    for frame in (own, opposite):
        frame["population_state"] = frame["crossdock_arm"].map(arm_to_state)
        frame["poses"] = frame["molecule_count"]
    denominators = {arm_to_state[arm]: count for arm, count in denominators_by_arm.items()}
    # Derive interaction-type annotations from the matched own-state poses,
    # rather than borrowing labels from the complete AHC population.
    interaction_columns = [
        "population_state", "pose_index", "molecule_id", "residue_name",
        "original_residue_number", "interaction_type",
    ]
    own_interactions = pd.read_csv(
        own_root / "interactions_long.csv", usecols=interaction_columns
    )
    own_interactions["molecule_id"] = own_interactions["molecule_id"].astype(str)
    typed_parts = []
    for arm, spec in crossdock.ARM_SPECS.items():
        keys = selected.loc[
            selected["crossdock_arm"].eq(arm),
            ["population_state", "pose_index", "molecule_id"],
        ]
        part = own_interactions.loc[
            own_interactions["population_state"].eq(spec["ligand"])
        ].merge(
            keys,
            on=["population_state", "pose_index", "molecule_id"],
            how="inner",
            validate="many_to_one",
        )
        part["residue"] = crossdock.residue_column(part)
        part = part.dropna(subset=["residue", "interaction_type"]).drop_duplicates(
            ["population_state", "molecule_id", "residue", "interaction_type"]
        )
        typed_parts.append(part)
    typed = pd.concat(typed_parts, ignore_index=True)
    typed_counts = (
        typed.groupby(["residue", "population_state", "interaction_type"])["molecule_id"]
        .nunique().rename("poses").reset_index()
    )
    typed_counts["prevalence"] = typed_counts.apply(
        lambda row: row["poses"] / denominators[row["population_state"]], axis=1
    )
    return cfg, own, opposite, typed_counts, denominators


def residue_sort_key(residue: str) -> tuple[int, str]:
    digits = "".join(character for character in residue if character.isdigit())
    name = "".join(character for character in residue if not character.isdigit())
    return int(digits), name


def prepare_wide(summary: pd.DataFrame, residues: list[str]) -> pd.DataFrame:
    wide = summary.pivot(index="residue", columns="population_state", values="prevalence")
    wide = wide.reindex(residues).fillna(0.0)
    for state in ("PPS", "PR"):
        if state not in wide:
            wide[state] = 0.0
    return wide


def draw_figure(
    full_summary: pd.DataFrame,
    last_summary: pd.DataFrame,
    typed_counts: pd.DataFrame,
    pocket_membership: dict[str, set[str]],
    reference: pd.DataFrame,
    denominators: dict[str, int],
    output_dir: Path,
    figure_title: str = "PR/PPS residue interaction hotspots: own vs opposite receptor state",
) -> pd.DataFrame:
    candidate = pd.concat([full_summary, last_summary], ignore_index=True)
    candidate_wide = candidate.groupby(["residue", "population_state"])["prevalence"].max().unstack(fill_value=0)
    selected = set(candidate_wide.index[candidate_wide.max(axis=1) >= base.MIN_PREVALENCE])
    for state in ("PPS", "PR"):
        selected.update(candidate_wide[state].nlargest(base.TOP_N_PER_STATE).index)
    selected.update(base.REFERENCE_RESIDUES.intersection(candidate_wide.index))
    residues = sorted(selected, key=residue_sort_key)

    full = prepare_wide(full_summary, residues)
    last = prepare_wide(last_summary, residues)
    dominant = base.dominant_interaction_labels(typed_counts).pivot(
        index="residue", columns="population_state", values="dominant_type"
    ).reindex(residues).fillna("-")
    reference = reference.reindex(residues).fillna("-")
    membership = pd.DataFrame(
        {name: [residue in members for residue in residues] for name, members in pocket_membership.items()},
        index=residues,
    )

    scale = 2
    width = 1800 * scale
    left_outer, left_inner = 130 * scale, 650 * scale
    backbone_left, backbone_right = 650 * scale, 1150 * scale
    right_inner, right_outer = 1150 * scale, 1670 * scale
    # Match generate_residue_hotspot_figure_full_vs_last50.py.
    top, row_h, bottom = 150 * scale, 32 * scale, 95 * scale
    height = top + len(residues) * row_h + bottom
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    fonts = Path(r"C:\Windows\Fonts")
    regular = ImageFont.truetype(str(fonts / "arial.ttf"), 13 * scale)
    small = ImageFont.truetype(str(fonts / "arial.ttf"), 9 * scale)
    bold = ImageFont.truetype(str(fonts / "arialbd.ttf"), 13 * scale)
    title = ImageFont.truetype(str(fonts / "arialbd.ttf"), 19 * scale)
    color = {
        "PPS_full": "#D9534F", "PPS_last": "#9A9A9A",
        "PR_full": "#2878B5", "PR_last": "#9A9A9A",
        "grid": "#E7E7E7", "text": "#222222", "muted": "#666666",
    }

    draw.text((width / 2, 30 * scale), figure_title, fill=color["text"], font=title, anchor="ma")
    panel_heading_y = 82 * scale
    draw.text(((left_outer + left_inner) / 2, panel_heading_y), "PPS prevalence", fill="#D9534F", font=bold, anchor="ma")
    draw.text(((backbone_left + backbone_right) / 2, panel_heading_y), "Residue", fill=color["text"], font=bold, anchor="ma")
    draw.text(((right_inner + right_outer) / 2, panel_heading_y), "PR prevalence", fill="#2878B5", font=bold, anchor="ma")

    pps_interaction_x = backbone_left + 35 * scale
    pps_pocket_x = backbone_left + 105 * scale
    pps_reference_x = backbone_left + 185 * scale
    residue_x = (backbone_left + backbone_right) / 2
    pr_reference_x = backbone_right - 185 * scale
    pr_pocket_x = backbone_right - 105 * scale
    pr_interaction_x = backbone_right - 35 * scale
    for x, label in (
        (pps_interaction_x, "Interaction"), (pps_pocket_x, "Pocket"),
        (pps_reference_x, "Reference"), (pr_reference_x, "Reference"),
        (pr_pocket_x, "Pocket"), (pr_interaction_x, "Interaction"),
    ):
        draw.text((x, 118 * scale), label, fill=color["muted"], font=small, anchor="mm")

    xmax = max(80, int(np.ceil(max(full.to_numpy().max(), last.to_numpy().max()) * 100 / 10) * 10))
    for tick in range(0, xmax + 1, 10):
        xl = left_inner - int((left_inner - left_outer) * tick / xmax)
        xr = right_inner + int((right_outer - right_inner) * tick / xmax)
        draw.line((xl, top - 8 * scale, xl, top + len(residues) * row_h), fill=color["grid"], width=scale)
        draw.line((xr, top - 8 * scale, xr, top + len(residues) * row_h), fill=color["grid"], width=scale)
        draw.text((xl, top - 28 * scale), str(tick), fill=color["muted"], font=small, anchor="mm")
        draw.text((xr, top - 28 * scale), str(tick), fill=color["muted"], font=small, anchor="mm")
    draw.line((left_inner, top - 8 * scale, left_inner, top + len(residues) * row_h), fill="#999999", width=2 * scale)
    draw.line((right_inner, top - 8 * scale, right_inner, top + len(residues) * row_h), fill="#999999", width=2 * scale)

    annotation_top = 100 * scale
    annotation_bottom = top + len(residues) * row_h
    annotation_outline = "#B8B8B8"
    residue_clearance = 30 * scale
    outer_padding = 45 * scale
    draw.rectangle(
        (left_outer - outer_padding, annotation_top,
         residue_x - residue_clearance, annotation_bottom),
        outline=annotation_outline, width=2 * scale,
    )
    draw.rectangle(
        (residue_x + residue_clearance, annotation_top,
         right_outer + outer_padding, annotation_bottom),
        outline=annotation_outline, width=2 * scale,
    )
    draw.rectangle(
        (residue_x - residue_clearance, annotation_top,
         residue_x + residue_clearance, annotation_bottom),
        outline=annotation_outline, width=2 * scale,
    )

    pocket_legend_y = annotation_bottom + 22 * scale
    pocket_legend = [
        ("PPS P2", pps_pocket_x - 75 * scale),
        ("PPS P38", pps_pocket_x - 25 * scale),
        ("PPS P51", pps_pocket_x + 30 * scale),
        ("PR P11", pr_pocket_x - 75 * scale),
        ("PR P15", pr_pocket_x - 25 * scale),
        ("PR P60", pr_pocket_x + 30 * scale),
    ]
    for pocket, x in pocket_legend:
        draw.rectangle(
            (x, pocket_legend_y, x + 11 * scale, pocket_legend_y + 11 * scale),
            fill=base.POCKET_COLORS[pocket],
        )
        draw.text(
            (x + 15 * scale, pocket_legend_y + 5 * scale), pocket.split()[-1],
            fill=color["muted"], font=small, anchor="lm",
        )

    population_legend_y = annotation_bottom + 22 * scale
    population_legend = [
        ("Own state", color["PPS_full"], 300 * scale),
        ("Opposite state", color["PPS_last"], 405 * scale),
        ("Own state", color["PR_full"], 1320 * scale),
        ("Opposite state", color["PR_last"], 1425 * scale),
    ]
    for label, fill, x in population_legend:
        draw.rectangle(
            (x, population_legend_y, x + 11 * scale, population_legend_y + 11 * scale),
            fill=fill,
        )
        draw.text(
            (x + 15 * scale, population_legend_y + 5 * scale), label,
            fill=color["muted"], font=small, anchor="lm",
        )

    pocket_names = list(pocket_membership)
    for row, residue in enumerate(residues):
        yc = top + row * row_h + row_h / 2
        draw.text((residue_x, yc), residue, fill=color["text"], font=regular, anchor="mm")
        pps_interaction = (
            dominant.loc[residue, "PPS"]
            if round(max(full.loc[residue, "PPS"], last.loc[residue, "PPS"]) * 100, 1) > 0
            else ""
        )
        pr_interaction = (
            dominant.loc[residue, "PR"]
            if round(max(full.loc[residue, "PR"], last.loc[residue, "PR"]) * 100, 1) > 0
            else ""
        )
        draw.text((pps_interaction_x, yc), pps_interaction, fill=color["muted"], font=small, anchor="mm")
        draw.text((pr_interaction_x, yc), pr_interaction, fill=color["muted"], font=small, anchor="mm")
        pps_reference = "" if reference.loc[residue, "PPS"] == "-" else reference.loc[residue, "PPS"]
        pr_reference = "" if reference.loc[residue, "PR"] == "-" else reference.loc[residue, "PR"]
        draw.text((pps_reference_x, yc), pps_reference, fill=color["text"], font=small, anchor="mm")
        draw.text((pr_reference_x, yc), pr_reference, fill=color["text"], font=small, anchor="mm")

        for state, px in (("PPS", pps_pocket_x), ("PR", pr_pocket_x)):
            hits = [name for name in pocket_names if name.startswith(state + " ") and membership.loc[residue, name]]
            if len(hits) == 1:
                pocket = hits[0]
                draw.rectangle((px - 24 * scale, yc - 10 * scale, px + 24 * scale, yc + 10 * scale), fill=base.POCKET_COLORS[pocket])
                draw.text((px, yc), pocket.split()[-1], fill="white", font=small, anchor="mm")
            else:
                if hits:
                    label = "/".join(item.split()[-1].removeprefix("P") for item in hits)
                    draw.text((px, yc), label, fill=color["muted"], font=small, anchor="mm")

        for state, inner, outer, direction in (
            ("PPS", left_inner, left_outer, -1), ("PR", right_inner, right_outer, 1),
        ):
            full_value = float(full.loc[residue, state] * 100)
            last_value = float(last.loc[residue, state] * 100)
            span = abs(outer - inner)
            x_full = inner + direction * span * full_value / xmax
            x_last = inner + direction * span * last_value / xmax
            full_fill = color[f"{state}_full"]
            last_fill = color[f"{state}_last"]
            y_full = yc - 6 * scale
            y_last = yc + 6 * scale
            half_bar = 4 * scale
            draw.rectangle((min(inner, x_full), y_full - half_bar, max(inner, x_full), y_full + half_bar), fill=full_fill)
            draw.rectangle((min(inner, x_last), y_last - half_bar, max(inner, x_last), y_last + half_bar), fill=last_fill)
            if round(full_value, 1) > 0:
                anchor = "rm" if state == "PPS" else "lm"
                offset = -5 * scale if state == "PPS" else 5 * scale
                draw.text((x_full + offset, y_full), f"{full_value:.1f}", fill=color["muted"], font=small, anchor=anchor)
            if round(last_value, 1) > 0:
                anchor = "rm" if state == "PPS" else "lm"
                offset = -5 * scale if state == "PPS" else 5 * scale
                draw.text((x_last + offset, y_last), f"{last_value:.1f}", fill=color["text"], font=small, anchor=anchor)

    # draw.text((width / 2, height - 55 * scale), "Prevalence (%): upper pale bars = full population; lower dark bars = last 50 iterations.", fill=color["muted"], font=small, anchor="ma")
    output_dir.mkdir(parents=True, exist_ok=True)
    image.save(output_dir / f"{OUTPUT_STEM}.png", dpi=(1200, 1200))

    output = pd.concat(
        [full.add_prefix("own_"), last.add_prefix("opposite_")], axis=1
    ).join(dominant.add_prefix("dominant_own_")).join(reference.add_prefix("reference_")).join(membership.add_prefix("member_"))
    output.to_csv(output_dir / f"{OUTPUT_STEM}_data.csv")
    return output


def main(config_path: Path) -> None:
    cfg, full_summary, last_summary, typed_counts, denominators = summarize_own_and_opposite(config_path)
    pocket_root = crossdock.resolve(cfg["fpocket_root"], config_path)
    pockets = {
        f"{state} {label}": crossdock.read_pocket_residues(pocket_root / relative)
        for state, specs in crossdock.POCKET_SPECS.items()
        for label, relative, _ in specs
    }
    reference = base.read_reference_interactions()
    output_dir = crossdock.resolve(cfg["output_dir"], config_path)
    output = draw_figure(full_summary, last_summary, typed_counts, pockets, reference, denominators, output_dir)
    print(f"Matched denominators: PPS={denominators['PPS']}, PR={denominators['PR']}")
    print(output.to_string(float_format=lambda value: f"{value:.4f}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    main(args.config)
