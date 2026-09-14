from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

from Modules.Analysis_Figen_block.common.output_naming import allocate_versioned_group


PPS = "#D9534F"
PR = "#2878B5"
NEUTRAL = "#9A9A9A"
GRID = "#E2E2E2"
TEXT = "#262626"


def rgba(hex_color: str, alpha: int) -> tuple[int, int, int, int]:
    value = hex_color.lstrip("#")
    return tuple(int(value[i:i + 2], 16) for i in (0, 2, 4)) + (alpha,)


def load_rows(path: Path) -> dict[str, list[dict]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    groups: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        for column in (
            "own_state_score", "opposite_state_score", "selectivity_margin",
            "ligand_efficiency_margin",
        ):
            row[column] = float(row[column])
        groups[row["crossdock_arm"]].append(row)
    for arm in ("L_PPS_to_R_PR", "L_PR_to_R_PPS"):
        if not groups[arm]:
            raise ValueError(f"Missing expected cross-docking arm: {arm}")
    return groups


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def font(name: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(Path(r"C:\Windows\Fonts") / name), size)


def draw_axes(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int],
              y_min: float, y_max: float, body_font, y_label: bool = True) -> None:
    left, top, right, bottom = box
    for tick in range(math.ceil(y_min), math.floor(y_max) + 1, 2):
        y = bottom - (tick - y_min) / (y_max - y_min) * (bottom - top)
        draw.line((left, y, right, y), fill=GRID, width=2)
        draw.text((left - 18, y), f"{tick}", font=body_font, fill=TEXT, anchor="rm")
    draw.line((left, top, left, bottom), fill="#666666", width=2)
    draw.line((left, bottom, right, bottom), fill="#666666", width=2)
    if y_label:
        label = Image.new("RGBA", (420, 55), (255, 255, 255, 0))
        label_draw = ImageDraw.Draw(label)
        label_draw.text((210, 27), "Docking score (kcal/mol)", font=body_font,
                        fill=TEXT, anchor="mm")
        label = label.rotate(90, expand=True)
        image.paste(label, (left - 100, int((top + bottom) / 2 - label.height / 2)), label)


def paired_panel(draw: ImageDraw.ImageDraw, overlay: Image.Image,
                 box: tuple[int, int, int, int], rows: list[dict],
                 own_label: str, opposite_label: str, own_color: str,
                 opposite_color: str, title: str, y_min: float, y_max: float,
                 heading_font, body_font, small_font) -> None:
    left, top, right, bottom = box
    draw_axes(draw, box, y_min, y_max, body_font)
    x_own = left + int((right - left) * 0.27)
    x_opp = left + int((right - left) * 0.73)

    def y_map(value: float) -> float:
        return bottom - (value - y_min) / (y_max - y_min) * (bottom - top)

    layer = ImageDraw.Draw(overlay, "RGBA")
    rng = random.Random(42)
    for row in rows:
        y1, y2 = y_map(row["own_state_score"]), y_map(row["opposite_state_score"])
        selective = row["selectivity_class"] == "selective_own_state"
        line_color = rgba(own_color if selective else NEUTRAL, 52 if selective else 30)
        layer.line((x_own, y1, x_opp, y2), fill=line_color, width=2)
        jitter_own = rng.uniform(-8, 8)
        jitter_opp = rng.uniform(-8, 8)
        layer.ellipse((x_own - 5 + jitter_own, y1 - 5,
                       x_own + 5 + jitter_own, y1 + 5),
                      fill=rgba(own_color, 120))
        layer.ellipse((x_opp - 5 + jitter_opp, y2 - 5,
                       x_opp + 5 + jitter_opp, y2 + 5),
                      fill=rgba(opposite_color, 120))

    own_median = statistics.median(row["own_state_score"] for row in rows)
    opposite_median = statistics.median(row["opposite_state_score"] for row in rows)
    y1, y2 = y_map(own_median), y_map(opposite_median)
    draw.line((x_own, y1, x_opp, y2), fill="#222222", width=6)
    for x, y, color in ((x_own, y1, own_color), (x_opp, y2, opposite_color)):
        draw.ellipse((x - 12, y - 12, x + 12, y + 12), fill=color,
                     outline="white", width=3)
    draw.text((x_own - 18, y1), f"{own_median:.2f}", font=small_font,
              fill=TEXT, anchor="rm")
    draw.text((x_opp + 18, y2), f"{opposite_median:.2f}", font=small_font,
              fill=TEXT, anchor="lm")
    draw.text(((left + right) / 2, top - 70), title, font=heading_font,
              fill=TEXT, anchor="mm")
    selected = sum(row["selectivity_class"] == "selective_own_state" for row in rows)
    draw.text(((left + right) / 2, top - 30),
              f"Selective: {selected}/{len(rows)} ({selected / len(rows):.1%})",
              font=small_font, fill=TEXT, anchor="mm")
    draw.text((x_own, bottom + 38), own_label, font=body_font, fill=TEXT, anchor="ma")
    draw.text((x_opp, bottom + 38), opposite_label, font=body_font, fill=TEXT, anchor="ma")


def margin_panel(draw: ImageDraw.ImageDraw, overlay: Image.Image,
                 box: tuple[int, int, int, int], groups: list[tuple[str, list[dict], str]],
                 heading_font, body_font, small_font) -> None:
    left, top, right, bottom = box
    all_values = [row["selectivity_margin"] for _, rows, _ in groups for row in rows]
    x_min = math.floor(min(all_values) - 0.5)
    x_max = math.ceil(max(all_values) + 0.5)

    def x_map(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * (right - left)

    for tick in range(x_min, x_max + 1, 2):
        x = x_map(tick)
        draw.line((x, top, x, bottom), fill=GRID, width=2)
        draw.text((x, bottom + 24), str(tick), font=small_font, fill=TEXT, anchor="ma")
    draw.line((left, bottom, right, bottom), fill="#666666", width=2)
    for value, dash in ((0, False), (2, True)):
        x = x_map(value)
        if dash:
            for y in range(top, bottom, 18):
                draw.line((x, y, x, min(y + 9, bottom)), fill="#555555", width=3)
        else:
            draw.line((x, top, x, bottom), fill="#555555", width=2)
    draw.text((x_map(2), top - 18), "Selectivity threshold = 2", font=small_font,
              fill=TEXT, anchor="ms")

    layer = ImageDraw.Draw(overlay, "RGBA")
    y_positions = [top + int((bottom - top) * 0.32), top + int((bottom - top) * 0.72)]
    rng = random.Random(2026)
    for (label, rows, color), center_y in zip(groups, y_positions):
        values = [row["selectivity_margin"] for row in rows]
        bins = 30
        counts = [0] * bins
        for value in values:
            index = min(bins - 1, max(0, int((value - x_min) / (x_max - x_min) * bins)))
            counts[index] += 1
        max_count = max(counts)
        points_top, points_bottom = [], []
        for index, count in enumerate(counts):
            x = x_map(x_min + (index + 0.5) / bins * (x_max - x_min))
            half = 72 * count / max_count
            points_top.append((x, center_y - half))
            points_bottom.append((x, center_y + half))
        polygon = points_top + list(reversed(points_bottom))
        layer.polygon(polygon, fill=rgba(color, 48))
        layer.line(points_top, fill=rgba(color, 180), width=3)
        layer.line(points_bottom, fill=rgba(color, 180), width=3)
        for value in values:
            x, y = x_map(value), center_y + rng.uniform(-34, 34)
            layer.ellipse((x - 4, y - 4, x + 4, y + 4), fill=rgba(color, 105))
        median = statistics.median(values)
        q25, q75 = quantile(values, 0.25), quantile(values, 0.75)
        draw.line((x_map(q25), center_y, x_map(q75), center_y), fill="#222222", width=9)
        draw.ellipse((x_map(median) - 11, center_y - 11, x_map(median) + 11, center_y + 11),
                     fill=color, outline="white", width=3)
        draw.text((left - 22, center_y), label, font=body_font, fill=TEXT, anchor="rm")
        draw.text((x_map(median), center_y - 92), f"median {median:.2f}",
                  font=small_font, fill=TEXT, anchor="ms")
    draw.text(((left + right) / 2, top - 65),
              "Distribution of paired state-selectivity margins",
              font=heading_font, fill=TEXT, anchor="mm")
    draw.text(((left + right) / 2, bottom + 62),
              "Selectivity margin: opposite-state score − own-state score (kcal/mol)",
              font=body_font, fill=TEXT, anchor="ma")


def write_summary(groups: dict[str, list[dict]], path: Path) -> None:
    fields = [
        "crossdock_arm", "molecule_count", "own_state_score_median",
        "opposite_state_score_median", "selectivity_margin_median",
        "ligand_efficiency_margin_median", "selective_own_state_count",
        "selective_own_state_fraction",
    ]
    output = []
    for arm in ("L_PPS_to_R_PR", "L_PR_to_R_PPS"):
        rows = groups[arm]
        selected = sum(row["selectivity_class"] == "selective_own_state" for row in rows)
        output.append({
            "crossdock_arm": arm,
            "molecule_count": len(rows),
            "own_state_score_median": statistics.median(row["own_state_score"] for row in rows),
            "opposite_state_score_median": statistics.median(row["opposite_state_score"] for row in rows),
            "selectivity_margin_median": statistics.median(row["selectivity_margin"] for row in rows),
            "ligand_efficiency_margin_median": statistics.median(row["ligand_efficiency_margin"] for row in rows),
            "selective_own_state_count": selected,
            "selective_own_state_fraction": selected / len(rows),
        })
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)


def main(input_csv: Path, output_dir: Path) -> None:
    groups = load_rows(input_csv)
    pps, pr = groups["L_PPS_to_R_PR"], groups["L_PR_to_R_PPS"]
    all_scores = [row[column] for rows in (pps, pr) for row in rows
                  for column in ("own_state_score", "opposite_state_score")]
    y_min = math.floor(min(all_scores) - 0.5)
    y_max = math.ceil(max(all_scores) + 0.5)

    global image
    image = Image.new("RGBA", (2400, 1800), "white")
    draw = ImageDraw.Draw(image)
    overlay = Image.new("RGBA", image.size, (255, 255, 255, 0))
    title_font = font("arialbd.ttf", 48)
    heading_font = font("arialbd.ttf", 31)
    body_font = font("arial.ttf", 25)
    small_font = font("arial.ttf", 22)

    draw.text((1200, 65), "Cross-docking reveals asymmetric receptor-state compatibility",
              font=title_font, fill=TEXT, anchor="ma")
    paired_panel(draw, overlay, (170, 260, 1120, 920), pps,
                 "Own: PPS", "Opposite: PR", PPS, PR,
                 "PPS-derived compounds", y_min, y_max,
                 heading_font, body_font, small_font)
    paired_panel(draw, overlay, (1320, 260, 2270, 920), pr,
                 "Own: PR", "Opposite: PPS", PR, PPS,
                 "PR-derived compounds", y_min, y_max,
                 heading_font, body_font, small_font)
    margin_panel(draw, overlay, (220, 1170, 2250, 1650),
                 [("PPS-derived", pps, PPS), ("PR-derived", pr, PR)],
                 heading_font, body_font, small_font)
    image = Image.alpha_composite(image, overlay).convert("RGB")
    image_path, summary_path = allocate_versioned_group(
        output_dir,
        ("crossdocking_selectivity.png", "crossdocking_selectivity_summary.csv"),
    )
    image.save(image_path, dpi=(300, 300))
    write_summary(groups, summary_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render cross-docking selectivity figure")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    main(args.input, args.output_dir)
