import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")


def parse_config(path):
    cfg = {}
    with open(path, encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(f"{path}:{lineno}: expected key=value")
            key, value = line.split("=", 1)
            cfg[key.strip().lower().replace("-", "_")] = value.strip()
    for key in ("prdir", "ppsdir", "outdir"):
        if not cfg.get(key):
            raise ValueError(f"Missing required config key: {key}")
    cfg["last_n"] = int(cfg.get("last_n", "50"))
    cfg["fp_radius"] = int(cfg.get("fp_radius", "2"))
    cfg["fp_size"] = int(cfg.get("fp_size", "16384"))
    cfg["scaffold_sample_size"] = int(cfg.get("scaffold_sample_size", "3000"))
    cfg["pair_plot_sample_size"] = int(cfg.get("pair_plot_sample_size", "200000"))
    cfg["random_state"] = int(cfg.get("random_state", "42"))
    return cfg


def resolve(value, config_path, source=None):
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    if not path.is_absolute():
        path = Path(config_path).resolve().parent / path
    if path.is_dir() and source:
        path = path / f"{source}_scores.csv"
    return path


def load_population(raw_path, source, cache_path):
    raw = pd.read_csv(raw_path, usecols=lambda column: column in {"step", "smiles"})
    membership = pd.read_csv(
        cache_path, usecols=["source", "source_index", "canon_smiles"]
    )
    membership = membership[membership["source"] == source].copy()
    membership["source_index"] = membership["source_index"].astype(int)
    frame = raw.iloc[membership["source_index"].to_numpy()].copy()
    frame["source_index"] = membership["source_index"].to_numpy()
    frame["canon_smiles"] = membership["canon_smiles"].to_numpy()
    frame["source"] = source
    frame["step"] = pd.to_numeric(frame["step"], errors="coerce")
    return frame.dropna(subset=["step", "canon_smiles"]).reset_index(drop=True)


def murcko_smiles(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    if scaffold.GetNumAtoms() == 0:
        return "[NO_SCAFFOLD]"
    return Chem.MolToSmiles(scaffold, canonical=True)


def deterministic_sample(values, size, seed):
    values = sorted(values)
    if len(values) <= size:
        return values
    rng = np.random.default_rng(seed)
    indices = np.sort(rng.choice(len(values), size=size, replace=False))
    return [values[index] for index in indices]


def fingerprint_scaffolds(scaffolds, generator):
    valid_scaffolds = []
    fingerprints = []
    for scaffold in scaffolds:
        molecule = Chem.MolFromSmiles(scaffold)
        if molecule is None or molecule.GetNumAtoms() == 0:
            continue
        valid_scaffolds.append(scaffold)
        fingerprints.append(generator.GetFingerprint(molecule))
    return valid_scaffolds, fingerprints


def pairwise_and_nearest(fingerprints):
    pairwise_blocks = []
    nearest = np.zeros(len(fingerprints), dtype=float)
    for index in range(1, len(fingerprints)):
        similarities = np.asarray(
            DataStructs.BulkTanimotoSimilarity(fingerprints[index], fingerprints[:index]),
            dtype=np.float32,
        )
        pairwise_blocks.append(similarities)
        if similarities.size:
            nearest[index] = max(nearest[index], float(similarities.max()))
            nearest[:index] = np.maximum(nearest[:index], similarities)
    pairwise = (
        np.concatenate(pairwise_blocks).astype(np.float32, copy=False)
        if pairwise_blocks else np.asarray([], dtype=np.float32)
    )
    return pairwise, nearest


def distribution_summary(source, scope, population_count, unique_count,
                         no_scaffold_count, sampled_count, pairwise, nearest):
    result = {
        "source": source,
        "scope": scope,
        "population_molecules": population_count,
        "unique_canonical_scaffolds": unique_count,
        "no_scaffold_molecules": no_scaffold_count,
        "sampled_unique_scaffolds": sampled_count,
        "pair_count": int(pairwise.size),
    }
    for name, values in (("pairwise", pairwise), ("nearest", nearest)):
        result.update({
            f"{name}_mean": float(values.mean()) if values.size else None,
            f"{name}_median": float(np.median(values)) if values.size else None,
            f"{name}_q05": float(np.quantile(values, 0.05)) if values.size else None,
            f"{name}_q25": float(np.quantile(values, 0.25)) if values.size else None,
            f"{name}_q75": float(np.quantile(values, 0.75)) if values.size else None,
            f"{name}_q95": float(np.quantile(values, 0.95)) if values.size else None,
            f"{name}_fraction_ge_0_4": float((values >= 0.4).mean()) if values.size else None,
            f"{name}_fraction_ge_0_6": float((values >= 0.6).mean()) if values.size else None,
            f"{name}_fraction_ge_0_8": float((values >= 0.8).mean()) if values.size else None,
        })
    return result


def ecdf(values):
    ordered = np.sort(values)
    y = np.arange(1, len(ordered) + 1, dtype=float) / len(ordered)
    return ordered, y


def make_figure(plot_values, output_path):
    colors = {
        ("PR", "full"): "#9CC4E4",
        ("PR", "last_window"): "#2878B5",
        ("PPS", "full"): "#F2AAA6",
        ("PPS", "last_window"): "#D9534F",
    }
    labels = {
        ("PR", "full"): "PR full",
        ("PR", "last_window"): "PR final 50",
        ("PPS", "full"): "PPS full",
        ("PPS", "last_window"): "PPS final 50",
    }
    width, height = 2400, 1080
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font_root = Path(r"C:\Windows\Fonts")
    title_font = ImageFont.truetype(str(font_root / "arialbd.ttf"), 38)
    heading_font = ImageFont.truetype(str(font_root / "arialbd.ttf"), 25)
    body_font = ImageFont.truetype(str(font_root / "arial.ttf"), 20)
    small_font = ImageFont.truetype(str(font_root / "arial.ttf"), 18)
    draw.text((width / 2, 35), "Canonical Murcko scaffold fingerprint similarity", fill="#222222", font=title_font, anchor="ma")

    panels = [(120, 160, 1110, 850), (1290, 160, 2280, 850)]
    panel_titles = ["Pairwise scaffold similarity", "Nearest-scaffold similarity"]
    for panel, panel_title in zip(panels, panel_titles):
        left, top, right, bottom = panel
        draw.text(((left + right) / 2, 115), panel_title, fill="#222222", font=heading_font, anchor="ma")
        for tick in np.linspace(0, 1, 6):
            x = left + tick * (right - left)
            y = bottom - tick * (bottom - top)
            draw.line((x, top, x, bottom), fill="#E6E6E6", width=2)
            draw.line((left, y, right, y), fill="#E6E6E6", width=2)
            draw.text((x, bottom + 20), f"{tick:.1f}", fill="#666666", font=small_font, anchor="ma")
            draw.text((left - 18, y), f"{tick:.1f}", fill="#666666", font=small_font, anchor="rm")
        draw.line((left, bottom, right, bottom), fill="#777777", width=3)
        draw.line((left, top, left, bottom), fill="#777777", width=3)
        draw.text(((left + right) / 2, bottom + 52), "Tanimoto similarity", fill="#333333", font=body_font, anchor="ma")

    y_label = Image.new("RGBA", (360, 60), (255, 255, 255, 0))
    y_draw = ImageDraw.Draw(y_label)
    y_draw.text((180, 30), "Cumulative fraction", fill="#333333", font=body_font, anchor="mm")
    y_label = y_label.rotate(90, expand=True)
    y_center = int((panels[0][1] + panels[0][3]) / 2)
    image.paste(y_label, (18, y_center - y_label.height // 2), y_label)

    for key in (("PR", "full"), ("PR", "last_window"),
                ("PPS", "full"), ("PPS", "last_window")):
        for panel, value_name in zip(panels, ("pairwise", "nearest")):
            x_values, y_values = ecdf(plot_values[key][value_name])
            if len(x_values) > 5000:
                indices = np.linspace(0, len(x_values) - 1, 5000).astype(int)
                x_values, y_values = x_values[indices], y_values[indices]
            left, top, right, bottom = panel
            points = [
                (left + float(x) * (right - left), bottom - float(y) * (bottom - top))
                for x, y in zip(x_values, y_values)
            ]
            if len(points) >= 2:
                draw.line(points, fill=colors[key], width=5, joint="curve")

    legend_y = 995
    x_cursor = 620
    for key in (("PR", "full"), ("PR", "last_window"),
                ("PPS", "full"), ("PPS", "last_window")):
        draw.line((x_cursor, legend_y, x_cursor + 42, legend_y), fill=colors[key], width=7)
        draw.text((x_cursor + 52, legend_y), labels[key], fill="#333333", font=body_font, anchor="lm")
        x_cursor += 300
    image.save(output_path, dpi=(300, 300))


def make_combined_summary(scaffold_summary, writeup_root):
    """Combine exact Murcko diversity, scaffold-family, and whole-molecule similarity."""
    aggregate_path = writeup_root / "scaffold_diversity_analysis" / "aggregate.csv"
    if not aggregate_path.exists():
        raise FileNotFoundError(
            f"Required molecular-similarity aggregate was not found: {aggregate_path}"
        )
    molecular = pd.read_csv(aggregate_path)
    required = {
        "source", "scope", "unique_scaffold_fraction", "median_similarity",
        "median_nearest_similarity_within_scope",
    }
    missing = required.difference(molecular.columns)
    if missing:
        raise ValueError(f"Missing columns in {aggregate_path}: {sorted(missing)}")

    scaffold = scaffold_summary[
        ["source", "scope", "pairwise_median", "nearest_median"]
    ].copy()
    combined = molecular[
        [
            "source", "scope", "unique_scaffold_fraction", "median_similarity",
            "median_nearest_similarity_within_scope",
        ]
    ].merge(scaffold, on=["source", "scope"], how="inner", validate="one_to_one")
    combined = combined.rename(columns={
        "unique_scaffold_fraction": "exact_scaffold_diversity",
        "pairwise_median": "scaffold_pairwise_similarity",
        "nearest_median": "scaffold_nearest_similarity",
        "median_similarity": "molecular_pairwise_similarity",
        "median_nearest_similarity_within_scope": "molecular_nearest_similarity",
    })
    scope_label = {"full": "full", "last_window": "final 50"}
    combined["Population"] = (
        combined["source"] + " " + combined["scope"].map(scope_label)
    )
    order = ["PR full", "PR final 50", "PPS full", "PPS final 50"]
    combined["_order"] = pd.Categorical(combined["Population"], order, ordered=True)
    combined = combined.sort_values("_order").drop(columns=["source", "scope", "_order"])
    combined = combined[[
        "Population", "exact_scaffold_diversity",
        "scaffold_pairwise_similarity", "scaffold_nearest_similarity",
        "molecular_pairwise_similarity", "molecular_nearest_similarity",
    ]]

    formatted = pd.DataFrame({
        "Population": combined["Population"],
        "Exact scaffold diversity¹": combined["exact_scaffold_diversity"].map(
            lambda value: f"{value:.1%}"
        ),
        "Scaffold pairwise similarity²": combined["scaffold_pairwise_similarity"].map(
            lambda value: f"{value:.3f}"
        ),
        "Scaffold nearest similarity³": combined["scaffold_nearest_similarity"].map(
            lambda value: f"{value:.3f}"
        ),
        "Molecular pairwise similarity⁴": combined["molecular_pairwise_similarity"].map(
            lambda value: f"{value:.3f}"
        ),
        "Molecular nearest similarity⁵": combined["molecular_nearest_similarity"].map(
            lambda value: f"{value:.3f}"
        ),
    })
    return combined, formatted


def main(config_path):
    cfg = parse_config(config_path)
    writeup_root = resolve(cfg["outdir"], config_path)
    outdir = writeup_root / "scaffold_similarity_analysis"
    outdir.mkdir(parents=True, exist_ok=True)
    cache_path = resolve(cfg["outdir"], config_path) / "embedding_cache" / "joint_umap_fp.csv.gz"
    if not cache_path.exists():
        raise FileNotFoundError(f"Population cache not found: {cache_path}")

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=cfg["fp_radius"], fpSize=cfg["fp_size"]
    )
    summary_rows = []
    nearest_rows = []
    histogram_rows = []
    sampled_scaffold_rows = []
    plot_values = {}
    bins = np.linspace(0, 1, 21)

    for source_index, (source, key) in enumerate((("PR", "prdir"), ("PPS", "ppsdir"))):
        population = load_population(resolve(cfg[key], config_path, source), source, cache_path)
        population["murcko_scaffold"] = population["canon_smiles"].map(murcko_smiles)
        population = population[population["murcko_scaffold"].notna()].copy()
        max_step = int(population["step"].max())
        first_last_step = max(int(population["step"].min()), max_step - cfg["last_n"] + 1)
        recent = population[population["step"].between(first_last_step, max_step)].copy()

        for scope_index, (scope, frame) in enumerate((("full", population), ("last_window", recent))):
            no_scaffold_count = int((frame["murcko_scaffold"] == "[NO_SCAFFOLD]").sum())
            unique_scaffolds = sorted(
                set(frame.loc[frame["murcko_scaffold"] != "[NO_SCAFFOLD]", "murcko_scaffold"])
            )
            seed = cfg["random_state"] + source_index * 10 + scope_index
            sampled = deterministic_sample(
                unique_scaffolds, cfg["scaffold_sample_size"], seed
            )
            sampled, fingerprints = fingerprint_scaffolds(sampled, generator)
            pairwise, nearest = pairwise_and_nearest(fingerprints)

            summary_rows.append(distribution_summary(
                source, scope, len(frame), len(unique_scaffolds),
                no_scaffold_count, len(sampled), pairwise, nearest,
            ))
            nearest_rows.extend({
                "source": source,
                "scope": scope,
                "murcko_scaffold": scaffold,
                "nearest_similarity": float(similarity),
            } for scaffold, similarity in zip(sampled, nearest))
            sampled_scaffold_rows.extend({
                "source": source,
                "scope": scope,
                "murcko_scaffold": scaffold,
            } for scaffold in sampled)

            counts, _ = np.histogram(pairwise, bins=bins)
            histogram_rows.extend({
                "source": source,
                "scope": scope,
                "bin_left": float(bins[index]),
                "bin_right": float(bins[index + 1]),
                "pair_count": int(count),
                "pair_fraction": float(count / pairwise.size) if pairwise.size else None,
            } for index, count in enumerate(counts))

            rng = np.random.default_rng(seed + 1000)
            if pairwise.size > cfg["pair_plot_sample_size"]:
                indices = rng.choice(pairwise.size, cfg["pair_plot_sample_size"], replace=False)
                pair_plot = pairwise[indices]
            else:
                pair_plot = pairwise
            plot_values[(source, scope)] = {
                "pairwise": pair_plot,
                "nearest": nearest,
            }
            print(
                f"[{source} {scope}] molecules={len(frame)} unique_scaffolds={len(unique_scaffolds)} "
                f"sampled={len(sampled)} pairs={pairwise.size}"
            )

    summary = pd.DataFrame(summary_rows)
    combined_summary, combined_formatted = make_combined_summary(summary, writeup_root)
    nearest_table = pd.DataFrame(nearest_rows)
    histogram = pd.DataFrame(histogram_rows)
    sampled_table = pd.DataFrame(sampled_scaffold_rows)
    summary.to_csv(outdir / "scaffold_similarity_summary.csv", index=False)
    combined_summary.to_csv(outdir / "combined_similarity_summary.csv", index=False)
    combined_formatted.to_csv(
        outdir / "combined_similarity_summary_formatted.csv", index=False,
        encoding="utf-8-sig",
    )
    nearest_table.to_csv(outdir / "nearest_scaffold_similarity.csv", index=False)
    histogram.to_csv(outdir / "pairwise_scaffold_similarity_histogram.csv", index=False)
    sampled_table.to_csv(outdir / "sampled_unique_scaffolds.csv", index=False)
    make_figure(plot_values, outdir / "scaffold_similarity_ecdf.png")
    (outdir / "manifest.json").write_text(json.dumps({
        "scientific_question": (
            "Are formally distinct canonical Murcko scaffolds genuinely diverse, "
            "or structurally similar variants of related scaffold families?"
        ),
        "analysis_unit": "Unique canonical Bemis-Murcko scaffold (unweighted)",
        "fingerprint": {
            "type": "Morgan bit fingerprint on the isolated scaffold",
            "radius": cfg["fp_radius"],
            "fp_size": cfg["fp_size"],
            "similarity": "Tanimoto",
        },
        "last_n": cfg["last_n"],
        "scaffold_sample_size": cfg["scaffold_sample_size"],
        "pair_plot_sample_size": cfg["pair_plot_sample_size"],
        "random_state": cfg["random_state"],
        "no_scaffold_handling": "Acyclic [NO_SCAFFOLD] records counted but excluded from similarity",
    }, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare fingerprint similarity among unique canonical Murcko scaffolds."
    )
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)
