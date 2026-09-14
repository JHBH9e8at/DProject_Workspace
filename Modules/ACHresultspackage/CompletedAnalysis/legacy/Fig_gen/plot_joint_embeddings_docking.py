import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import pandas as pd


COLORS = {"PR": "#2878B5", "PPS": "#D9534F"}
REF_COLOR = "#7F3C8D"
SCORE_CMAP = "viridis"


def parse_config(path):
    cfg = {}
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line and not line.startswith("#") and "=" in line:
                key, value = line.split("=", 1)
                cfg[key.strip().lower().replace("-", "_")] = value.strip()
    if not cfg.get("outdir"):
        raise ValueError("Missing required config key: outdir")
    cfg.setdefault("umap_mode", "on")
    cfg.setdefault("tsne_mode", "on")
    cfg.setdefault("descriptormode", cfg.get("descriptor_mode", "both"))
    cfg.setdefault("featuremode", cfg.get("feature_mode", "both"))
    return cfg


def resolve_outdir(value, config_path):
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    if not path.is_absolute():
        path = Path(config_path).resolve().parent / path
    return path


def requested_caches(cfg, cache_dir):
    methods = [
        method
        for method in ("umap", "tsne")
        if cfg[f"{method}_mode"].lower() == "on"
    ]
    features = []
    if cfg["featuremode"].lower() in ("fp", "both"):
        features.append("fp")
    if cfg["featuremode"].lower() in ("desc", "both"):
        if cfg["descriptormode"].lower() in ("all", "both"):
            features.append("all_desc")
        if cfg["descriptormode"].lower() in ("select", "both"):
            features.append("select_desc")
    return [
        cache_dir / f"joint_{method}_{feature}.csv.gz"
        for method in methods
        for feature in features
    ]


def plot_cache(path, figure_dir):
    df = pd.read_csv(path)
    df["docking_score"] = pd.to_numeric(df["docking_score"], errors="coerce")
    references = df[df["source"] == "REF"].copy()
    if references.empty:
        raise ValueError(
            f"Embedding cache contains no reference compound: {path}. "
            "Re-run build_joint_embeddings.py to rebuild the cache."
        )
    population = df[
        df["source"].isin(["PR", "PPS"])
    ].dropna(subset=["docking_score"]).copy()
    method = df["method"].iloc[0].upper()
    feature = df["feature"].iloc[0]
    score_norm = Normalize(
        vmin=population["docking_score"].min(),
        vmax=population["docking_score"].max(),
    )
    fig, axes = plt.subplots(
        1, 3, figsize=(19, 7), sharex=True, sharey=True,
        constrained_layout=True,
    )

    for source in ("PR", "PPS"):
        part = population[population["source"] == source]
        axes[0].scatter(
            part["embedding_x"], part["embedding_y"], s=5, alpha=0.35,
            color=COLORS[source], linewidths=0, rasterized=True, label=source,
        )
    for _, ref in references.iterrows():
        label = ref.get("label", "Reference")
        for ax in axes:
            ax.scatter(
                ref["embedding_x"], ref["embedding_y"], marker="*", s=190,
                color=REF_COLOR, edgecolor="white", linewidth=0.9,
                zorder=10, label=label if ax is axes[0] else "_nolegend_",
            )
            ax.annotate(
                label, (ref["embedding_x"], ref["embedding_y"]),
                xytext=(7, 6), textcoords="offset points", fontsize=10,
                fontweight="bold", zorder=11,
            )
    axes[0].set_title("PR + PPS", fontsize=16, fontweight="normal")
    axes[0].legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=3,
        frameon=False,
        fontsize=12,
        markerscale=2.5,
    )

    for ax, source in zip(axes[1:], ("PR", "PPS")):
        part = population[population["source"] == source].sort_values(
            "docking_score", ascending=False
        )
        ax.scatter(
            part["embedding_x"], part["embedding_y"],
            c=part["docking_score"],
            cmap=SCORE_CMAP,
            norm=score_norm,
            s=6, alpha=0.55, linewidths=0, rasterized=True,
        )
        ax.set_title(source, fontsize=16, fontweight="normal")

        score_mappable = ScalarMappable(
            norm=score_norm, cmap=SCORE_CMAP
        )
        score_mappable.set_array([])
        colorbar = fig.colorbar(
            score_mappable,
            ax=ax,
            orientation="horizontal",
            location="bottom",
            shrink=0.82,
            pad=-0.1,
        )
        colorbar.set_label("Docking score", labelpad=8)

    axes[0].set_xlabel(f"{method} 1")
    axes[0].set_ylabel(f"{method} 2")
    for ax in axes:
        ax.grid(alpha=0.1)
    fig.suptitle(
        f"Joint PR/PPS Chemical Space by Docking Score — "
        f"{method} ({feature})",
        fontsize=15,
    )
    figure_dir.mkdir(parents=True, exist_ok=True)
    output = figure_dir / f"joint_{method.lower()}_{feature}_docking.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {output}")

def main(config_path):
    cfg = parse_config(config_path)
    outdir = resolve_outdir(cfg["outdir"], config_path)
    cache_dir = outdir / "embedding_cache"
    figure_dir = outdir / "joint_embedding_docking_figures"
    for cache in requested_caches(cfg, cache_dir):
        if not cache.exists():
            raise FileNotFoundError(
                f"Embedding cache not found: {cache}. "
                "Run build_joint_embeddings.py first."
            )
        plot_cache(cache, figure_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot cached joint PR/PPS UMAP/t-SNE maps by docking score."
    )
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)
