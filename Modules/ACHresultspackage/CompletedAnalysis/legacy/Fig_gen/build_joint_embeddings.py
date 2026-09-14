import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")

# Increment when the cached row composition or embedding inputs change.
# Version 3 fits UMAP and descriptor scaling on the generated PR/PPS
# populations only, then projects the configured reference with transform().
# Scikit-learn t-SNE has no transform(), so its reference remains part of the
# joint fit.
CACHE_SCHEMA_VERSION = 3

DESC_ALL = [
    "desc_QED", "desc_SAscore", "desc_CLogP", "desc_MolWt",
    "desc_HeavyAtomCount", "desc_NumHAcceptors", "desc_NumHDonors",
    "desc_NumHeteroatoms", "desc_NumRotatableBonds",
    "desc_NumAromaticRings", "desc_NumAliphaticRings", "desc_RingCount",
    "desc_TPSA", "desc_PenLogP", "desc_FormalCharge", "desc_Bertz",
    "desc_MaxConsecutiveRotatableBonds", "desc_FlourineCount",
]
DESC_SELECT = [
    "desc_QED", "desc_SAscore", "desc_CLogP", "desc_MolWt",
    "desc_TPSA", "desc_NumHAcceptors", "desc_NumHDonors",
]


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
    cfg.setdefault("umap_mode", "on")
    cfg.setdefault("tsne_mode", "on")
    cfg.setdefault("descriptormode", cfg.get("descriptor_mode", "both"))
    cfg.setdefault("featuremode", cfg.get("feature_mode", "both"))
    cfg.setdefault("cache_mode", "use")
    cfg.setdefault("random_state", "42")
    return cfg


def resolve_path(value, config_path, job=None):
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    if not path.is_absolute():
        path = Path(config_path).resolve().parent / path
    if path.is_dir() and job:
        path = path / f"{job}_scores.csv"
    return path


def canonicalize(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    return Chem.MolToSmiles(mol, canonical=True) if mol else None


def load_source(path, source):
    df = pd.read_csv(path)
    score_col = f"{source}_r_i_docking_score"
    required = ["smiles", "step", score_col]
    missing = [c for c in required if c not in df]
    if missing:
        raise ValueError(f"{source}: missing columns {missing}")
    if "valid" in df:
        valid = df["valid"].astype(str).str.lower().eq("true")
        df = df[valid].copy()
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df[df[score_col].notna() & df[score_col].ne(0)].copy()
    df["source_index"] = df.index
    df["canon_smiles"] = df["smiles"].map(canonicalize)
    df = df.dropna(subset=["canon_smiles"])
    # Retain the best score for duplicates within each source. Cross-source
    # duplicates remain present so PR/PPS overlap is visible.
    best = df.groupby("canon_smiles")[score_col].idxmin()
    df = df.loc[best].copy()
    df["source"] = source
    df["docking_score"] = df[score_col]
    df["label"] = ""
    print(f"[{source}] retained {len(df)} molecules")
    return df


def load_reference(path, reference_name="OmecamtivMecarbil", label="OM"):
    df = pd.read_csv(path)
    required = ["smiles", "name"]
    missing = [c for c in required if c not in df]
    if missing:
        raise ValueError(f"Reference descriptors: missing columns {missing}")

    matches = df[df["name"].astype(str).str.strip() == reference_name]
    if matches.empty:
        raise ValueError(f"Reference '{reference_name}' not found in {path}")

    row = matches.iloc[[0]].copy()
    row["source_index"] = row.index
    row["canon_smiles"] = row["smiles"].map(canonicalize)
    if row["canon_smiles"].isna().any():
        raise ValueError(f"Reference '{reference_name}' has an invalid SMILES")
    row["source"] = "REF"
    row["step"] = np.nan
    row["docking_score"] = np.nan
    row["label"] = label
    print(f"[REF] included {reference_name} as {label}")
    return row


def fingerprint_matrix(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    matrix = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for row, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        fp = generator.GetFingerprint(mol)
        DataStructs.ConvertToNumpyArray(fp, matrix[row])
    return matrix


def cache_signature(paths, method, feature, random_state):
    payload = {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "inputs": [
            {"path": str(p), "size": p.stat().st_size, "mtime": p.stat().st_mtime_ns}
            for p in paths
        ],
        "method": method,
        "feature": feature,
        "random_state": random_state,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def build_one(df, paths, cache_dir, method, feature, random_state, cache_mode):
    output = cache_dir / f"joint_{method}_{feature}.csv.gz"
    metadata = cache_dir / f"joint_{method}_{feature}.json"
    signature = cache_signature(paths, method, feature, random_state)
    if cache_mode == "use" and output.exists() and metadata.exists():
        saved = json.loads(metadata.read_text(encoding="utf-8"))
        if saved.get("signature") == signature:
            print(f"[cache hit] {output}")
            return

    if feature == "fp":
        selected = df.copy()
        population = selected[selected["source"].isin(["PR", "PPS"])].copy()
        references = selected[selected["source"] == "REF"].copy()
        population_matrix = fingerprint_matrix(
            population["canon_smiles"].tolist()
        )
        reference_matrix = fingerprint_matrix(
            references["canon_smiles"].tolist()
        )
    else:
        cols = DESC_ALL if feature == "all_desc" else DESC_SELECT
        missing = [c for c in cols if c not in df]
        if missing:
            raise ValueError(f"{feature}: missing descriptor columns {missing}")
        selected = df.dropna(subset=cols).copy()
        population = selected[selected["source"].isin(["PR", "PPS"])].copy()
        references = selected[selected["source"] == "REF"].copy()
        scaler = StandardScaler()
        population_matrix = scaler.fit_transform(
            population[cols].astype(float)
        )
        reference_matrix = scaler.transform(references[cols].astype(float))

    if references.empty:
        raise ValueError(f"No reference rows available for {method}/{feature}")

    selected = pd.concat([population, references], ignore_index=True)
    if method == "umap":
        import umap
        reducer = umap.UMAP(
            n_neighbors=15,
            min_dist=0.1,
            metric=(
                "jaccard"
                if population_matrix.dtype == np.uint8
                else "euclidean"
            ),
            random_state=random_state,
        )
        print(
            f"[calculate] {method}/{feature}: fit {population_matrix.shape}; "
            f"transform references {reference_matrix.shape}"
        )
        population_embedding = reducer.fit_transform(population_matrix)
        reference_embedding = reducer.transform(reference_matrix)
        embedding = np.vstack([population_embedding, reference_embedding])
        reference_projection = "transform"
    else:
        matrix = np.vstack([population_matrix, reference_matrix])
        print(f"[calculate] {method}/{feature}: joint fit {matrix.shape}")
        reducer = TSNE(
            n_components=2,
            perplexity=30,
            init="pca",
            learning_rate="auto",
            random_state=random_state,
        )
        embedding = reducer.fit_transform(matrix)
        reference_projection = "joint_fit"

    result = selected[
        [
            "source", "source_index", "smiles", "canon_smiles", "step",
            "docking_score", "label",
        ]
    ].copy()
    result["embedding_x"] = embedding[:, 0]
    result["embedding_y"] = embedding[:, 1]
    result["method"] = method
    result["feature"] = feature
    result.to_csv(output, index=False, compression="gzip")
    metadata.write_text(json.dumps({
        "signature": signature, "rows": len(result), "method": method,
        "feature": feature, "random_state": random_state,
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "reference_rows": int(result["source"].eq("REF").sum()),
        "reference_projection": reference_projection,
    }, indent=2), encoding="utf-8")
    print(f"[saved] {output}")


def main(config_path):
    cfg = parse_config(config_path)
    pr_path = resolve_path(cfg["prdir"], config_path, "PR")
    pps_path = resolve_path(cfg["ppsdir"], config_path, "PPS")
    ref_path = resolve_path(
        cfg.get(
            "ref_desc",
            r"Q:\coding_dir\701_Project\workspace\AHC_Related"
            r"\ref_structures\ref_desc.csv",
        ),
        config_path,
    )
    outdir = resolve_path(cfg["outdir"], config_path)
    cache_dir = outdir / "embedding_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    combined = pd.concat([
        load_source(pr_path, "PR"),
        load_source(pps_path, "PPS"),
        load_reference(ref_path),
    ], ignore_index=True)
    methods = [m for m in ("umap", "tsne") if cfg[f"{m}_mode"].lower() == "on"]
    feature_mode = cfg["featuremode"].lower()
    descriptor_mode = cfg["descriptormode"].lower()
    features = []
    if feature_mode in ("fp", "both"):
        features.append("fp")
    if feature_mode in ("desc", "both"):
        if descriptor_mode in ("all", "both"):
            features.append("all_desc")
        if descriptor_mode in ("select", "both"):
            features.append("select_desc")
    for method in methods:
        for feature in features:
            build_one(
                combined, [pr_path, pps_path, ref_path], cache_dir, method, feature,
                int(cfg["random_state"]), cfg["cache_mode"].lower(),
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build joint PR/PPS embeddings.")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)
