"""Stable joint PR/PPS chemical-space embedding calculations."""

from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler


RDLogger.DisableLog("rdApp.*")
SCRIPT_ID = "PY-136"
CACHE_SCHEMA_VERSION = 3

DESC_ALL = [
    "desc_QED", "desc_SAscore", "desc_CLogP", "desc_MolWt",
    "desc_HeavyAtomCount", "desc_NumHAcceptors", "desc_NumHDonors",
    "desc_NumHeteroatoms", "desc_NumRotatableBonds", "desc_NumAromaticRings",
    "desc_NumAliphaticRings", "desc_RingCount", "desc_TPSA", "desc_PenLogP",
    "desc_FormalCharge", "desc_Bertz", "desc_MaxConsecutiveRotatableBonds",
    "desc_FlourineCount",
]
DESC_SELECT = ["desc_QED", "desc_SAscore", "desc_CLogP", "desc_MolWt",
               "desc_TPSA", "desc_NumHAcceptors", "desc_NumHDonors"]


def canonicalize(smiles):
    molecule = Chem.MolFromSmiles(str(smiles))
    return Chem.MolToSmiles(molecule, canonical=True) if molecule else None


def load_source(path, source):
    frame = pd.read_csv(path)
    score_column = f"{source}_r_i_docking_score"
    missing = [column for column in ("smiles", "step", score_column) if column not in frame]
    if missing:
        raise ValueError(f"{source}: missing columns {missing}")
    if "valid" in frame:
        frame = frame[frame.valid.astype(str).str.lower().eq("true")].copy()
    frame[score_column] = pd.to_numeric(frame[score_column], errors="coerce")
    frame = frame[frame[score_column].notna() & frame[score_column].ne(0)].copy()
    frame["source_index"] = frame.index
    frame["canon_smiles"] = frame.smiles.map(canonicalize)
    frame = frame.dropna(subset=["canon_smiles"])
    frame = frame.loc[frame.groupby("canon_smiles")[score_column].idxmin()].copy()
    frame["source"] = source
    frame["docking_score"] = frame[score_column]
    frame["label"] = ""
    return frame


def load_reference(path, reference_name="OmecamtivMecarbil", label="OM"):
    frame = pd.read_csv(path)
    missing = [column for column in ("smiles", "name") if column not in frame]
    if missing:
        raise ValueError(f"Reference descriptors: missing columns {missing}")
    matches = frame[frame.name.astype(str).str.strip().eq(reference_name)]
    if matches.empty:
        raise ValueError(f"Reference '{reference_name}' not found in {path}")
    row = matches.iloc[[0]].copy()
    row["source_index"] = row.index
    row["canon_smiles"] = row.smiles.map(canonicalize)
    if row.canon_smiles.isna().any():
        raise ValueError(f"Reference '{reference_name}' has an invalid SMILES")
    row["source"] = "REF"; row["step"] = np.nan
    row["docking_score"] = np.nan; row["label"] = label
    return row


def prepare_joint(pr_path, pps_path, reference_path, **reference_options):
    return pd.concat((load_source(pr_path, "PR"), load_source(pps_path, "PPS"),
                      load_reference(reference_path, **reference_options)), ignore_index=True)


def fingerprint_matrix(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    matrix = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for row, value in enumerate(smiles):
        DataStructs.ConvertToNumpyArray(
            generator.GetFingerprint(Chem.MolFromSmiles(value)), matrix[row]
        )
    return matrix


def calculate_embedding(frame, *, method, feature, random_state=42):
    if feature == "fp":
        selected = frame.copy()
        population = selected[selected.source.isin(("PR", "PPS"))].copy()
        references = selected[selected.source.eq("REF")].copy()
        population_matrix = fingerprint_matrix(population.canon_smiles.tolist())
        reference_matrix = fingerprint_matrix(references.canon_smiles.tolist())
    else:
        columns = DESC_ALL if feature == "all_desc" else DESC_SELECT
        missing = [column for column in columns if column not in frame]
        if missing:
            raise ValueError(f"{feature}: missing descriptor columns {missing}")
        selected = frame.dropna(subset=columns).copy()
        population = selected[selected.source.isin(("PR", "PPS"))].copy()
        references = selected[selected.source.eq("REF")].copy()
        scaler = StandardScaler()
        population_matrix = scaler.fit_transform(population[columns].astype(float))
        reference_matrix = scaler.transform(references[columns].astype(float))
    if references.empty:
        raise ValueError(f"No reference rows available for {method}/{feature}")
    selected = pd.concat((population, references), ignore_index=True)
    if method == "umap":
        import umap
        reducer = umap.UMAP(n_neighbors=15, min_dist=.1,
                            metric="jaccard" if population_matrix.dtype == np.uint8 else "euclidean",
                            random_state=random_state)
        embedding = np.vstack((reducer.fit_transform(population_matrix),
                               reducer.transform(reference_matrix)))
        projection = "transform"
    elif method == "tsne":
        matrix = np.vstack((population_matrix, reference_matrix))
        embedding = TSNE(n_components=2, perplexity=30, init="pca",
                         learning_rate="auto", random_state=random_state).fit_transform(matrix)
        projection = "joint_fit"
    else:
        raise ValueError(f"Unknown embedding method: {method}")
    result = selected[["source", "source_index", "smiles", "canon_smiles", "step",
                       "docking_score", "label"]].copy()
    result["embedding_x"] = embedding[:, 0]; result["embedding_y"] = embedding[:, 1]
    result["method"] = method; result["feature"] = feature
    metadata = {"rows": len(result), "method": method, "feature": feature,
                "random_state": random_state, "cache_schema_version": CACHE_SCHEMA_VERSION,
                "reference_rows": int(result.source.eq("REF").sum()),
                "reference_projection": projection}
    return result, metadata
