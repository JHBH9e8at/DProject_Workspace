"""Extract Glide docking scores from an SDF or SDFGZ pose file."""

import argparse
import gzip
from pathlib import Path

import pandas as pd
from rdkit import Chem


def extract_scores(pose_file, output_csv):
    pose_file = Path(pose_file).resolve()
    output_csv = Path(output_csv).resolve()

    if not pose_file.is_file():
        raise FileNotFoundError(f"Pose file not found: {pose_file}")

    if pose_file.suffix.lower() == ".sdfgz":
        handle = gzip.open(pose_file, "rb")
    else:
        handle = pose_file.open("rb")

    rows = []
    with handle:
        supplier = Chem.ForwardSDMolSupplier(handle)
        for pose_number, mol in enumerate(supplier, start=1):
            if mol is None:
                continue
            properties = mol.GetPropsAsDict()
            if "r_i_docking_score" not in properties:
                continue
            rows.append({
                "pose": pose_number,
                "title": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
                "smiles": Chem.MolToSmiles(mol),
                **properties,
            })

    if not rows:
        raise RuntimeError(f"No pose with r_i_docking_score found in {pose_file}")

    scores = pd.DataFrame(rows).sort_values("r_i_docking_score")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(output_csv, index=False)

    best = scores.iloc[0]
    print(f"Best docking score: {best['r_i_docking_score']}")
    print(f"Scores CSV: {output_csv}")
    return scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Glide scores")
    parser.add_argument("--poses", required=True, help="Glide SDF/SDFGZ output")
    parser.add_argument("--output", required=True, help="Output score CSV")
    args = parser.parse_args()

    extract_scores(args.poses, args.output)

