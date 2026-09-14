from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


ARMS = ("L_PPS_to_R_PR", "L_PR_to_R_PPS")


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        yield from csv.DictReader(handle)


def as_float(value):
    if value in (None, "", "nan", "NaN"):
        return None
    return float(value)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None):
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def candidate_key(state: str, molecule_id: str):
    return state, molecule_id


def select_candidates(score_path: Path, top_n: int):
    grouped = defaultdict(list)
    for row in read_csv(score_path):
        row["selectivity_margin"] = float(row["selectivity_margin"])
        row["own_state_score"] = float(row["own_state_score"])
        row["opposite_state_score"] = float(row["opposite_state_score"])
        row["ligand_efficiency_margin"] = as_float(row.get("ligand_efficiency_margin"))
        grouped[row["own_state"]].append(row)
    selected = []
    for state in ("PPS", "PR"):
        eligible = [row for row in grouped[state] if row["selectivity_margin"] > 0]
        eligible.sort(key=lambda row: (-row["selectivity_margin"], row["own_state_score"]))
        for rank, row in enumerate(eligible[:top_n], 1):
            row = dict(row)
            row["candidate_rank"] = rank
            selected.append(row)
    return selected


def load_alerts(pps_path: Path, pr_path: Path):
    output = {}
    for state, path, variant_column in (
        ("PPS", pps_path, "PPS_best_variant"),
        ("PR", pr_path, "PR_best_variant"),
    ):
        for row in read_csv(path):
            molecule_id = row.get(variant_column, "")
            if molecule_id:
                output[candidate_key(state, molecule_id)] = row
    return output


def load_pose_metadata(path: Path, key_state_column: str, selected_keys: set):
    output = {}
    for row in read_csv(path):
        key = candidate_key(row[key_state_column], row["molecule_id"])
        if key in selected_keys:
            output[key] = row
    return output


def load_quality(path: Path, key_state_column: str, selected_keys: set):
    output = {}
    for row in read_csv(path):
        key = candidate_key(row[key_state_column], row["molecule_id"])
        if key in selected_keys:
            output[key] = row
    return output


def load_interactions(path: Path, key_state_column: str, selected_keys: set):
    features = defaultdict(set)
    residues = defaultdict(set)
    residue_types = defaultdict(lambda: defaultdict(set))
    for row in read_csv(path):
        key = candidate_key(row[key_state_column], row["molecule_id"])
        if key not in selected_keys or row.get("interaction_present", "True") != "True":
            continue
        residue = row["original_residue_id"]
        interaction = row["interaction_type"]
        features[key].add((residue, interaction))
        residues[key].add(residue)
        residue_types[key][residue].add(interaction)
    return features, residues, residue_types


def jaccard(left: set, right: set):
    union = left | right
    return len(left & right) / len(union) if union else None


def transition(own: bool, opposite: bool):
    if own and opposite:
        return "retained"
    if own:
        return "lost"
    if opposite:
        return "gained"
    return "absent"


def sorted_join(values):
    return ";".join(sorted(values))


def build_outputs(candidates, alerts, own_meta, opposite_meta, own_quality,
                  opposite_quality, own_features, opposite_features,
                  own_residues, opposite_residues, own_types, opposite_types):
    candidate_rows = []
    feature_rows = []
    residue_rows = []
    for candidate in candidates:
        key = candidate_key(candidate["own_state"], candidate["molecule_id"])
        own_f = own_features.get(key, set())
        opp_f = opposite_features.get(key, set())
        own_r = own_residues.get(key, set())
        opp_r = opposite_residues.get(key, set())
        own_q = own_quality[key]
        opp_q = opposite_quality[key]
        own_m = own_meta[key]
        opp_m = opposite_meta[key]
        alert = alerts.get(key, {})
        candidate_rows.append({
            "own_state": candidate["own_state"],
            "candidate_rank": candidate["candidate_rank"],
            "molecule_id": candidate["molecule_id"],
            "crossdock_arm": candidate["crossdock_arm"],
            "own_state_score": candidate["own_state_score"],
            "opposite_state_score": candidate["opposite_state_score"],
            "selectivity_margin": candidate["selectivity_margin"],
            "ligand_efficiency_margin": candidate["ligand_efficiency_margin"],
            "selectivity_class": candidate["selectivity_class"],
            "own_glide_variant": own_m.get("glide_variant", ""),
            "opposite_glide_variant": opp_m.get("glide_variant", ""),
            "best_crossdock_smiles": candidate.get("best_crossdock_smiles", ""),
            "molecular_weight": alert.get("desc_MolWt", ""),
            "heavy_atom_count": alert.get("desc_HeavyAtomCount", own_m.get("heavy_atom_count", "")),
            "clogp": alert.get("desc_CLogP", ""),
            "tpsa": alert.get("desc_TPSA", ""),
            "rotatable_bonds": alert.get("desc_NumRotatableBonds", ""),
            "PAINS": alert.get("PAINS", ""),
            "BRENK": alert.get("BRENK", ""),
            "own_clash_count": own_q.get("clash_count", ""),
            "opposite_clash_count": opp_q.get("clash_count", ""),
            "own_clashes_per_heavy_atom": own_q.get("clashes_per_heavy_atom", ""),
            "opposite_clashes_per_heavy_atom": opp_q.get("clashes_per_heavy_atom", ""),
            "own_strain_energy": own_q.get("strain_energy", ""),
            "opposite_strain_energy": opp_q.get("strain_energy", ""),
            "own_interaction_feature_count": len(own_f),
            "opposite_interaction_feature_count": len(opp_f),
            "retained_interaction_feature_count": len(own_f & opp_f),
            "lost_interaction_feature_count": len(own_f - opp_f),
            "gained_interaction_feature_count": len(opp_f - own_f),
            "interaction_feature_jaccard": jaccard(own_f, opp_f),
            "own_interacting_residue_count": len(own_r),
            "opposite_interacting_residue_count": len(opp_r),
            "retained_residue_count": len(own_r & opp_r),
            "lost_residue_count": len(own_r - opp_r),
            "gained_residue_count": len(opp_r - own_r),
            "residue_jaccard": jaccard(own_r, opp_r),
            "retained_residues": sorted_join(own_r & opp_r),
            "lost_residues": sorted_join(own_r - opp_r),
            "gained_residues": sorted_join(opp_r - own_r),
        })
        for residue, interaction in sorted(own_f | opp_f):
            feature_rows.append({
                "own_state": candidate["own_state"],
                "candidate_rank": candidate["candidate_rank"],
                "molecule_id": candidate["molecule_id"],
                "residue": residue,
                "interaction_type": interaction,
                "own_present": (residue, interaction) in own_f,
                "opposite_present": (residue, interaction) in opp_f,
                "transition": transition((residue, interaction) in own_f,
                                         (residue, interaction) in opp_f),
            })
        for residue in sorted(own_r | opp_r):
            residue_rows.append({
                "own_state": candidate["own_state"],
                "candidate_rank": candidate["candidate_rank"],
                "molecule_id": candidate["molecule_id"],
                "residue": residue,
                "own_interaction_types": sorted_join(own_types[key].get(residue, set())),
                "opposite_interaction_types": sorted_join(opposite_types[key].get(residue, set())),
                "own_present": residue in own_r,
                "opposite_present": residue in opp_r,
                "transition": transition(residue in own_r, residue in opp_r),
            })
    return candidate_rows, feature_rows, residue_rows


def render_report(candidates: list[dict], output_path: Path):
    lines = [
        "# Top 5 state-selective candidates per directed population",
        "",
        "Candidates are ranked by positive selectivity margin (opposite-state score minus own-state score), with own-state docking score used as the secondary sort key.",
        "",
    ]
    for state in ("PPS", "PR"):
        lines += [f"## {state}-directed candidates", "",
                  "| Rank | Molecule | Own score | Opposite score | Margin | Residue Jaccard | Retained | Lost | Gained |",
                  "|---:|---|---:|---:|---:|---:|---:|---:|---:|"]
        for row in [item for item in candidates if item["own_state"] == state]:
            lines.append(
                f"| {row['candidate_rank']} | {row['molecule_id']} | "
                f"{float(row['own_state_score']):.3f} | {float(row['opposite_state_score']):.3f} | "
                f"{float(row['selectivity_margin']):.3f} | {float(row['residue_jaccard']):.3f} | "
                f"{row['retained_residue_count']} | {row['lost_residue_count']} | {row['gained_residue_count']} |"
            )
        lines.append("")
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main(args):
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    candidates = select_candidates(args.paired_scores, args.top_n)
    selected_keys = {candidate_key(row["own_state"], row["molecule_id"]) for row in candidates}
    alerts = load_alerts(args.pps_filtered, args.pr_filtered)
    own_meta = load_pose_metadata(args.own_pose_metadata, "population_state", selected_keys)
    opposite_meta = load_pose_metadata(args.opposite_pose_metadata, "ligand_state", selected_keys)
    own_quality = load_quality(args.own_pose_quality, "population_state", selected_keys)
    opposite_quality = load_quality(args.opposite_pose_quality, "ligand_state", selected_keys)
    own_f, own_r, own_t = load_interactions(args.own_interactions, "population_state", selected_keys)
    opp_f, opp_r, opp_t = load_interactions(args.opposite_interactions, "ligand_state", selected_keys)

    for label, mapping in (
        ("own metadata", own_meta), ("opposite metadata", opposite_meta),
        ("own quality", own_quality), ("opposite quality", opposite_quality),
    ):
        missing = selected_keys - set(mapping)
        if missing:
            raise ValueError(f"Missing {label} for {sorted(missing)}")

    summary, features, residues = build_outputs(
        candidates, alerts, own_meta, opposite_meta, own_quality, opposite_quality,
        own_f, opp_f, own_r, opp_r, own_t, opp_t,
    )
    summary.sort(key=lambda row: ((0 if row["own_state"] == "PPS" else 1), row["candidate_rank"]))
    write_csv(output_dir / "top5_candidate_summary.csv", summary)
    write_csv(output_dir / "top5_interaction_transitions.csv", features)
    write_csv(output_dir / "top5_residue_transitions.csv", residues)
    render_report(summary, output_dir / "top5_candidate_report.md")
    (output_dir / "manifest.json").write_text(json.dumps({
        "selection": "Top positive selectivity margin per own state; own-state score tie-break",
        "top_n_per_state": args.top_n,
        "candidate_count": len(summary),
        "interaction_identity": "original_residue_id + interaction_type",
        "residue_identity": "original_residue_id",
        "inputs": {key: str(value) for key, value in vars(args).items() if isinstance(value, Path)},
    }, indent=2), encoding="utf-8")
    print(f"Saved {len(summary)} candidates, {len(features)} interaction transitions, and {len(residues)} residue transitions")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--paired-scores", type=Path, required=True)
    parser.add_argument("--pps-filtered", type=Path, required=True)
    parser.add_argument("--pr-filtered", type=Path, required=True)
    parser.add_argument("--own-pose-metadata", type=Path, required=True)
    parser.add_argument("--opposite-pose-metadata", type=Path, required=True)
    parser.add_argument("--own-pose-quality", type=Path, required=True)
    parser.add_argument("--opposite-pose-quality", type=Path, required=True)
    parser.add_argument("--own-interactions", type=Path, required=True)
    parser.add_argument("--opposite-interactions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=5)
    main(parser.parse_args())
