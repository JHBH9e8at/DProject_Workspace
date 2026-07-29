"""Tests for the config-driven PPS/PR interaction batch runner."""

import tempfile
import unittest
import gzip
from pathlib import Path

import pandas as pd

from analysis.interaction.batch_runner import (
    build_population_manifest,
    run_batch_interactions,
)
from analysis.interaction.interaction_runner import parse_config


def _scores(state, rows):
    score_col = f"{state}_r_i_docking_score"
    variant_col = f"{state}_best_variant"
    return pd.DataFrame(
        [
            {
                "step": index,
                "smiles": f"C{'C' * index}",
                "valid": True,
                "unique": True,
                score_col: score,
                variant_col: f"{index}_0-1",
            }
            for index, score in enumerate(rows, start=1)
        ]
    )


def _write_pose_files(run_dir, state, count):
    for index in range(1, count + 1):
        pose = (
            run_dir
            / f"{state}_GlideDock"
            / str(index)
            / f"{index}_0-1_lib.sdfgz"
        )
        pose.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(pose, "wb") as handle:
            handle.write(
                f"pose-{index}\n  test\n\nM  END\n$$$$\n".encode("utf-8")
            )


def _fake_analyzer(
    protein_file,
    pose_file,
    receptor_state,
    output_dir,
    analysis,
    **_,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [{"receptor_state": receptor_state, "pose_index": 0, "docking_score": -8.0}]
    ).to_csv(output_dir / "pose_metadata.csv", index=False)
    if analysis in {"posecheck", "both"}:
        pd.DataFrame(
            [{"receptor_state": receptor_state, "pose_index": 0, "clash_count": 0}]
        ).to_csv(output_dir / "pose_quality.csv", index=False)
    if analysis in {"prolif", "both"}:
        pd.DataFrame(
            [{
                "receptor_state": receptor_state,
                "pose_index": 0,
                "interaction_type": "HBDonor",
            }]
        ).to_csv(output_dir / "interactions_long.csv", index=False)
    pd.DataFrame(
        [{"receptor_state": receptor_state, "analysis_component": analysis}]
    ).to_csv(output_dir / "run_summary.csv", index=False)


class InteractionBatchTests(unittest.TestCase):
    def test_config_defaults_and_blank_test_mode(self):
        with tempfile.TemporaryDirectory() as temp:
            config = Path(temp) / "interaction.in"
            config.write_text(
                "\n".join(
                    [
                        "PPS_results=pps.csv",
                        "PR_results=pr.csv",
                        "pps_run_dir=pps_run",
                        "pr_run_dir=pr_run",
                        "pps_receptor=pps.pdb",
                        "pr_receptor=pr.pdb",
                        "output_dir=out",
                        "testmode_top_n=",
                    ]
                ),
                encoding="utf-8",
            )
            parsed = parse_config(config)
            self.assertEqual(parsed["pps_run_name"], "PPS")
            self.assertEqual(parsed["pr_run_name"], "PR")
            self.assertEqual(parsed["pps_variant_col"], "PPS_best_variant")
            self.assertEqual(parsed["testmode_top_n"], None)

    def test_manifest_test_mode_selects_best_scores(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            scores = root / "pps.csv"
            run_dir = root / "pps_run"
            receptor = root / "pps.pdb"
            _scores("PPS", [-8.0, -10.0, -9.0]).to_csv(scores, index=False)
            _write_pose_files(run_dir, "PPS", 3)
            receptor.write_text("ATOM\n", encoding="utf-8")

            manifest = build_population_manifest(
                results_csv=scores,
                run_dir=run_dir,
                run_name="PPS",
                receptor_file=receptor,
                state="PPS",
                testmode_top_n=2,
            )
            self.assertEqual(manifest["molecule_id"].tolist(), ["2_0-1", "3_0-1"])
            self.assertTrue(manifest["pose_file_exists"].all())

    def test_batch_writes_two_state_aggregates(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            pps_scores = root / "pps.csv"
            pr_scores = root / "pr.csv"
            pps_run = root / "pps_run"
            pr_run = root / "pr_run"
            pps_receptor = root / "pps.pdb"
            pr_receptor = root / "pr.pdb"
            output = root / "output"
            _scores("PPS", [-8.0, -10.0]).to_csv(pps_scores, index=False)
            _scores("PR", [-7.0, -9.0]).to_csv(pr_scores, index=False)
            _write_pose_files(pps_run, "PPS", 2)
            _write_pose_files(pr_run, "PR", 2)
            pps_receptor.write_text("ATOM\n", encoding="utf-8")
            pr_receptor.write_text("ATOM\n", encoding="utf-8")

            manifest, summary = run_batch_interactions(
                pps_results=pps_scores,
                pr_results=pr_scores,
                pps_run_dir=pps_run,
                pr_run_dir=pr_run,
                pps_receptor=pps_receptor,
                pr_receptor=pr_receptor,
                output_dir=output,
                testmode_top_n=1,
                analysis="both",
                analyzer=_fake_analyzer,
            )
            self.assertEqual(len(manifest), 2)
            self.assertEqual(summary["selection_mode"], "testmode_top_n")
            self.assertEqual(summary["status_counts"], {"completed": 2})
            interactions = pd.read_csv(
                output / "aggregate" / "interactions_long.csv"
            )
            self.assertEqual(set(interactions["population_state"]), {"PPS", "PR"})
            self.assertEqual(len(interactions), 2)


if __name__ == "__main__":
    unittest.main()
