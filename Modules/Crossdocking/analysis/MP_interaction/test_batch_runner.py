"""Tests for the config-driven PPS/PR interaction batch runner."""

import tempfile
import unittest
import gzip
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from analysis.MP_interaction.batch_runner import (
    build_population_manifest,
    run_batch_interactions,
)
from analysis.MP_interaction.interaction_runner import parse_config
from analysis.MP_interaction.pose_quality import (
    run_posecheck_parallel,
    split_sdf_chunks,
)


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


def _fake_posecheck_worker(protein_pdb, chunk_spec, result_file):
    start = int(chunk_spec["start_index"])
    count = int(chunk_spec["pose_count"])
    table = pd.DataFrame(
        {
            "chunk_index": int(chunk_spec["chunk_index"]),
            "pose_index": range(start, start + count),
            "clash_count": range(start, start + count),
            "strain_energy": [float(index) + 0.5 for index in range(start, start + count)],
            "chunk_sha256": chunk_spec["chunk_sha256"],
            "worker_pid": 1,
            "chunk_elapsed_seconds": 0.01,
        }
    )
    table.to_csv(result_file, index=False)
    return {
        "chunk_index": int(chunk_spec["chunk_index"]),
        "pose_count": count,
        "result_file": str(result_file),
        "elapsed_seconds": 0.01,
        "worker_pid": 1,
    }


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
            self.assertEqual(parsed["posecheck_workers"], 1)
            self.assertEqual(parsed["posecheck_chunk_size"], 500)
            self.assertIsNone(parsed["prolif_workers"])

    def test_config_parallel_worker_settings(self):
        with tempfile.TemporaryDirectory() as temp:
            config = Path(temp) / "interaction.in"
            config.write_text(
                "\n".join(
                    [
                        "pps_results=pps.csv",
                        "pr_results=pr.csv",
                        "pps_run_dir=pps_run",
                        "pr_run_dir=pr_run",
                        "pps_receptor=pps.pdb",
                        "pr_receptor=pr.pdb",
                        "output_dir=out",
                        "posecheck_workers=16",
                        "posecheck_chunk_size=250",
                        "prolif_workers=16",
                    ]
                ),
                encoding="utf-8",
            )
            parsed = parse_config(config)
            self.assertEqual(parsed["posecheck_workers"], 16)
            self.assertEqual(parsed["posecheck_chunk_size"], 250)
            self.assertEqual(parsed["prolif_workers"], 16)

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

    def test_parallel_posecheck_chunks_preserve_global_pose_indices(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            receptor = root / "protein.pdb"
            ligand = root / "poses.sdf"
            receptor.write_text("ATOM\n", encoding="utf-8")
            ligand.write_text(
                "".join(
                    f"pose-{index}\n  test\n\nM  END\n$$$$\n"
                    for index in range(5)
                ),
                encoding="utf-8",
            )
            metadata = pd.DataFrame(
                {
                    "pose_index": range(5),
                    "pose_read_status": ["valid"] * 5,
                    "heavy_atom_count": [10] * 5,
                    "rotatable_bond_count": [2] * 5,
                }
            )
            quality = run_posecheck_parallel(
                receptor,
                ligand,
                metadata,
                root / "chunks",
                workers=2,
                chunk_size=2,
                executor_factory=lambda workers: ThreadPoolExecutor(
                    max_workers=workers
                ),
                worker_function=_fake_posecheck_worker,
            )
            self.assertEqual(quality["pose_index"].tolist(), list(range(5)))
            self.assertEqual(quality["clash_count"].tolist(), list(range(5)))
            self.assertTrue(quality["posecheck_status"].eq("completed").all())
            self.assertEqual(
                len(list((root / "chunks" / "sdf_chunks").glob("chunk_*.sdf"))),
                3,
            )

    def test_parallel_posecheck_resume_reuses_valid_chunks(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            receptor = root / "protein.pdb"
            ligand = root / "poses.sdf"
            receptor.write_text("ATOM\n", encoding="utf-8")
            ligand.write_text(
                "".join(
                    f"pose-{index}\n  test\n\nM  END\n$$$$\n"
                    for index in range(4)
                ),
                encoding="utf-8",
            )
            metadata = pd.DataFrame(
                {
                    "pose_index": range(4),
                    "pose_read_status": ["valid"] * 4,
                }
            )
            kwargs = dict(
                protein_pdb=receptor,
                ligand_sdf=ligand,
                pose_metadata=metadata,
                work_dir=root / "chunks",
                workers=2,
                chunk_size=2,
                executor_factory=lambda workers: ThreadPoolExecutor(
                    max_workers=workers
                ),
                worker_function=_fake_posecheck_worker,
            )
            first = run_posecheck_parallel(**kwargs)

            def fail_if_called(*_args, **_kwargs):
                raise AssertionError("valid checkpoint should have been reused")

            kwargs["worker_function"] = fail_if_called
            kwargs["resume"] = True
            second = run_posecheck_parallel(**kwargs)
            pd.testing.assert_frame_equal(first, second)


if __name__ == "__main__":
    unittest.main()
