"""Run the combined PR/PPS last-window comparison."""

from pathlib import Path
from Modules.Analysis_Figen_block.figures.blocks.physchem_last_window_workflow import run_last_window


def run(*, pr_input: Path, pps_input: Path, results_root: Path,
        run_id: str, last_n: int, resume: bool):
    return run_last_window(
        pr_csv=pr_input, pps_csv=pps_input, last_n=last_n,
        pr_threshold=-6.0, pps_threshold=-9.0, run_id=f"{run_id}_last_window",
        results_root=results_root, resume=resume,
    )
