"""Run consolidated PR/PPS scaffold-diversity analysis."""

from pathlib import Path
from Modules.Analysis_Figen_block.ahc_analysis.blocks.scaffold_diversity_runner import run_scaffold_diversity


def run(*, pr_input: Path, pps_input: Path, membership: Path,
        results_root: Path, run_id: str, last_n: int, resume: bool):
    return run_scaffold_diversity(
        populations={"PR": pr_input, "PPS": pps_input},
        membership_csv=membership, run_id=f"{run_id}_scaffold",
        results_root=results_root, resume=resume, last_n=last_n,
    )
