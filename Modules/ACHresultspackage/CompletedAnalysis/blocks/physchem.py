"""Run PR and PPS physicochemical/trajectory workflows with shared axes."""

from pathlib import Path
import pandas as pd

from Modules.Analysis_Figen_block.figures.blocks.ahc_physchem_trajectory_workflow import run_workflow


DESCRIPTORS = ["desc_MolWt", "desc_NumRotatableBonds", "desc_CLogP"]


def run(*, pr_input: Path, pps_input: Path, reference: Path, results_root: Path,
        run_id: str, last_n: int, resume: bool) -> dict[str, object]:
    maxima = []
    for source, path in (("PR", pr_input), ("PPS", pps_input)):
        score = f"{source}_r_i_docking_score"
        frame = pd.read_csv(path, usecols=[score, "desc_MolWt"])
        maxima.append(float(frame.loc[frame[score] != 0, "desc_MolWt"].max()))
    shared_mw_upper = max(maxima) + 100.0
    outputs = {}
    for source, path, threshold in (("PR", pr_input, -6.0), ("PPS", pps_input, -9.0)):
        outputs[source] = run_workflow(
            input_csv=path, job_name=source, threshold=threshold,
            run_id=f"{run_id}_{source}", results_root=results_root, resume=resume,
            ref_csv=reference, ref_keys=["OM"], descriptors=DESCRIPTORS,
            mw_upper_bound=shared_mw_upper,
        )
    return {"shared_mw_upper": shared_mw_upper, "runs": outputs}
