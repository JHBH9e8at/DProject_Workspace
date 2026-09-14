"""Analysis/Figure workflow for the PY-009 last-N comparison."""

import argparse
from pathlib import Path
from ...ahc_analysis.blocks.physchem_last_window import prepare_pr_pps
from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import create_run_context, generate_run_id
from .physchem_last_window import render_last_window


SCRIPT_ID = "PY-149"


def run_last_window(*, pr_csv, pps_csv, last_n=50, pr_threshold=None,
                    pps_threshold=None, run_id, results_root=None, resume=False,
                    upstream_run_ids=()):
    analysis = create_run_context(("analysis", "physchem_last_window"), run_id,
                                  results_root=results_root, resume=resume)
    figures = create_run_context(("Figs", "physchem_last_window"), run_id,
                                 results_root=results_root, resume=resume)
    config = {"last_n": last_n, "pr_threshold": pr_threshold,
              "pps_threshold": pps_threshold}
    rows, summary = prepare_pr_pps(pr_csv, pps_csv, last_n=last_n)
    rows_path, summary_path = allocate_versioned_group(
        analysis.tables_dir, ("physchem_last_window_rows.csv", "physchem_last_window_summary.csv"))
    rows.to_csv(rows_path, index=False); summary.to_csv(summary_path, index=False)
    update_manifest(analysis, script_id=SCRIPT_ID, status="completed", config=config,
                    inputs=(Path(pr_csv), Path(pps_csv)), outputs=(rows_path, summary_path), upstream_run_ids=upstream_run_ids)
    (png,) = allocate_versioned_group(figures.run_dir,
                                      (f"PR_PPS_physchem_docking_scatter_last{last_n}.png",))
    render_last_window(rows_csv=rows_path, summary_csv=summary_path, output_png=png,
                       last_n=last_n, pr_threshold=pr_threshold,
                       pps_threshold=pps_threshold)
    update_manifest(figures, script_id="PY-148", status="completed", config=config,
                    inputs=(rows_path, summary_path), outputs=(png,), upstream_run_ids=(*upstream_run_ids,analysis.run_id))
    return analysis, figures, (rows_path, summary_path), png

def main(argv=None):
    p=argparse.ArgumentParser(description="Generate PR/PPS last-window physicochemical comparison"); p.add_argument("--pr-input",required=True,type=Path); p.add_argument("--pps-input",required=True,type=Path); p.add_argument("--last-n",type=int,default=50); p.add_argument("--pr-threshold",type=float); p.add_argument("--pps-threshold",type=float); p.add_argument("--run-id"); p.add_argument("--results-root",type=Path); p.add_argument("--resume",action="store_true"); p.add_argument("--upstream-run-id",action="append",default=[]); a=p.parse_args(argv)
    x,*_=run_last_window(pr_csv=a.pr_input,pps_csv=a.pps_input,last_n=a.last_n,pr_threshold=a.pr_threshold,pps_threshold=a.pps_threshold,run_id=a.run_id or generate_run_id("PR_PPS","physchem_last_window"),results_root=validate_results_root(a.results_root) if a.results_root else None,resume=a.resume,upstream_run_ids=tuple(a.upstream_run_id)); print(f"Analysis run: {x.run_dir}"); return 0

if __name__ == "__main__": raise SystemExit(main())
