"""Versioned full-style own/opposite hotspot Figure workflow."""

import argparse
from pathlib import Path
from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import create_run_context, generate_run_id
from ...interaction_result_analysis.blocks.config_helpers import parse_int_assignments, parse_path_assignments
from .own_opposite_hotspot import render_fullstyle

SCRIPT_ID="PY-129"

def run_figure(*, prevalence_csv, typed_csv, pocket_files, reference_csv, denominators,
               run_id, results_root=None, resume=False, upstream_run_ids=()):
    context=create_run_context(("Figs","own_opposite_hotspot"),run_id,results_root=results_root,resume=resume)
    png,data=allocate_versioned_group(context.run_dir,("own_opposite_hotspot_fullstyle.png","own_opposite_hotspot_fullstyle_data.csv"))
    inputs=[Path(prevalence_csv),Path(typed_csv),Path(reference_csv),*(Path(p) for p in pocket_files.values())]
    try:
        _,generated,table=render_fullstyle(prevalence_csv=prevalence_csv,typed_csv=typed_csv,
            pocket_files=pocket_files,reference_csv=reference_csv,denominators=denominators,output_png=png)
        if generated != data: generated.replace(data)
        update_manifest(context,script_id=SCRIPT_ID,status="completed",config={"denominators":denominators},inputs=inputs,outputs=(png,data),upstream_run_ids=upstream_run_ids)
    except Exception:
        update_manifest(context,script_id=SCRIPT_ID,status="failed",config={},upstream_run_ids=upstream_run_ids); raise
    return context,png,data,table

def main(argv=None):
    p=argparse.ArgumentParser(description="Render own/opposite hotspot figure"); p.add_argument("--prevalence",required=True,type=Path); p.add_argument("--typed",required=True,type=Path); p.add_argument("--pocket-file",action="append",required=True); p.add_argument("--reference",required=True,type=Path); p.add_argument("--denominator",action="append",required=True); p.add_argument("--run-id"); p.add_argument("--results-root",type=Path); p.add_argument("--resume",action="store_true"); p.add_argument("--upstream-run-id",action="append",default=[]); a=p.parse_args(argv)
    c,*_=run_figure(prevalence_csv=a.prevalence,typed_csv=a.typed,pocket_files=parse_path_assignments(a.pocket_file,label="pocket-file"),reference_csv=a.reference,denominators=parse_int_assignments(a.denominator,label="denominator"),run_id=a.run_id or generate_run_id("PR_PPS","own_opposite_hotspot_figure"),results_root=validate_results_root(a.results_root) if a.results_root else None,resume=a.resume,upstream_run_ids=tuple(a.upstream_run_id)); print(f"Figure run: {c.run_dir}"); return 0

if __name__ == "__main__": raise SystemExit(main())
