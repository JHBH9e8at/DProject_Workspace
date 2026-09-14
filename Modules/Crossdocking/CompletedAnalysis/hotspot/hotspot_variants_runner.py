"""Version stable hotspot-variant tables under work_results."""
import argparse
from pathlib import Path
from Modules.Analysis_Figen_block.common.output_naming import allocate_versioned_group
from Modules.Analysis_Figen_block.common.paths import validate_results_root
from Modules.Analysis_Figen_block.common.provenance import update_manifest
from Modules.Analysis_Figen_block.common.run_context import create_run_context, generate_run_id
from .config import parse_int_assignments
from .hotspot_variants import summarize_full_last,load_top_candidates
SCRIPT_ID="PY-131"

def _write(feature,run_id,results_root,resume,tables,names,inputs,config,upstream_run_ids=()):
    ctx=create_run_context(("analysis",feature),run_id,results_root=results_root,resume=resume)
    paths=allocate_versioned_group(ctx.tables_dir,names)
    for table,path in zip(tables,paths): table.to_csv(path,index=False)
    update_manifest(ctx,script_id=SCRIPT_ID,status="completed",config=config,inputs=tuple(Path(p) for p in inputs),outputs=paths,upstream_run_ids=upstream_run_ids)
    return ctx,paths

def run_full_last(*,interactions_csv,metadata_csv,total_poses,last_start,run_id,results_root=None,resume=False,upstream_run_ids=()):
    full,last,denoms=summarize_full_last(interactions_csv,metadata_csv,total_poses,last_start)
    ctx,paths=_write("residue_hotspot_full_last",run_id,results_root,resume,(full,last),
        ("full_prevalence.csv","last_window_prevalence.csv"),(interactions_csv,metadata_csv),{"total_poses":total_poses,"last_start":last_start,"last_denominators":denoms},upstream_run_ids)
    return ctx,full,last,denoms,paths

def run_top_candidates(*,summary_csv,transitions_csv,run_id,results_root=None,resume=False,upstream_run_ids=()):
    own,opposite,typed,denoms=load_top_candidates(summary_csv,transitions_csv)
    ctx,paths=_write("top_candidate_hotspot",run_id,results_root,resume,(own,opposite,typed),
        ("own_prevalence.csv","opposite_prevalence.csv","own_type_prevalence.csv"),(summary_csv,transitions_csv),{"denominators":denoms},upstream_run_ids)
    return ctx,own,opposite,typed,denoms,paths

def main_full_last(argv=None):
    parser=argparse.ArgumentParser(description="Build full/last-window hotspot tables")
    parser.add_argument("--interactions",required=True,type=Path); parser.add_argument("--metadata",required=True,type=Path)
    parser.add_argument("--total-poses",action="append",required=True); parser.add_argument("--last-start",action="append",required=True)
    parser.add_argument("--run-id"); parser.add_argument("--results-root",type=Path); parser.add_argument("--resume",action="store_true")
    parser.add_argument("--upstream-run-id",action="append",default=[]); args=parser.parse_args(argv)
    root=validate_results_root(args.results_root) if args.results_root else None
    context,*_,paths=run_full_last(interactions_csv=args.interactions,metadata_csv=args.metadata,
        total_poses=parse_int_assignments(args.total_poses,label="total-poses"),
        last_start=parse_int_assignments(args.last_start,label="last-start"),
        run_id=args.run_id or generate_run_id("PR_PPS","hotspot_full_last"),results_root=root,
        resume=args.resume,upstream_run_ids=tuple(args.upstream_run_id))
    print(f"Hotspot run: {context.run_dir}"); [print(f"Table: {path}") for path in paths]; return 0

def main_top_candidates(argv=None):
    parser=argparse.ArgumentParser(description="Build top-candidate hotspot tables")
    parser.add_argument("--summary",required=True,type=Path); parser.add_argument("--transitions",required=True,type=Path)
    parser.add_argument("--run-id"); parser.add_argument("--results-root",type=Path); parser.add_argument("--resume",action="store_true")
    parser.add_argument("--upstream-run-id",action="append",default=[]); args=parser.parse_args(argv)
    root=validate_results_root(args.results_root) if args.results_root else None
    context,*_,paths=run_top_candidates(summary_csv=args.summary,transitions_csv=args.transitions,
        run_id=args.run_id or generate_run_id("PR_PPS","top_candidate_hotspot"),results_root=root,
        resume=args.resume,upstream_run_ids=tuple(args.upstream_run_id))
    print(f"Hotspot run: {context.run_dir}"); [print(f"Table: {path}") for path in paths]; return 0
