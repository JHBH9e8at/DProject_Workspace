"""Write own/opposite hotspot tables below work_results."""

import argparse
from pathlib import Path
from Modules.Analysis_Figen_block.common.output_naming import allocate_versioned_group
from Modules.Analysis_Figen_block.common.paths import validate_results_root
from Modules.Analysis_Figen_block.common.provenance import update_manifest
from Modules.Analysis_Figen_block.common.run_context import create_run_context, generate_run_id
from .own_opposite_hotspot import build_own_opposite_tables

SCRIPT_ID = "PY-126"

def run_own_opposite_hotspot(*, own_metadata, own_interactions, opposite_metadata,
                             opposite_interactions, run_id, results_root=None, resume=False,
                             upstream_run_ids=()):
    context = create_run_context(("analysis", "own_opposite_hotspot"), run_id,
                                 results_root=results_root, resume=resume)
    inputs = tuple(Path(p) for p in (own_metadata, own_interactions, opposite_metadata, opposite_interactions))
    try:
        selected, summary, typed, denominators = build_own_opposite_tables(own_metadata=own_metadata,
            own_interactions=own_interactions, opposite_metadata=opposite_metadata,
            opposite_interactions=opposite_interactions)
        paths = allocate_versioned_group(context.tables_dir, ("matched_best_own_poses.csv", "own_opposite_residue_prevalence.csv", "own_interaction_type_prevalence.csv"))
        selected.to_csv(paths[0], index=False); summary.to_csv(paths[1], index=False); typed.to_csv(paths[2], index=False)
        update_manifest(context, script_id=SCRIPT_ID, status="completed",
                        config={"matched_denominators": denominators}, inputs=inputs,
                        outputs=paths, upstream_run_ids=upstream_run_ids)
    except Exception:
        update_manifest(context, script_id=SCRIPT_ID, status="failed", config={},
                        upstream_run_ids=upstream_run_ids); raise
    return context, selected, summary, typed, denominators, paths

def main(argv=None):
    parser=argparse.ArgumentParser(description="Build own/opposite hotspot tables")
    parser.add_argument("--own-metadata",required=True,type=Path); parser.add_argument("--own-interactions",required=True,type=Path)
    parser.add_argument("--opposite-metadata",required=True,type=Path); parser.add_argument("--opposite-interactions",required=True,type=Path)
    parser.add_argument("--run-id"); parser.add_argument("--results-root",type=Path); parser.add_argument("--resume",action="store_true")
    parser.add_argument("--upstream-run-id",action="append",default=[]); args=parser.parse_args(argv)
    root=validate_results_root(args.results_root) if args.results_root else None
    context,*_,paths=run_own_opposite_hotspot(own_metadata=args.own_metadata,own_interactions=args.own_interactions,
        opposite_metadata=args.opposite_metadata,opposite_interactions=args.opposite_interactions,
        run_id=args.run_id or generate_run_id("PR_PPS","own_opposite_hotspot"),results_root=root,
        resume=args.resume,upstream_run_ids=tuple(args.upstream_run_id))
    print(f"Hotspot run: {context.run_dir}"); [print(f"Table: {path}") for path in paths]; return 0

if __name__ == "__main__": raise SystemExit(main())
