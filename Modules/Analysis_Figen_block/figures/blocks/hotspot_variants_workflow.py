"""Versioned Figure workflows for PY-015/PY-016 hotspot variants."""

import argparse

from pathlib import Path

from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import create_run_context, generate_run_id
from ...interaction_result_analysis.blocks.config_helpers import parse_int_assignments, parse_path_assignments
from .hotspot_variants import render_full_last, render_top_candidates


SCRIPT_ID = "PY-134"


def _run(*, feature, stem, renderer, renderer_kwargs, inputs, config,
         run_id, results_root, resume, upstream_run_ids=()):
    context = create_run_context(("Figs", feature), run_id, results_root=results_root, resume=resume)
    png, data = allocate_versioned_group(
        context.run_dir, (f"{stem}.png", f"{stem}_data.csv")
    )
    try:
        _, generated_data, table = renderer(output_png=png, **renderer_kwargs)
        if generated_data != data:
            generated_data.replace(data)
        update_manifest(
            context, script_id=SCRIPT_ID, status="completed", config=config,
            inputs=tuple(Path(path) for path in inputs), outputs=(png, data),
            upstream_run_ids=upstream_run_ids,
        )
    except Exception:
        update_manifest(context, script_id=SCRIPT_ID, status="failed", config=config, upstream_run_ids=upstream_run_ids)
        raise
    return context, png, data, table


def run_full_last_figure(*, full_csv, last_csv, typed_csv, pocket_files,
                         reference_csv, denominators, run_id,
                         results_root=None, resume=False, upstream_run_ids=()):
    return _run(
        feature="residue_hotspot_full_last", stem="residue_hotspot_full_last",
        renderer=render_full_last,
        renderer_kwargs=dict(full_csv=full_csv, last_csv=last_csv, typed_csv=typed_csv,
                             pocket_files=pocket_files, reference_csv=reference_csv,
                             denominators=denominators),
        inputs=(full_csv, last_csv, typed_csv, reference_csv, *pocket_files.values()),
        config={"denominators": denominators}, run_id=run_id,
        results_root=results_root, resume=resume, upstream_run_ids=upstream_run_ids,
    )


def run_top_candidate_figure(*, own_csv, opposite_csv, typed_csv, pocket_files,
                             reference_csv, denominators, run_id,
                             results_root=None, resume=False, upstream_run_ids=()):
    return _run(
        feature="top_candidate_hotspot", stem="top_candidate_hotspot",
        renderer=render_top_candidates,
        renderer_kwargs=dict(own_csv=own_csv, opposite_csv=opposite_csv,
                             typed_csv=typed_csv, pocket_files=pocket_files,
                             reference_csv=reference_csv, denominators=denominators),
        inputs=(own_csv, opposite_csv, typed_csv, reference_csv, *pocket_files.values()),
        config={"denominators": denominators}, run_id=run_id,
        results_root=results_root, resume=resume, upstream_run_ids=upstream_run_ids,
    )


def _common_parser(description):
    parser=argparse.ArgumentParser(description=description)
    parser.add_argument("--typed",required=True,type=Path); parser.add_argument("--pocket-file",action="append",required=True)
    parser.add_argument("--reference",required=True,type=Path); parser.add_argument("--denominator",action="append",required=True)
    parser.add_argument("--run-id"); parser.add_argument("--results-root",type=Path); parser.add_argument("--resume",action="store_true")
    parser.add_argument("--upstream-run-id",action="append",default=[]); return parser


def main_full_last(argv=None):
    p=_common_parser("Render full/last hotspot figure"); p.add_argument("--full",required=True,type=Path); p.add_argument("--last",required=True,type=Path); a=p.parse_args(argv)
    c,*_=run_full_last_figure(full_csv=a.full,last_csv=a.last,typed_csv=a.typed,pocket_files=parse_path_assignments(a.pocket_file,label="pocket-file"),reference_csv=a.reference,denominators=parse_int_assignments(a.denominator,label="denominator"),run_id=a.run_id or generate_run_id("PR_PPS","hotspot_full_last_figure"),results_root=validate_results_root(a.results_root) if a.results_root else None,resume=a.resume,upstream_run_ids=tuple(a.upstream_run_id)); print(f"Figure run: {c.run_dir}"); return 0


def main_top_candidates(argv=None):
    p=_common_parser("Render top-candidate hotspot figure"); p.add_argument("--own",required=True,type=Path); p.add_argument("--opposite",required=True,type=Path); a=p.parse_args(argv)
    c,*_=run_top_candidate_figure(own_csv=a.own,opposite_csv=a.opposite,typed_csv=a.typed,pocket_files=parse_path_assignments(a.pocket_file,label="pocket-file"),reference_csv=a.reference,denominators=parse_int_assignments(a.denominator,label="denominator"),run_id=a.run_id or generate_run_id("PR_PPS","top_candidate_hotspot_figure"),results_root=validate_results_root(a.results_root) if a.results_root else None,resume=a.resume,upstream_run_ids=tuple(a.upstream_run_id)); print(f"Figure run: {c.run_dir}"); return 0
