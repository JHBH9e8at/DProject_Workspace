"""Generate versioned chemical-space figures from existing embedding products."""

from __future__ import annotations

import argparse
from pathlib import Path

from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import create_run_context, generate_run_id
from .chemical_space import load_embedding_cache, render_docking, render_iteration


SCRIPT_ID = "PY-113"


def run_chemical_space_figures(*, cache_path: str | Path, run_id: str,
                               results_root: str | Path | None = None, resume: bool = False,
                               upstream_run_ids=()):
    frame = load_embedding_cache(cache_path)
    method = str(frame.method.iloc[0]).lower(); feature = str(frame.feature.iloc[0])
    context = create_run_context(("Figs", "chemical_space", f"{method}_{feature}"), run_id,
                                 results_root=results_root, resume=resume)
    names = (f"joint_{method}_{feature}.png", f"joint_{method}_{feature}_docking.png")
    outputs = allocate_versioned_group(context.run_dir, names)
    config = {"method": method, "feature": feature}
    try:
        render_iteration(cache_path, outputs[0]); render_docking(cache_path, outputs[1])
        update_manifest(context, script_id=SCRIPT_ID, status="completed", config=config,
                        inputs=(Path(cache_path),), outputs=outputs, upstream_run_ids=upstream_run_ids)
    except Exception:
        update_manifest(context, script_id=SCRIPT_ID, status="failed", config=config, upstream_run_ids=upstream_run_ids); raise
    return context, outputs

def main(argv=None):
    p=argparse.ArgumentParser(description="Render chemical-space figures"); p.add_argument("--cache",required=True,type=Path); p.add_argument("--run-id"); p.add_argument("--results-root",type=Path); p.add_argument("--resume",action="store_true"); p.add_argument("--upstream-run-id",action="append",default=[]); a=p.parse_args(argv)
    c,_=run_chemical_space_figures(cache_path=a.cache,run_id=a.run_id or generate_run_id("PR_PPS","chemical_space_figures"),results_root=validate_results_root(a.results_root) if a.results_root else None,resume=a.resume,upstream_run_ids=tuple(a.upstream_run_id)); print(f"Figure run: {c.run_dir}"); return 0

if __name__ == "__main__": raise SystemExit(main())
