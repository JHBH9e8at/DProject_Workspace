"""Run joint chemical-space calculation below work_results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ...common.output_naming import allocate_versioned_group
from ...common.paths import validate_results_root
from ...common.provenance import update_manifest
from ...common.run_context import create_run_context, generate_run_id
from .chemical_space import prepare_joint, calculate_embedding


SCRIPT_ID = "PY-137"


def run_chemical_space(*, pr_path, pps_path, reference_path, methods=("umap", "tsne"),
                       features=("fp", "all_desc", "select_desc"), random_state=42,
                       run_id, results_root=None, resume=False, upstream_run_ids=()):
    context = create_run_context(("analysis", "chemical_space"), run_id,
                                 results_root=results_root, resume=resume)
    joint = prepare_joint(pr_path, pps_path, reference_path)
    names = tuple(name for method in methods for feature in features
                  for name in (f"joint_{method}_{feature}.csv.gz", f"joint_{method}_{feature}.json"))
    paths = allocate_versioned_group(context.tables_dir, names)
    outputs = []; summaries = {}
    try:
        iterator = iter(paths)
        for method in methods:
            for feature in features:
                csv_path, json_path = next(iterator), next(iterator)
                table, metadata = calculate_embedding(joint, method=method, feature=feature,
                                                      random_state=random_state)
                table.to_csv(csv_path, index=False, compression="gzip")
                json_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
                outputs.extend((csv_path, json_path)); summaries[f"{method}/{feature}"] = metadata
        update_manifest(context, script_id=SCRIPT_ID, status="completed",
                        config={"methods": list(methods), "features": list(features),
                                "random_state": random_state, "summaries": summaries},
                        inputs=(pr_path, pps_path, reference_path), outputs=outputs, upstream_run_ids=upstream_run_ids)
    except Exception:
        update_manifest(context, script_id=SCRIPT_ID, status="failed",
                        config={"methods": list(methods), "features": list(features),
                                "random_state": random_state}, upstream_run_ids=upstream_run_ids)
        raise
    return context, tuple(outputs)

def main(argv=None):
    p=argparse.ArgumentParser(description="Calculate joint AHC chemical-space embeddings")
    p.add_argument("--pr-input",required=True,type=Path); p.add_argument("--pps-input",required=True,type=Path); p.add_argument("--reference",required=True,type=Path)
    p.add_argument("--method",action="append",choices=("umap","tsne")); p.add_argument("--feature",action="append",choices=("fp","all_desc","select_desc")); p.add_argument("--random-state",type=int,default=42)
    p.add_argument("--run-id"); p.add_argument("--results-root",type=Path); p.add_argument("--resume",action="store_true"); p.add_argument("--upstream-run-id",action="append",default=[]); a=p.parse_args(argv)
    c,_=run_chemical_space(pr_path=a.pr_input,pps_path=a.pps_input,reference_path=a.reference,methods=tuple(a.method or ("umap","tsne")),features=tuple(a.feature or ("fp","all_desc","select_desc")),random_state=a.random_state,run_id=a.run_id or generate_run_id("PR_PPS","chemical_space"),results_root=validate_results_root(a.results_root) if a.results_root else None,resume=a.resume,upstream_run_ids=tuple(a.upstream_run_id)); print(f"Chemical-space run: {c.run_dir}"); return 0

if __name__ == "__main__": raise SystemExit(main())
