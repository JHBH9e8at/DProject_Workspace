"""Coarse-grained entry point for the complete Crossdocking calculation block."""
from __future__ import annotations
import argparse
from pathlib import Path
from Modules.Analysis_Figen_block.common.paths import assert_within, get_results_root, validate_results_root
from .crossdocking_runner import main as calculation_main, parse_config

def main(argv=None):
    parser=argparse.ArgumentParser(description="Run the complete Crossdocking calculation block")
    parser.add_argument("--config",required=True,type=Path)
    parser.add_argument("--results-root",type=Path,help="Allowed work_results boundary")
    args=parser.parse_args(argv)
    boundary=validate_results_root(args.results_root) if args.results_root else get_results_root(create=True)
    cfg=parse_config(args.config)
    assert_within(cfg["output_dir"],boundary)
    calculation_main(["--config",str(args.config)])
    return 0

if __name__=="__main__": raise SystemExit(main())
