"""Single top-level entry point for the repository's executable blocks."""
from __future__ import annotations
import argparse
from .Analysis_Figen_block.run_analysis import main as analysis_main
from .Crossdocking.run_crossdocking import main as crossdocking_main

BLOCKS={"analysis":analysis_main,"crossdocking":crossdocking_main}

def main(argv=None):
    parser=argparse.ArgumentParser(description="Run a complete repository module block")
    parser.add_argument("block",choices=sorted(BLOCKS)); parser.add_argument("block_args",nargs=argparse.REMAINDER)
    args=parser.parse_args(argv)
    block_args=args.block_args[1:] if args.block_args[:1]==["--"] else args.block_args
    return BLOCKS[args.block](block_args)

if __name__=="__main__": raise SystemExit(main())
