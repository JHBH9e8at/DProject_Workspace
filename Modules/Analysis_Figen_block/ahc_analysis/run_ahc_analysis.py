"""User-facing configuration entry point for completed AHC result analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKING_ROOT = HERE.parents[2]
if str(WORKING_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKING_ROOT))

from Modules.ACHresultspackage.CompletedAnalysis.run_analysis import (  # noqa: E402
    parse_config,
    run,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run all enabled completed AHC analysis blocks"
    )
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args(argv)
    run(args.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
