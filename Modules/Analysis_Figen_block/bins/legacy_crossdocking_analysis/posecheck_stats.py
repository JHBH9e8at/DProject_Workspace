"""Compatibility wrapper for the canonical Crossdocking implementation."""

from Modules.Crossdocking.CompletedAnalysis.score.posecheck_stats import *  # noqa: F401,F403
from Modules.Crossdocking.CompletedAnalysis.score.posecheck_stats import main


if __name__ == "__main__":
    raise SystemExit(main())
