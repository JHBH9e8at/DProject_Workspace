"""Compatibility entry point for relocated PoseCheck cohort statistics."""
from Modules.Crossdocking.CompletedAnalysis.score.posecheck_stats import *  # noqa: F401,F403

if __name__ == "__main__":
    from Modules.Crossdocking.CompletedAnalysis.score.posecheck_stats import main
    raise SystemExit(main())
