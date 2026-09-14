"""Compatibility entry point for relocated exact-pose extraction."""
from Modules.Crossdocking.CompletedAnalysis.candidate_interpretation.extract_top5_state_selective_poses import *  # noqa: F401,F403

if __name__ == "__main__":
    import runpy
    runpy.run_module(
        "Modules.Crossdocking.CompletedAnalysis.candidate_interpretation.extract_top5_state_selective_poses",
        run_name="__main__",
    )
