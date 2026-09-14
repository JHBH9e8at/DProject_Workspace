"""Compatibility entry point for relocated candidate interpretation."""
from Modules.Crossdocking.CompletedAnalysis.candidate_interpretation.analyze_top5_state_selective_candidates import *  # noqa: F401,F403

if __name__ == "__main__":
    import runpy
    runpy.run_module(
        "Modules.Crossdocking.CompletedAnalysis.candidate_interpretation.analyze_top5_state_selective_candidates",
        run_name="__main__",
    )
