"""Compatibility import for relocated population-overlap analysis."""
from Modules.Crossdocking.CompletedAnalysis.score.population_overlap import *  # noqa: F401,F403

if __name__ == "__main__":
    import runpy
    runpy.run_module("Modules.Crossdocking.CompletedAnalysis.score.population_overlap", run_name="__main__")
