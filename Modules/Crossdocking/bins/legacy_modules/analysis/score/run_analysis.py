"""Compatibility import for relocated selectivity analysis CLI."""
from Modules.Crossdocking.CompletedAnalysis.score.run_analysis import *  # noqa: F401,F403

if __name__ == "__main__":
    import runpy
    runpy.run_module("Modules.Crossdocking.CompletedAnalysis.score.run_analysis", run_name="__main__")
