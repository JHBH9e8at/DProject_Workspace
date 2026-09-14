"""Compatibility entry point for the relocated selectivity figure."""
from Modules.Crossdocking.CompletedAnalysis.figures.plot_crossdocking_selectivity import *  # noqa: F401,F403

if __name__ == "__main__":
    import runpy
    runpy.run_module(
        "Modules.Crossdocking.CompletedAnalysis.figures.plot_crossdocking_selectivity",
        run_name="__main__",
    )
