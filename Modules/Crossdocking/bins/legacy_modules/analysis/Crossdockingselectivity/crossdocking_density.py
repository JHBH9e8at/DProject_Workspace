"""Compatibility entry point for the relocated cross-docking density figure."""
from Modules.Crossdocking.CompletedAnalysis.figures.crossdocking_density import *  # noqa: F401,F403

if __name__ == "__main__":
    import runpy
    runpy.run_module(
        "Modules.Crossdocking.CompletedAnalysis.figures.crossdocking_density",
        run_name="__main__",
    )
