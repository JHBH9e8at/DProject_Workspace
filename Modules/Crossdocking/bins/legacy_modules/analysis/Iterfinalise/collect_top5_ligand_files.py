"""Compatibility entry point for relocated candidate-file collection."""
from Modules.Crossdocking.CompletedAnalysis.candidate_interpretation.collect_top5_ligand_files import *  # noqa: F401,F403

if __name__ == "__main__":
    import runpy
    runpy.run_module(
        "Modules.Crossdocking.CompletedAnalysis.candidate_interpretation.collect_top5_ligand_files",
        run_name="__main__",
    )
