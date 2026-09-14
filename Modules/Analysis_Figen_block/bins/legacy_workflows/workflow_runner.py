"""Stable package entry point for registered end-to-end workflows."""

from __future__ import annotations

import argparse
from collections.abc import Callable

from ..figures.blocks.ahc_physchem_trajectory_workflow import main as ahc_physchem_main
from ..structural_analysis.blocks.redocking_workflow import main as redocking_main
from ..figures.blocks.chemical_space_workflow import main as chemical_space_figures_main
from ..figures.blocks.hotspot_variants_workflow import main_full_last as hotspot_full_last_figure_main, main_top_candidates as hotspot_top_candidates_figure_main
from ..figures.blocks.own_opposite_hotspot_workflow import main as own_opposite_hotspot_figure_main
from ..figures.blocks.physchem_last_window_workflow import main as physchem_last_window_main
from ..ahc_analysis.blocks.chemical_space_runner import main as chemical_space_main
from ..ahc_analysis.blocks.population_runner import main as ahc_population_main
from ..ahc_analysis.blocks.scaffold_diversity_runner import main as scaffold_diversity_main
from ..structural_analysis.blocks.pocket_volume_runner import main as pocket_volume_main
from ..crossdocking_analysis.score_workflow import main as crossdocking_score_main
from ..crossdocking_analysis.posecheck_stats import main as posecheck_stats_main
from ..interaction_result_analysis.blocks.hotspot_variants_runner import (
    main_full_last as hotspot_full_last_main,
    main_top_candidates as hotspot_top_candidates_main,
)
from ..interaction_result_analysis.blocks.own_opposite_hotspot_runner import main as own_opposite_hotspot_main
from ..interaction_result_analysis.blocks.residue_hotspot_runner import main as residue_hotspot_main


WORKFLOWS: dict[str, Callable[[list[str] | None], int]] = {
    "ahc-physchem-trajectory": ahc_physchem_main,
    "ahc-population": ahc_population_main,
    "scaffold-diversity": scaffold_diversity_main,
    "chemical-space": chemical_space_main,
    "chemical-space-figures": chemical_space_figures_main,
    "physchem-last-window": physchem_last_window_main,
    "pocket-volume": pocket_volume_main,
    "own-opposite-hotspot-figure": own_opposite_hotspot_figure_main,
    "hotspot-full-last-figure": hotspot_full_last_figure_main,
    "hotspot-top-candidates-figure": hotspot_top_candidates_figure_main,
    "crossdocking-score-analysis": crossdocking_score_main,
    "posecheck-stats": posecheck_stats_main,
    "residue-hotspot": residue_hotspot_main,
    "own-opposite-hotspot": own_opposite_hotspot_main,
    "hotspot-full-last": hotspot_full_last_main,
    "hotspot-top-candidates": hotspot_top_candidates_main,
    "redocking-validation": redocking_main,
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a registered Analysis_Figen_block workflow."
    )
    parser.add_argument("workflow", choices=sorted(WORKFLOWS))
    parser.add_argument(
        "workflow_args",
        nargs=argparse.REMAINDER,
        help="Arguments passed unchanged to the selected workflow.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return WORKFLOWS[args.workflow](args.workflow_args)


if __name__ == "__main__":
    raise SystemExit(main())
