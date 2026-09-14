"""Sequential interaction-analysis adapter."""

from .common.batch_core import (
    _collect_population_outputs, _job_is_complete, _materialize_population_sdf,
    build_population_manifest, run_batch_interactions, select_population,
)

__all__ = [
    "select_population", "build_population_manifest", "run_batch_interactions",
    "_collect_population_outputs", "_job_is_complete", "_materialize_population_sdf",
]
