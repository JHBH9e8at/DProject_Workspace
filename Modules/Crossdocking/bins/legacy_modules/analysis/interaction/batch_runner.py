"""Compatibility boundary for shared serial batch orchestration."""

from ..interaction_common.batch_core import (
    STATE_SPECS,
    _collect_population_outputs,
    _job_is_complete,
    _materialize_population_sdf,
    _safe_name,
    _true_mask,
    _utc_now,
    build_population_manifest,
    run_batch_interactions,
    select_population,
)

__all__ = [
    "STATE_SPECS", "select_population", "build_population_manifest",
    "run_batch_interactions", "_collect_population_outputs", "_job_is_complete",
    "_materialize_population_sdf", "_safe_name", "_true_mask", "_utc_now",
]
