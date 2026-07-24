"""Run the complete production cross-docking and score-analysis workflow."""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


BASE_REQUIRED_KEYS = [
    "pps_input",
    "pps_run_dir",
    "pr_input",
    "pr_run_dir",
    "pps_grid",
    "pr_grid",
    "output_dir",
]
VALID_ONOFF = {"on", "off"}
VALID_PRECISION = {"HTVS", "SP", "XP"}
VALID_LIGAND_PREP_MODE = {"on", "off"}


def _utc_now():
    return datetime.now(timezone.utc).isoformat()


def _write_run_summary(summary, summary_file):
    summary_file.parent.mkdir(parents=True, exist_ok=True)
    summary_file.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def parse_config(path):
    """Parse and validate a flat key=value cross-docking config."""
    cfg = {}
    with open(path, "r", encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(
                    f"{path}:{lineno}: invalid line (expected key=value): {raw!r}"
                )
            key, _, value = line.partition("=")
            key = key.strip()
            if key in cfg:
                raise ValueError(f"{path}:{lineno}: duplicate key: {key}")
            cfg[key] = value.strip()

    missing = [key for key in BASE_REQUIRED_KEYS if not cfg.get(key)]
    if missing:
        raise ValueError(f"{path}: missing required key(s): {missing}")

    precision = cfg.get("precision", "SP").strip().upper()
    if precision not in VALID_PRECISION:
        raise ValueError(
            f"{path}: invalid precision {precision!r}, "
            f"must be from {sorted(VALID_PRECISION)}"
        )
    cfg["precision"] = precision

    ligand_prep_mode = cfg.get("ligand_prep_mode", "off").strip().lower()
    if ligand_prep_mode not in VALID_LIGAND_PREP_MODE:
        raise ValueError(
            f"{path}: invalid ligand_prep_mode {ligand_prep_mode!r}, "
            "must be 'on' or 'off'"
        )
    cfg["ligand_prep_mode"] = ligand_prep_mode

    for key in ("allow_unavailable", "resume", "fail_fast"):
        value = cfg.get(key, "off").strip().lower()
        if value not in VALID_ONOFF:
            raise ValueError(
                f"{path}: invalid {key} value {value!r}, must be 'on' or 'off'"
            )
        cfg[key] = value == "on"

    cfg["best_dscore_top_per"] = float(
        cfg.get("best_dscore_top_per", "1.0").strip()
    )
    if not 0 < cfg["best_dscore_top_per"] <= 100:
        raise ValueError(
            f"{path}: best_dscore_top_per must be greater than 0 and at most 100"
        )
    testmode_top_n = cfg.get("testmode_top_n", "").strip()
    cfg["testmode_top_n"] = int(testmode_top_n) if testmode_top_n else None
    if cfg["testmode_top_n"] is not None and cfg["testmode_top_n"] <= 0:
        raise ValueError(f"{path}: testmode_top_n must be a positive integer")
    cfg["selectivity_threshold"] = float(cfg.get("selectivity_threshold", "2.0"))
    cfg["strong_score_threshold"] = float(
        cfg.get("strong_score_threshold", "-8.0")
    )

    for key in (
        "pps_variant_col",
        "pr_variant_col",
        "pps_scores",
        "pr_scores",
    ):
        cfg[key] = cfg.get(key, "").strip() or None

    if bool(cfg["pps_scores"]) != bool(cfg["pr_scores"]):
        raise ValueError(f"{path}: pps_scores and pr_scores must be provided together")

    return cfg


def run_full_pipeline(
    pps_input,
    pps_run_dir,
    pr_input,
    pr_run_dir,
    pps_grid,
    pr_grid,
    output_dir,
    pps_run_name="PPS",
    pr_run_name="PR",
    step_col="step",
    pps_variant_col=None,
    pr_variant_col=None,
    smiles_col="smiles",
    best_dscore_top_per=1.0,
    testmode_top_n=None,
    precision="SP",
    allow_unavailable=False,
    resume=False,
    fail_fast=False,
    ligand_prep_mode="off",
    selectivity_threshold=2.0,
    strong_score_threshold=-8.0,
    pps_scores=None,
    pr_scores=None,
    pps_score_col="PPS_r_i_docking_score",
    pr_score_col="PR_r_i_docking_score",
    valid_col="valid",
    config_file=None,
):
    """Run both cross-docking arms followed by the complete score block."""
    from analysis.score.score_runner import run_score_block
    from pipelines.double_arm import run_double_arm
    from prep.filter_population import filter_population

    if bool(pps_scores) != bool(pr_scores):
        raise ValueError("--pps-scores and --pr-scores must be provided together")

    output_dir = Path(output_dir).resolve()
    filtering_dir = output_dir / "filtering"
    docking_dir = output_dir / "crossdocking"
    analysis_dir = output_dir / "analysis"
    summary_file = output_dir / "pipeline_run_summary.json"
    output_dir.mkdir(parents=True, exist_ok=True)
    pps_filtered = filtering_dir / "PPS_filtered.csv"
    pr_filtered = filtering_dir / "PR_filtered.csv"

    summary = {
        "status": "running",
        "started_at_utc": _utc_now(),
        "finished_at_utc": None,
        "current_block": "filtering",
        "pps_raw_input": str(pps_input),
        "pr_raw_input": str(pr_input),
        "pps_filtered_input": str(pps_filtered),
        "pr_filtered_input": str(pr_filtered),
        "best_dscore_top_per": best_dscore_top_per,
        "testmode_top_n": testmode_top_n,
        "selection_mode": (
            "testmode_top_n"
            if testmode_top_n is not None
            else "best_dscore_top_per"
        ),
        "crossdocking_output": str(docking_dir),
        "analysis_output": str(analysis_dir),
        "interaction_analysis": "not_run_single_complex_block",
        "config_file": str(config_file) if config_file else None,
        "error": None,
    }
    _write_run_summary(summary, summary_file)

    try:
        print("\n" + "=" * 72)
        print("FILTERING PPS POPULATION")
        print("=" * 72)
        filter_population(
            input_csv=pps_input,
            output_csv=pps_filtered,
            run_name=pps_run_name,
            best_dscore_top_per=best_dscore_top_per,
            testmode_top_n=testmode_top_n,
        )

        print("\n" + "=" * 72)
        print("FILTERING PR POPULATION")
        print("=" * 72)
        filter_population(
            input_csv=pr_input,
            output_csv=pr_filtered,
            run_name=pr_run_name,
            best_dscore_top_per=best_dscore_top_per,
            testmode_top_n=testmode_top_n,
        )

        summary["current_block"] = "crossdocking"
        _write_run_summary(summary, summary_file)
        run_double_arm(
            pps_input=pps_filtered,
            pps_run_dir=pps_run_dir,
            pr_input=pr_filtered,
            pr_run_dir=pr_run_dir,
            pps_grid=pps_grid,
            pr_grid=pr_grid,
            output_dir=docking_dir,
            pps_run_name=pps_run_name,
            pr_run_name=pr_run_name,
            step_col=step_col,
            pps_variant_col=pps_variant_col,
            pr_variant_col=pr_variant_col,
            # Top-N selection was already applied after filtering and sorting.
            limit=None,
            precision=precision,
            allow_unavailable=allow_unavailable,
            resume=resume,
            fail_fast=fail_fast,
            ligand_prep_mode=ligand_prep_mode,
            smiles_col=smiles_col,
        )

        summary["current_block"] = "score_analysis"
        _write_run_summary(summary, summary_file)
        run_score_block(
            input_dir=docking_dir,
            output_dir=analysis_dir,
            selectivity_threshold=selectivity_threshold,
            strong_score_threshold=strong_score_threshold,
            pps_scores=pps_scores,
            pr_scores=pr_scores,
            pps_score_col=pps_score_col,
            pr_score_col=pr_score_col,
            smiles_col=smiles_col,
            valid_col=valid_col,
        )
    except Exception as error:
        summary["status"] = "failed"
        summary["finished_at_utc"] = _utc_now()
        summary["error"] = f"{type(error).__name__}: {error}"
        _write_run_summary(summary, summary_file)
        raise

    summary["status"] = "completed"
    summary["current_block"] = None
    summary["finished_at_utc"] = _utc_now()
    _write_run_summary(summary, summary_file)

    print("\nFull cross-docking workflow completed")
    print(f"Cross-docking output: {docking_dir}")
    print(f"Score-analysis output: {analysis_dir}")
    print(f"Run summary: {summary_file}")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Full cross-docking runner\n"
            "1. PPS -> PR cross-docking\n"
            "2. PR -> PPS cross-docking\n"
            "3. score QC and selectivity analysis\n"
            "4. optional raw-population overlap\n\n"
            "Config format: flat key=value lines (see t_crossdock.in)"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to cross-docking config (.in) file",
    )
    args = parser.parse_args(argv)
    cfg = parse_config(args.config)

    return run_full_pipeline(
        pps_input=cfg["pps_input"],
        pps_run_dir=cfg["pps_run_dir"],
        pr_input=cfg["pr_input"],
        pr_run_dir=cfg["pr_run_dir"],
        pps_grid=cfg["pps_grid"],
        pr_grid=cfg["pr_grid"],
        output_dir=cfg["output_dir"],
        pps_run_name=cfg.get("pps_run_name", "PPS"),
        pr_run_name=cfg.get("pr_run_name", "PR"),
        step_col=cfg.get("step_col", "step"),
        pps_variant_col=cfg["pps_variant_col"],
        pr_variant_col=cfg["pr_variant_col"],
        smiles_col=cfg.get("smiles_col", "smiles"),
        best_dscore_top_per=cfg["best_dscore_top_per"],
        testmode_top_n=cfg["testmode_top_n"],
        precision=cfg["precision"],
        allow_unavailable=cfg["allow_unavailable"],
        resume=cfg["resume"],
        fail_fast=cfg["fail_fast"],
        ligand_prep_mode=cfg["ligand_prep_mode"],
        selectivity_threshold=cfg["selectivity_threshold"],
        strong_score_threshold=cfg["strong_score_threshold"],
        pps_scores=cfg["pps_scores"],
        pr_scores=cfg["pr_scores"],
        pps_score_col=cfg.get("pps_score_col", "PPS_r_i_docking_score"),
        pr_score_col=cfg.get("pr_score_col", "PR_r_i_docking_score"),
        valid_col=cfg.get("valid_col", "valid"),
        config_file=Path(args.config).resolve(),
    )


if __name__ == "__main__":
    main()
