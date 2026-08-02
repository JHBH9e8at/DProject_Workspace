"""Run PPS/PR interaction analysis from a flat key=value configuration file."""

import argparse


REQUIRED_KEYS = {
    "pps_results",
    "pr_results",
    "pps_run_dir",
    "pr_run_dir",
    "pps_receptor",
    "pr_receptor",
    "output_dir",
}
VALID_ONOFF = {"on", "off"}
VALID_ANALYSIS = {"posecheck", "prolif", "both"}


def _optional_float(cfg, key):
    value = cfg.get(key, "").strip()
    return float(value) if value else None


def parse_config(path):
    """Parse and validate a flat interaction-analysis configuration."""
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
            key = key.strip().lower()
            if key in cfg:
                raise ValueError(f"{path}:{lineno}: duplicate key: {key}")
            cfg[key] = value.strip()

    missing = sorted(key for key in REQUIRED_KEYS if not cfg.get(key))
    if missing:
        raise ValueError(f"{path}: missing required key(s): {missing}")

    defaults = {
        "pps_run_name": "PPS",
        "pr_run_name": "PR",
        "step_col": "step",
        "pps_variant_col": "PPS_best_variant",
        "pr_variant_col": "PR_best_variant",
        "smiles_col": "smiles",
        "analysis": "both",
    }
    for key, value in defaults.items():
        cfg[key] = cfg.get(key, "").strip() or value

    cfg["analysis"] = cfg["analysis"].lower()
    if cfg["analysis"] not in VALID_ANALYSIS:
        raise ValueError(
            f"{path}: analysis must be one of {sorted(VALID_ANALYSIS)}"
        )

    for key in (
        "include_secondary_interactions",
        "allow_unavailable",
        "resume",
        "fail_fast",
    ):
        value = cfg.get(key, "off").strip().lower()
        if value not in VALID_ONOFF:
            raise ValueError(f"{path}: {key} must be 'on' or 'off'")
        cfg[key] = value == "on"

    top_n = cfg.get("testmode_top_n", "").strip()
    cfg["testmode_top_n"] = int(top_n) if top_n else None
    if cfg["testmode_top_n"] is not None and cfg["testmode_top_n"] <= 0:
        raise ValueError(f"{path}: testmode_top_n must be a positive integer")

    cfg["posecheck_workers"] = int(
        cfg.get("posecheck_workers", "").strip() or "1"
    )
    cfg["posecheck_chunk_size"] = int(
        cfg.get("posecheck_chunk_size", "").strip() or "500"
    )
    if cfg["posecheck_workers"] <= 0:
        raise ValueError(f"{path}: posecheck_workers must be a positive integer")
    if cfg["posecheck_chunk_size"] <= 0:
        raise ValueError(f"{path}: posecheck_chunk_size must be a positive integer")
    prolif_workers = cfg.get("prolif_workers", "").strip()
    cfg["prolif_workers"] = int(prolif_workers) if prolif_workers else None
    if cfg["prolif_workers"] is not None and cfg["prolif_workers"] <= 0:
        raise ValueError(f"{path}: prolif_workers must be a positive integer")

    cfg["residue_map"] = cfg.get("residue_map", "").strip() or None
    cfg["max_clashes"] = _optional_float(cfg, "max_clashes")
    cfg["max_clashes_per_heavy_atom"] = _optional_float(
        cfg, "max_clashes_per_heavy_atom"
    )
    cfg["max_strain_energy"] = _optional_float(cfg, "max_strain_energy")
    return cfg


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Interaction .in configuration")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    cfg = parse_config(args.config)
    try:
        from .batch_runner import run_batch_interactions
    except ImportError:
        from batch_runner import run_batch_interactions

    return run_batch_interactions(
        pps_results=cfg["pps_results"],
        pr_results=cfg["pr_results"],
        pps_run_dir=cfg["pps_run_dir"],
        pr_run_dir=cfg["pr_run_dir"],
        pps_receptor=cfg["pps_receptor"],
        pr_receptor=cfg["pr_receptor"],
        output_dir=cfg["output_dir"],
        pps_run_name=cfg["pps_run_name"],
        pr_run_name=cfg["pr_run_name"],
        step_col=cfg["step_col"],
        pps_variant_col=cfg["pps_variant_col"],
        pr_variant_col=cfg["pr_variant_col"],
        smiles_col=cfg["smiles_col"],
        testmode_top_n=cfg["testmode_top_n"],
        analysis=cfg["analysis"],
        residue_map=cfg["residue_map"],
        include_secondary_interactions=cfg["include_secondary_interactions"],
        max_clashes=cfg["max_clashes"],
        max_clashes_per_heavy_atom=cfg["max_clashes_per_heavy_atom"],
        max_strain_energy=cfg["max_strain_energy"],
        posecheck_workers=cfg["posecheck_workers"],
        posecheck_chunk_size=cfg["posecheck_chunk_size"],
        prolif_workers=cfg["prolif_workers"],
        allow_unavailable=cfg["allow_unavailable"],
        resume=cfg["resume"],
        fail_fast=cfg["fail_fast"],
    )


if __name__ == "__main__":
    main()
