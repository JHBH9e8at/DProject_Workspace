"""Single configuration-driven entry point for in-run AHC monitoring."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
from types import ModuleType


HERE = Path(__file__).resolve().parent
BLOCKS_DIR = HERE / "blocks"
DEFAULT_REF_CSV = HERE.parents[2] / "AHC_Related" / "ref_structures" / "ref_desc.csv"

BASE_REQUIRED_KEYS = ("analysismode", "odir", "jn", "threshold")
VALID_ANALYSISMODE = {"mp", "fa"}
VALID_REF = {"OM", "MAV", "AFI"}
VALID_ONOFF = {"on", "off"}
VALID_MODE = {"all", "select", "both"}
VALID_FEAT = {"fp", "desc", "both"}
BLOCK_FILES = {
    "merge": "0_0_maketemp.py", "clean": "0_0_dclean.py",
    "physchem": "0_1_Physioprop.py", "trajectory": "0_2_dcurve.py",
    "umap": "1_1_UMAP.py", "tsne": "1_2_tsne.py",
}


def _load_block(key: str) -> ModuleType:
    path = BLOCKS_DIR / BLOCK_FILES[key]
    spec = importlib.util.spec_from_file_location(f"ahc_indirun_{key}", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load monitoring block: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _resolve_path(value: str, config_dir: Path) -> Path:
    path = Path(value).expanduser()
    return (config_dir / path).resolve(strict=False) if not path.is_absolute() else path.resolve(strict=False)


def parse_config(path: str | Path) -> dict[str, object]:
    """Read and validate a monitoring key=value file."""
    config_path = Path(path).expanduser().resolve(strict=False)
    cfg: dict[str, object] = {}
    with config_path.open("r", encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(f"{config_path}:{lineno}: expected key=value")
            key, value = line.split("=", 1)
            key = key.strip().lower()
            if key in cfg:
                raise ValueError(f"{config_path}:{lineno}: duplicate key {key!r}")
            cfg[key] = value.strip()

    missing = [key for key in BASE_REQUIRED_KEYS if not cfg.get(key)]
    if missing:
        raise ValueError(f"{config_path}: missing required keys: {missing}")
    analysis_mode = str(cfg["analysismode"]).lower()
    if analysis_mode not in VALID_ANALYSISMODE:
        raise ValueError(f"analysismode must be one of {sorted(VALID_ANALYSISMODE)}")
    cfg["analysismode"] = analysis_mode
    selected_input = "basedir" if analysis_mode == "mp" else "indir"
    if not cfg.get(selected_input):
        raise ValueError(f"analysismode={analysis_mode} requires {selected_input!r}")

    job_name = str(cfg["jn"]).upper()
    if job_name not in {"PR", "PPS"}:
        raise ValueError("jn must be PR or PPS")
    cfg["jn"] = job_name
    cfg["threshold"] = float(str(cfg["threshold"]))
    refs = [item.strip().upper() for item in str(cfg.get("ref", "OM")).split(",") if item.strip()]
    bad_refs = [item for item in refs if item not in VALID_REF]
    if bad_refs:
        raise ValueError(f"ref contains unsupported values: {bad_refs}")
    cfg["ref"] = refs

    for key in ("umap", "tsne"):
        value = str(cfg.get(key, "off")).lower()
        if value not in VALID_ONOFF:
            raise ValueError(f"{key} must be on or off")
        cfg[key] = value == "on"
    mode, feat = str(cfg.get("mode", "both")).lower(), str(cfg.get("feat", "both")).lower()
    if mode not in VALID_MODE:
        raise ValueError(f"mode must be one of {sorted(VALID_MODE)}")
    if feat not in VALID_FEAT:
        raise ValueError(f"feat must be one of {sorted(VALID_FEAT)}")
    cfg["mode"], cfg["feat"] = mode, feat

    config_dir = config_path.parent
    for key in ("odir", selected_input):
        cfg[key] = _resolve_path(str(cfg[key]), config_dir)
    cfg["ref_csv"] = _resolve_path(str(cfg.get("ref_csv", DEFAULT_REF_CSV)), config_dir)
    cfg["config_path"] = config_path
    return cfg


def run(config_path: str | Path) -> dict[str, Path]:
    """Execute the complete monitoring pipeline described by one config file."""
    cfg = parse_config(config_path)
    output_dir, job_name = Path(cfg["odir"]), str(cfg["jn"])
    figure_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    if cfg["analysismode"] == "mp":
        _, input_csv = _load_block("merge").gettempres(str(cfg["basedir"]), str(output_dir), job_name)
        if input_csv is None:
            raise RuntimeError("The merge block did not produce an input CSV")
        input_path = Path(input_csv)
    else:
        input_path = Path(cfg["indir"])

    _load_block("clean").opcleaning(str(input_path), None, job_name)
    cleaned_path = input_path.with_name(f"{input_path.stem}{job_name}_cleaned.csv")
    if not cleaned_path.is_file():
        raise RuntimeError(f"The cleaning block did not produce {cleaned_path}")

    _load_block("physchem").mainr(
        str(cleaned_path), str(figure_dir), job_name, float(cfg["threshold"]),
        list(cfg["ref"]), str(cfg["ref_csv"]),
    )
    _load_block("trajectory").gen_learn(str(cleaned_path), str(figure_dir), job_name)
    if cfg["umap"]:
        _load_block("umap").run_umap(
            str(cleaned_path), str(figure_dir), job_name, str(cfg["ref_csv"]),
            str(cfg["mode"]), list(cfg["ref"]), str(cfg["feat"]),
        )
    if cfg["tsne"]:
        _load_block("tsne").run_tsne(
            str(cleaned_path), str(figure_dir), job_name, str(cfg["ref_csv"]),
            str(cfg["mode"]), list(cfg["ref"]), str(cfg["feat"]),
        )
    return {"input_csv": input_path, "cleaned_csv": cleaned_path, "figures": figure_dir}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the complete in-run AHC monitoring block")
    parser.add_argument("--config", required=True, type=Path, help="Monitoring .in file")
    args = parser.parse_args(argv)
    run(args.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
