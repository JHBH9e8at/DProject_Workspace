import argparse
import importlib.util
import os


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

base_dir = os.path.dirname(os.path.abspath(__file__))

BASE_REQUIRED_KEYS = ["analysismode", "odir", "jn", "threshold"]
VALID_ANALYSISMODE = {"mp", "fa"}
VALID_REF = {"OM", "MAV", "AFI"}
VALID_ONOFF = {"on", "off"}
VALID_MODE = {"all", "select", "both"}
VALID_FEAT = {"fp", "desc", "both"}

REF_IND_DEFAULT = "/home/andy/proj/701/Analysis/ref_desc.csv"  # hardcoded for now is for ref thingys if we change targtet this line needs to change!



def parse_config(path):
    cfg = {}
    with open(path, "r") as f:
        for lineno, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(f"{path}:{lineno}: invalid line (expected key=value): {raw!r}")
            key, _, val = line.partition("=")
            cfg[key.strip()] = val.strip()
 
    missing = [k for k in BASE_REQUIRED_KEYS if k not in cfg or cfg[k] == ""]
    if missing:
        raise ValueError(f"{path}: missing required key(s): {missing}")
 
    analysismode = cfg["analysismode"].strip().lower()
    if analysismode not in VALID_ANALYSISMODE:
        raise ValueError(f"{path}: invalid analysismode value {analysismode!r}, must be from {sorted(VALID_ANALYSISMODE)}")
    cfg["analysismode"] = analysismode
 
    if analysismode == "mp":
        if not cfg.get("basedir"):
            raise ValueError(f"{path}: analysismode=mp requires 'basedir' (path to iterations dir)")
    else:  # fa
        if not cfg.get("indir"):
            raise ValueError(f"{path}: analysismode=fa requires 'indir' (path to final scores.csv)")
 
    cfg["threshold"] = float(cfg["threshold"])
 
    ref_raw = cfg.get("ref", "")
    ref_list = [r.strip() for r in ref_raw.split(",") if r.strip()]
    bad_ref = [r for r in ref_list if r not in VALID_REF]
    if bad_ref:
        raise ValueError(f"{path}: invalid ref value(s) {bad_ref}, must be from {sorted(VALID_REF)}")
    cfg["ref"] = ref_list
 
    for key in ("umap", "tsne"):
        val = cfg.get(key, "off").strip().lower()
        if val not in VALID_ONOFF:
            raise ValueError(f"{path}: invalid {key} value {val!r}, must be 'on' or 'off'")
        cfg[key] = (val == "on")
 
    mode = cfg.get("mode", "both").strip().lower()
    if mode not in VALID_MODE:
        raise ValueError(f"{path}: invalid mode value {mode!r}, must be from {sorted(VALID_MODE)}")
    cfg["mode"] = mode
 
    feat = cfg.get("feat", "both").strip().lower()
    if feat not in VALID_FEAT:
        raise ValueError(f"{path}: invalid feat value {feat!r}, must be from {sorted(VALID_FEAT)}")
    cfg["feat"] = feat
 
    return cfg


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
"""initial analysis runner..,
Always check Job name (jn) in your config file !!!!!
0. generate temp   : 0_0_maketemp.py
1. Dataset cleaing : 0_0_dclean.py
2. physchem distri : 0_1_Physioprop.py
3. Learning curve  : 0_2_dcurve.py
4. umap            : 1_1_UMAP.py
4. tsne            : 1_2_tsne.py

config file(config.in) format::
                required arg list [ cannot be ommited ]
    analysis mode       : mp[midpoint] / fa[finished]
    if fa - indir       : input [[csv location]]
    if mp - basedir     : iteration(csv) directory

    odir                : output location
    jn                  : job name!!!!
    threshold           :-9.0 or -6
                optional arg list [ can be ommited ]
    ref=OM,MAV,AFI      : deafult   = none, chose from OM,MAV,AFI
    umap=on             : deafult   = OFF, [ON/OFF]
    tsne=on             : deafult   = OFF, [ON/OFF]
    mode=both   descriptor selection : default = Both [all/select/both]
        selects which descriptors to use for generation.
    feat=both   Feature selector     : default = Both [fp/desc/both]
        Selects which features to use to generate umap or tsne [finger print or descriptors]
""",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--config", type=str, required=True, help="Path to config (.in/.cfg) file")
    args = parser.parse_args()
    cfg = parse_config(args.config)

    proxycsv = load_module(os.path.join(base_dir, "0_0_maketemp.py"), "proxycsv")
    dclean = load_module(os.path.join(base_dir, "0_0_dclean.py"), "dclean")
    physio = load_module(os.path.join(base_dir, "0_1_Physioprop.py"), "physio")
    dcurve = load_module(os.path.join(base_dir, "0_2_dcurve.py"), "dcurve")
    umap_mod = load_module(os.path.join(base_dir, "1_1_UMAP.py"), "umap_mod")
    tsne_mod = load_module(os.path.join(base_dir, "1_2_tsne.py"), "tsne_mod")


    odir = cfg["odir"]
    jn = cfg["jn"]
    threshold = cfg["threshold"]
    ref = cfg["ref"]
    figdir = os.path.join(odir,"figures")

 
    if not os.path.isdir(odir):
        os.makedirs(odir, exist_ok=True)
    if not os.path.isdir(figdir):
        os.makedirs(figdir, exist_ok=True)


            ###### this block is only for mp analysis !!!
    if cfg["analysismode"] == "mp":
        _, indir = proxycsv.gettempres(cfg["basedir"], odir, jn)
        if indir is None:
            raise RuntimeError("gettempres failed to produce a merged csv -- check basedir contents")
    else:
        indir = cfg["indir"]
 
    


    df_clean = dclean.opcleaning(indir, None, jn)
    cleaned_path = indir.replace(".csv", f"{jn}_cleaned.csv")
    physio.mainr(cleaned_path, figdir, jn, threshold, ref)
    dcurve.gen_learn(cleaned_path, figdir, jn)
 
    # optional route 
    if cfg["umap"]:
        umap_mod.run_umap(cleaned_path, figdir, jn, REF_IND_DEFAULT, cfg["mode"], ref, cfg["feat"])
    if cfg["tsne"]:
        tsne_mod.run_tsne(cleaned_path, figdir, jn, REF_IND_DEFAULT, cfg["mode"], ref, cfg["feat"])


# python 0_2_runner.py --config /home/andy/proj/701/code/developing_dir/tempdevtest/run.in
# python runanaly.py --config /home/andy/proj/701/code/developing_dir/t_mp.in
# see what file name.. many patches.


# example of run.in

    # analysismode=fa
    # indir=/home/andy/proj/701/code/developing_dir/tempdevtest/scores.csv
    # odir=/home/andy/proj/701/code/developing_dir/tempdevtest
    # jn=PPS
    # threshold=-9.0
    # ref=OM,MAV,AFI
    # umap=on
    # tsne=off
    # mode=select
    # feat=fp

    # analysismode=mp
    # basedir=/home/andy/proj/701/output/2026_06_02_SMILES-RNN_PPS_tester_run_2/iterations
    # odir=/home/andy/proj/701/code/developing_dir/tempdevtest
    # jn=PPS
    # threshold=-9.0
    # ref=OM,MAV,AFI
    # umap=on
    # tsne=off
    # mode=select
    # feat=fp