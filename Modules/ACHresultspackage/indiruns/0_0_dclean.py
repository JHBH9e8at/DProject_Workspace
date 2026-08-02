import os
import argparse
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

def convertcano(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None, None

        canon = Chem.MolToSmiles(mol, canonical=True)
        canon_mol = Chem.MolFromSmiles(canon)
        return canon, canon_mol

    except:
        return None, None

def make_mask(df, mol_col="_mol"):
    pains_params = FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains_catalog = FilterCatalog(pains_params)
    brenk_params = FilterCatalogParams()
    brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    brenk_catalog = FilterCatalog(brenk_params)

    def check_pains(mol):
        if mol is None:
            return None
        return pains_catalog.HasMatch(mol)
    def check_brenk(mol):
        if mol is None:
            return None
        return brenk_catalog.HasMatch(mol)

    df["PAINS"] = df[mol_col].apply(check_pains)
    df["BRENK"] = df[mol_col].apply(check_brenk)
    return df

def opcleaning(indir: str, sd: str, rn: str):
    if rn is None:
        rn = "PPS"

    steps = []  # (step_name, count_before, count_after)
    def track(step_name, before, after):
        removed = before - after
        steps.append((step_name, before, after, removed))
        print(f"{step_name}: {before} -> {after} (removed {removed})")

    print("reprot")
    df = pd.read_csv(indir)
    n0 = len(df)
    print(f"ini:: {n0}")

    dfval = df[df["valid"] == True].copy()
    track("valid filter", n0, len(dfval))

    dfvalu = dfval[dfval["unique"] == True].copy()
    track("unique filter", len(dfval), len(dfvalu))

    dfcano = dfvalu.copy()
    cano_results = dfcano["smiles"].apply(convertcano)
    dfcano["canon_smiles"] = cano_results.apply(lambda x: x[0])
    dfcano["_mol"] = cano_results.apply(lambda x: x[1])
    dfcanofr = dfcano.dropna(subset=["canon_smiles"]).copy()
    track("canonical conversion filter", len(dfvalu), len(dfcanofr))

    score_col = f"{rn}_r_i_docking_score"
    idx = dfcanofr.groupby("canon_smiles")[score_col].idxmin()

    df_best = dfcanofr.loc[idx].copy()
    df_best = df_best.reset_index(drop=True)
    track("duplicate removal (best per canon_smiles)", len(dfcanofr), len(df_best))

    print(f"tot_removed (ini -> final): {n0 - len(df_best)}")

    df_best = make_mask(df_best, mol_col="_mol")
    df_best = df_best.drop(columns=["_mol"])

    if sd is None:
        outpa = indir.replace(".csv", f"{rn}_cleaned.csv")
    else:
        outpa = sd

    df_best.to_csv(outpa, index=False)
    print(f"saved: {outpa}")

    logpa = os.path.splitext(outpa)[0] + "_cleanlog.txt"
    with open(logpa, "w") as f:
        f.write(f"input: {indir}\n")
        f.write(f"output: {outpa}\n")
        f.write(f"initial instance count: {n0}\n\n")
        for step_name, before, after, removed in steps:
            f.write(f"[{step_name}]\n")
            f.write(f"  before: {before}\n")
            f.write(f"  after:  {after}\n")
            f.write(f"  removed: {removed}\n\n")
        f.write(f"final instance count: {len(df_best)}\n")
        f.write(f"total removed: {n0 - len(df_best)}\n")
    print(f"log saved: {logpa}")

    return df_best

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
""" cleans dset steps... 
1. valid molecule filter 
2. unique molecule filter
3. canonical conversion filter(should be o)
4. dupli filter(canosmi based top dockscore)
5. compute violation mask
""",formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--indir", required=True, help="input CSV file path")
    parser.add_argument("--sd", default=None, help="save dir path Default will generate file in same dir as input+prefix")
    parser.add_argument("--rn", default=None, help="The string of Run name_def set to 'PPS'")
    args = parser.parse_args()
    df = opcleaning(args.indir, args.sd, args.rn)

# python 0_0_dclean.py --indir /home/andy/proj/701/code/developing_dir/tempdevtest/scores.csv


# patch for later::: addd loger in txt format porobs ...

