import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
import numpy as np
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import os
import argparse
from datetime import datetime



DESC_COLS_ALL = [
    "desc_QED", "desc_SAscore", "desc_CLogP", "desc_MolWt",
    "desc_HeavyAtomCount", "desc_NumHAcceptors", "desc_NumHDonors",
    "desc_NumHeteroatoms", "desc_NumRotatableBonds", "desc_NumAromaticRings",
    "desc_NumAliphaticRings", "desc_RingCount", "desc_TPSA", "desc_PenLogP",
    "desc_FormalCharge", "desc_Bertz", "desc_MaxConsecutiveRotatableBonds",
    "desc_FlourineCount"
]

DESC_COLS_SELECT = [
    "desc_QED", "desc_SAscore", "desc_CLogP", "desc_MolWt",
    "desc_TPSA", "desc_NumHAcceptors", "desc_NumHDonors"
]

REF_MAP = {"OM": "OmecamtivMecarbil", "MAV": "Mavacamten", "AFI": "Aficamten"}

REF_MOLS = {
    "OmecamtivMecarbil": "COC(O)N1CCN(C[C@H]2CCC[C@@H](NC(O)N[C@@H]3CC[C@H](C)NC3)C2F)CC1",
    "Mavacamten": "C[C@@H](C1=CC=CC=C1)NC2=CC(=O)N(C(=O)N2)C(C)C",
    "Aficamten": "CCC1=NC(=NO1)C2=CC3=C(C=C2)[C@@H](CC3)NC(=O)C4=CN(N=C4)C"
}

REF_COLORS = {"OmecamtivMecarbil": "red", "Mavacamten": "orange", "Aficamten": "green"}
REF_MARKERS = {"OmecamtivMecarbil": "*", "Mavacamten": "^", "Aficamten": "D"}        #change show setting
REF_SIZE = 100


def resolve_refs(ref_keys):
    """ref_keys: list of 'OM'/'MAV'/'AFI' (any subset, any order). Empty/None -> no refs."""
    if not ref_keys:
        return []
    return [REF_MAP[k] for k in ref_keys]


def fpgen(sm):
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=16384)
    mol = Chem.MolFromSmiles(sm)
    if mol is None:
        return None
    fp = mfpgen.GetFingerprint(mol)
    arr = np.zeros((16384,), dtype=np.int8)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def _plot_tsne(em, color_vals, ref_ems, active_refs, title, odir, fname):
    fig, ax = plt.subplots(figsize=(8, 8))
    sc = ax.scatter(em[:, 0], em[:, 1], c=color_vals, s=4, alpha=0.5)
    plt.colorbar(sc, ax=ax)
    if ref_ems is not None:
        for name in active_refs:
            em_ref = ref_ems[name]
            ax.scatter(em_ref[0], em_ref[1],
                       c=REF_COLORS[name], s=REF_SIZE,
                       marker=REF_MARKERS[name], label=name, zorder=5)
        ax.legend()
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(os.path.join(odir, fname), dpi=150, bbox_inches="tight")
    # plt.show()


def gen_tsne_fp(ind: str, odir: str, jn: str, ref_ind: str = None, ref_keys=None):
    base = os.path.splitext(os.path.basename(ind))[0]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    active_refs = resolve_refs(ref_keys)

    df = pd.read_csv(ind)
    dfc = df[df[f"{jn}_r_i_docking_score"] != 0].copy()
    print(f"filtered: {len(dfc)}")

    fpdf = dfc.copy()
    fpdf["fp"] = dfc["canon_smiles"].apply(fpgen)
    fpdfc = fpdf.dropna(subset=["fp"]).copy()
    print(f"invalid conversions: {len(fpdf) - len(fpdfc)}")

    X = np.stack(fpdfc["fp"])

    ref_ems = None
    if ref_ind is not None and active_refs:
        ref_fps = []
        for name in active_refs:
            canon = Chem.MolToSmiles(Chem.MolFromSmiles(REF_MOLS[name]), canonical=True)
            ref_fps.append(fpgen(canon))

        X_combined = np.vstack([X, np.array(ref_fps)])

        tsne = TSNE(n_components=2, perplexity=30, random_state=42)
        em_combined = tsne.fit_transform(X_combined)

        em = em_combined[:-len(active_refs)]
        ref_ems = {name: em_combined[-(len(active_refs)-i)] for i, name in enumerate(active_refs)}
        for name, em_ref in ref_ems.items():
            print(f"{name}: {em_ref}")
    else:
        tsne = TSNE(n_components=2, perplexity=30, random_state=42)
        em = tsne.fit_transform(X)

    _plot_tsne(em, fpdfc[f"{jn}_r_i_docking_score"], ref_ems, active_refs,
               f"t-SNE (FP) - {jn} Docking Score (w/ ref)",
               odir, f"{base}_{jn}_tsne_fp_score_{timestamp}.png")

    _plot_tsne(em, fpdfc["step"], ref_ems, active_refs,
               "t-SNE (FP) - Step (w/ ref)",
               odir, f"{base}_{jn}_tsne_fp_step_{timestamp}.png")

def _tsne_desc_core(dfc, jn, odir, base, timestamp, desc_cols, label, active_refs, ref_desc_df=None):
    X = dfc[desc_cols].values

    ref_ems = None
    if ref_desc_df is not None and active_refs:
        ref_rows = ref_desc_df[ref_desc_df["name"].isin(active_refs)]
        ref_rows = ref_rows.set_index("name").loc[active_refs]
        ref_X = ref_rows[desc_cols].values

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        ref_X_scaled = scaler.transform(ref_X)
        X_combined = np.vstack([X_scaled, ref_X_scaled])
    else:
        scaler = StandardScaler()
        X_combined = scaler.fit_transform(X)

    tsne = TSNE(n_components=2, perplexity=30, random_state=42)
    em_combined = tsne.fit_transform(X_combined)

    if ref_desc_df is not None and active_refs:
        em = em_combined[:-len(active_refs)]
        ref_ems = {name: em_combined[-(len(active_refs)-i)] for i, name in enumerate(active_refs)}
        for name, em_ref in ref_ems.items():
            print(f"{name}: {em_ref}")
    else:
        em = em_combined

    print(f"embedding shape ({label}): {em.shape}")

    _plot_tsne(em, dfc[f"{jn}_r_i_docking_score"], ref_ems, active_refs,
               f"t-SNE ({label}) - {jn} Docking Score (w/ ref)",
               odir, f"{base}_{jn}_tsne_{label}_score_{timestamp}.png")

    _plot_tsne(em, dfc["step"], ref_ems, active_refs,
               f"t-SNE ({label}) - Step (w/ ref)",
               odir, f"{base}_{jn}_tsne_{label}_step_{timestamp}.png")


def gen_tsne_desc(ind: str, odir: str, jn: str, ref_ind: str = None, mode: str = "both", ref_keys=None):
    base = os.path.splitext(os.path.basename(ind))[0]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    active_refs = resolve_refs(ref_keys)

    df = pd.read_csv(ind)
    dfc = df[df[f"{jn}_r_i_docking_score"] != 0].copy()
    print(f"filtered: {len(dfc)}")

    ref_desc_df = pd.read_csv(ref_ind) if (ref_ind is not None and active_refs) else None

    if mode in ("all", "both"):
        _tsne_desc_core(dfc, jn, odir, base, timestamp, DESC_COLS_ALL, "all_desc", active_refs, ref_desc_df)
    if mode in ("select", "both"):
        _tsne_desc_core(dfc, jn, odir, base, timestamp, DESC_COLS_SELECT, "sel_desc", active_refs, ref_desc_df)

def run_tsne(ind: str, odir: str, jn: str, ref_ind: str = None, mode: str = "both",
             ref_keys=None, feat: str = "both"):
    """feat: 'fp', 'desc', or 'both' -- selects which t-SNE variant(s) to run."""
    if feat in ("fp", "both"):
        gen_tsne_fp(ind, odir, jn, ref_ind, ref_keys)
    if feat in ("desc", "both"):
        gen_tsne_desc(ind, odir, jn, ref_ind, mode, ref_keys)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate t-SNE plots. Please check jn !!!!!! Default set to PPS")
    parser.add_argument("--ind", type=str, required=True,  help="Input CSV path")
    parser.add_argument("--odir", type=str, required=True,  help="Output directory, use /home/andy/proj/701/Analysis/Figures")
    parser.add_argument("--jn", type=str, default="PPS",  help="Jobname ie 'PPS'")
    parser.add_argument("--ref_ind", type=str, default="/home/andy/proj/701/Analysis/ref_desc.csv", help="ref calced: /home/andy/proj/701/Analysis/ref_desc.csv")
    parser.add_argument("--mode", type=str, default="both", choices=["all", "select", "both"],
                        help="Descriptor mode (default: both)... all = desc all, select = selected_desc, both=plot both")
    parser.add_argument("--ref", type=str, nargs="*", choices=["OM", "MAV", "AFI"],
                        default=[], help="chose any subset of ref mol OM MAV AFI (default: none)")
    parser.add_argument("--feat", type=str, default="both", choices=["fp", "desc", "both"],
                        help="Which t-SNE variant(s) to run: fp, desc, or both (default: both)")
    args = parser.parse_args()

    run_tsne(args.ind, args.odir, args.jn, args.ref_ind, args.mode, args.ref, args.feat)



# python 1_2_tsne.py --ind /home/andy/proj/701/code/developing_dir/tempdevtest/scoresPPS_cleaned.csv --odir /home/andy/proj/701/code/developing_dir/tempdevtest/figs --ref OM --feat fp