import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

desc = [
    "desc_QED", "desc_SAscore", "desc_CLogP", "desc_MolWt",
    "desc_TPSA", "desc_NumHAcceptors", "desc_NumHDonors"
]

REF_MAP = {"OM":  "OmecamtivMecarbil", "MAV": "Mavacamten", "AFI": "Aficamten"}
REF_COLORS = {"OM": "red", "MAV": "orange", "AFI": "green"}

def load_refs(ref_keys):
    if not ref_keys:
        return {}
    ref_df = pd.read_csv(r"/home/andy/proj/701/Analysis/ref_desc.csv")
    refs = {}
    for key in ref_keys:
        name = REF_MAP[key]
        row = ref_df[ref_df["name"] == name]
        refs[key] = row.iloc[0]
    return refs

def plot_distribution(df, jn, odir, refs={}):
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    for i, col in enumerate(desc):
        ax = axes[i]
        ax.hist(df[col].dropna(), bins=50, color="steelblue", edgecolor="none", alpha=0.7)

        for key, row in refs.items():
            if col in row and pd.notna(row[col]):
                ax.axvline(x=row[col], color=REF_COLORS[key], linewidth=1.5,
                           linestyle="--", label=REF_MAP[key])

        ax.set_title(col, fontsize=12)
        ax.set_xlabel("")
        ax.set_ylabel("Count")

    axes[-1].set_visible(False)

    if refs:
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(0.98, 0.02), fontsize=11)

    plt.suptitle(f"[{jn}] Physicochemical Properties Distribution", fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(f"{odir}/{jn}_physchem_distribution.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {odir}/{jn}_physchem_distribution.png")



def plot_scatter(df, jn, odir, threshold, refs={}):
    dc = f"{jn}_r_i_docking_score"

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    max_step = df["step"].max()
    min_last50 = max_step - 50
    last_50_mask = df["step"] >= min_last50

    for i, col in enumerate(desc):
        ax = axes[i]

        sample_all = df[[col, dc, "step"]].dropna()
        last50_mask_local = last_50_mask.reindex(sample_all.index, fill_value=False)

        good_mask = sample_all[dc] <= threshold
        weak_mask = sample_all[dc] > threshold

        for data, color, label in [
            (sample_all[~last50_mask_local], "steelblue", f"steps 0 ~ {min_last50 - 1}"),
            (sample_all[last50_mask_local], "tomato", f"steps {min_last50} ~ {max_step}"),
        ]:
            idx = data.index
            ax.scatter(
                data.loc[good_mask.reindex(idx, fill_value=False), col],
                data.loc[good_mask.reindex(idx, fill_value=False), dc],
                alpha=0.5, s=5, color=color, label=label
            )
            ax.scatter(
                data.loc[weak_mask.reindex(idx, fill_value=False), col],
                data.loc[weak_mask.reindex(idx, fill_value=False), dc],
                alpha=0.05, s=5, color=color
            )

        ax.axhline(y=threshold, color="black", linewidth=1, linestyle="--", alpha=0.7, label=f"docking = {threshold}")
        ax.set_title(col, fontsize=11)
        ax.set_xlabel(col)
        ax.set_ylabel(dc)

    axes[-1].set_visible(False)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(0.98, 0.02), fontsize=11)

    plt.suptitle(f"[{jn}] Physicochemical Properties vs Docking Score", fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(f"{odir}/{jn}_physchem_docking_scatter.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {odir}/{jn}_physchem_docking_scatter.png")



def mainr(ind, odir, jn, threshold, ref_keys):
    print(f"Loading {ind}")
    df = pd.read_csv(ind)
    refs = load_refs(ref_keys)

    plot_distribution(df, jn, odir, refs)
    plot_scatter(df, jn, odir, threshold)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""basic physiochem analysis..,
Always check Job name --jn !!!!!""",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--ind",  type=str, required=True, help="Input CSV path")
    parser.add_argument("--odir", type=str, required=True, help="Output directory")
    parser.add_argument("--jn",   type=str, required=True, help="Job name e.g. 'PPS'")
    parser.add_argument("--threshold", type=float, required=True, help="Docking score threshold (e.g. -9.0). good/weak split + axhline marker, used for scatter plot only")
    parser.add_argument("--ref",  type=str, nargs="*", choices=["OM", "MAV", "AFI"],
                        default=[], help="chose from ref mol OM MAV AFI")
    args = parser.parse_args()

    mainr(args.ind, args.odir, args.jn, args.threshold, args.ref)


#python 1_0_Physioprop.py --ind /home/andy/proj/701/code/developing_dir/tempdevtest/scores_cleaned.csv --odir /home/andy/proj/701/code/developing_dir/tempdevtest/figs --jn PPS --threshold -9.0 --ref OM



