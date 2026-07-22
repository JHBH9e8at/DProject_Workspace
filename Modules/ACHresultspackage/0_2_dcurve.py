import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
import numpy as np
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import os
import argparse
from datetime import datetime

def gen_learn(ind:str, odir:str, jn:str):
    base = os.path.splitext(os.path.basename(ind))[0]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = os.path.join(odir, f"Docking_Score_Trajectory_{jn}_{base}_{timestamp}.png")
    
    df = pd.read_csv(ind)
    zdf=df[df[f"{jn}_r_i_docking_score"] == 0].groupby("step").size()
    dff=df[df[f"{jn}_r_i_docking_score"] != 0].copy()
    print(f"tot docking fail: {zdf.sum()}")

    gped=dff.groupby("step")[f"{jn}_r_i_docking_score"]
    xbar=gped.mean()
    std=gped.std()
    sdu=xbar+std
    sdd=xbar-std

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    ax1.plot(xbar.index, xbar.values, linewidth=1.5, label="mean")
    ax1.fill_between(xbar.index, sdu, sdd, alpha=0.3, label="±1 std")
    ax1.set_xlabel("Step")
    ax1.tick_params(labelbottom=True)
    ax1.set_ylabel("Mean Doccking score")

    ax2.bar(zdf.index, zdf.values, color="tomato", alpha=0.6, label="zero count")
    ax2.set_ylabel("Zero count")
    ax2.set_xlabel("Step")
    ax2.legend()

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2)


    fig.suptitle("Docking Score Trajectory", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.show()

            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate docking score trajectory plot.")
    parser.add_argument("--ind", type=str, required=True, help="Input CSV path")
    parser.add_argument("--odir", type=str, required=True, help="Output [[directory]] use /home/andy/proj/701/Analysis/Figures")
    parser.add_argument("--jn", type=str, required=True, help="Jobname ie 'PPS'")
    args = parser.parse_args()

    gen_learn(args.ind, args.odir, args.jn)


    #example
    #python 0_2_dcurve.py --ind /home/andy/proj/701/code/developing_dir/tempdevtest/scoresPPS_cleaned.csv --odir /home/andy/proj/701/code/developing_dir/tempdevtest/figs --jn PPS