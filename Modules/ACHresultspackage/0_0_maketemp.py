import os
import glob
import pandas as pd
import argparse
from datetime import datetime


def gettempres(base_dir: str, odir: str, jn: str):
    pattern = os.path.join(base_dir, "*_scores.csv")
    csv_files = sorted(glob.glob(pattern))

    if len(csv_files) == 0:
        print("NOcsvrecheck_targ_basedir")
        return None, None

    if len(csv_files) == 1:
        print("Recheck_Targbasedir")
        return None, None

    files_to_merge = csv_files[:-1]
    max_iter = len(files_to_merge)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = os.path.join(odir, f"{jn}_{max_iter}_{timestamp}.csv")

    print("proxyrep")
    print(f"tot {len(csv_files)} ::::::: last file ({os.path.basename(csv_files[-1])}) not included")
    print(f"... merging {len(files_to_merge)}")

    dfs = [pd.read_csv(f) for f in files_to_merge]
    merged = pd.concat(dfs, ignore_index=True)

    if not os.path.isdir(odir):
        os.makedirs(odir, exist_ok=True)
    merged.to_csv(output_path, index=False)

    print(f"{output_path}::::({len(merged)} rows)")
    return merged, output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="***Merging iters in basedir_excluding last file***")
    parser.add_argument("--base_dir", type=str, required=True, help="source dir:souldbe smthing likie .../output/2026_06_02_SMILES-RNN_PPS_tester_run_2/iterations")
    parser.add_argument("--odir", type=str, required=True, help="Output directory for the merged temp csv")
    parser.add_argument("--jn", type=str, required=True, help="Job name e.g. 'PPS', used in output filename")
    args = parser.parse_args()

    gettempres(args.base_dir, args.odir, args.jn)