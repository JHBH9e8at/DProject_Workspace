import argparse

from dask.distributed import Client

from fetch_target import build_target_records, inspect_target, write_manifest


"""Dask-compatible target-file inspection using fetch_target.py's core logic.
currently under development and testing.
Likely it will contain some buggs ! So dont use this fiel f now!
"""

def fetch_targets_dask(
    input_csv,
    run_dir,
    run_name,
    output_csv,
    cluster,
    step_col="step",
    variant_col=None,
    limit=None,
    allow_unavailable=False,
):
    records = build_target_records(
        input_csv=input_csv,
        run_dir=run_dir,
        run_name=run_name,
        step_col=step_col,
        variant_col=variant_col,
        limit=limit,
    )

    client = Client(cluster)
    print(f"Dask dashboard: {client.dashboard_link}")

    try:
        futures = client.map(inspect_target, records)
        checked_records = client.gather(futures)
    finally:
        client.close()

    return write_manifest(checked_records, output_csv, allow_unavailable)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Dask version of fetch_target.py"
    )
    parser.add_argument("--input", required=True, help="Filtered/selected CSV")
    parser.add_argument("--run-dir", required=True, help="AHC result directory")
    parser.add_argument("--run-name", required=True, help="Run prefix, e.g. PPS or PR")
    parser.add_argument("--output", required=True, help="Output manifest CSV")
    parser.add_argument("--cluster", default="tcp://138.37.52.153:8786")
    parser.add_argument("--step-col", default="step")
    parser.add_argument("--variant-col", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--allow-unavailable", action="store_true")
    args = parser.parse_args()

    fetch_targets_dask(
        input_csv=args.input,
        run_dir=args.run_dir,
        run_name=args.run_name,
        output_csv=args.output,
        cluster=args.cluster,
        step_col=args.step_col,
        variant_col=args.variant_col,
        limit=args.limit,
        allow_unavailable=args.allow_unavailable,
    )


# # exmpl

# python fetch_target_dask.py \
#   --input /path/to/PPS_filtered_selected.csv \
#   --run-dir /home/andy/proj/701/output/2026_06_02_SMILES-RNN_PPS_tester_run_2 \
#   --run-name PPS \
#   --output /path/to/PPS_to_PR_manifest.csv \
#   --cluster tcp://138.37.52.153:8786          #<- sure to check 