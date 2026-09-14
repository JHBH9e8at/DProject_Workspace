# Completed cross-docking analysis

Run completed score and optional interaction analysis with:

```powershell
python -m Modules.Crossdocking.CompletedAnalysis.crossdocking_analysis --config ..\crossdocking.in
```

`multiprocessing` affects only PoseCheck and ProLIF. Glide cross-docking remains
sequential by design because concurrent external Glide jobs may exhaust licence
slots or destabilize the execution host.

With `posecheck_chunk_size=auto`, each receptor-state population uses:

```text
ceil(valid_pose_count / posecheck_workers)
```

The minimum chunk size is one. A positive integer can be supplied instead.
