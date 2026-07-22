# Single-ligand Glide docking

This directory contains a small, server-side docking workflow. A prepared
Glide receptor grid replaces the original protein file at this stage.

## Files

- `prepare_ligand.py`: ligand file to LigPrep SDF
- `make_glide_input.py`: write the Glide `.in` file
- `run_glide.py`: execute `$SCHRODINGER/glide`
- `extract_scores.py`: extract pose properties and docking scores
- `single_dock.py`: run all four stages in order

## Complete run

Run inside the same `molscore` conda environment used for AHC:


# Run format
python single_dock.py \
  --ligand /path/to/ligand.smi \
  --grid /path/to/PR_Auto_Grid.zip \
  --output-dir /home/andy/proj/701/output/single_dock_test \
  --precision SP


Expected outputs include:
```text
prepared_ligand.sdf
single_docking.in
glide_console.log
single_docking_lib.sdfgz
docking_scores.csv
```

# Run one stage at a time
python prepare_ligand.py --ligand ligand.smi --output prepared_ligand.sdf

python make_glide_input.py \
  --ligand prepared_ligand.sdf \
  --grid PR_Auto_Grid.zip \
  --output single_docking.in

python run_glide.py --input single_docking.in --output-dir output

python extract_scores.py \
  --poses output/single_docking_lib.sdfgz \
  --output output/docking_scores.csv


This single-ligand workflow runs Glide directly. Dask can later call
`single_dock()` in parallel when scaling to many ligands.

