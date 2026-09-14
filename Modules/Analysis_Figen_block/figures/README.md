# Standalone figure regeneration

```powershell
python run_figures.py --config figures.in
```

Use this runner when figures need to be regenerated independently. Normal AHC
and structural module runners already generate the figures coupled to their
analysis products, so running this module again is optional.

Available blocks are AHC physicochemical/trajectory, chemical-space,
physicochemical last-window, full/last hotspot, top-candidate hotspot, and
own/opposite hotspot. At least one must be ON. Only inputs for enabled blocks
are required.

Crossdocking density and selectivity figures are not duplicated here. Their
canonical implementation is `Modules/Crossdocking/CompletedAnalysis/figures`.
The old copies are preserved under `Analysis_Figen_block/bins/legacy_crossdocking_figures`.

All generated figures use collision-safe names: an existing file is retained
and the new figure receives `_1`, `_2`, and so on. Relative paths resolve from
`figures.in`.
