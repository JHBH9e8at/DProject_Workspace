# Institution

Queen Mary University of London School of Physical and Chemical Sciences

# Programme

MSc FT Chemistry Artificial Intelligence for Drug Discovery

# Description

MSc AI for Drug Discovery dissertation: structure-based generative modelling
of selective, state-specific cardiac myosin ligands using AHC-style RNN–RL and
MolScore.

## User-facing runners

Run each workflow through its own runner and `.in` configuration file. Scripts
inside `blocks/` are implementation details and should not normally be invoked
directly.

### Monitor an active AHC run

```bash
python Modules/ACHresultspackage/Monitoring/run_monitor.py \
  --config Modules/ACHresultspackage/Monitoring/t_mp.in
```

### Analyse completed AHC results

```bash
python Modules/ACHresultspackage/CompletedAnalysis/run_analysis.py \
  --config Modules/ACHresultspackage/CompletedAnalysis/analysis.in
```

### Run cross-docking

```bash
python Modules/Crossdocking/run_crossdocking.py \
  --config Modules/Crossdocking/crossdocking.in
```

### Run analysis and figure-generation modules

```bash
python Modules/Analysis_Figen_block/run_analysis.py \
  --config Modules/Analysis_Figen_block/analysis.in
```

Configure input data and output locations in the corresponding `.in` file
before running. New analysis outputs are written below the configured
`work_results` directory. The external AHC generation engine itself is not
bundled in this repository.
