# Institution:
Queen Mary University of London School of Physical and Chemical Sciences

# Programme:
MSc FT Chemistry Artificial Intelligence for Drug Discovery

# Description
 MSc AI for Drug Discovery dissertation: structure-based generative modelling of selective, state-specific cardiac myosin ligands using AHC style RNN–RL and MolScore.

## Block runners

Use one repository entry point instead of invoking individual scripts directly:

```bash
python -m Modules.runner crossdocking --config crossdocking.in
python -m Modules.runner analysis --workflow redocking-validation -- --help
python -m Modules.runner analysis --plan analysis_plan.json
```

The Crossdocking block runs the complete calculation pipeline and rejects an
output directory outside `work_results`. The Analysis block can run one of the
registered workflows or execute an ordered JSON plan. Existing detailed
entrypoints remain available for compatibility.

The external AHC generation engine is not bundled in this repository; AHC
result analysis is available through the Analysis block.
