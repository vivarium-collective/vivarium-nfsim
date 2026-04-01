# vivarium-nfsim

A [process-bigraph](https://github.com/vivarium-collective/process-bigraph) wrapper for [BioNetGen](https://bionetgen.org/)/[NFSim](https://github.com/RuleWorld/nfsim).

Wraps NFSim rule-based simulations as composable processes, enabling their wiring with other processes in the process-bigraph framework.

## Installation

```bash
pip install vivarium-nfsim
```

Requires [bionetgen](https://pypi.org/project/bionetgen/) (installs NFSim binaries automatically).

## Quick Start

```python
from bigraph_schema import allocate_core
from process_bigraph.composite import Composite
from process_bigraph.types.process import register_types
from vivarium_nfsim.process import NFSimProcess

core = allocate_core()
core = register_types(core)
core.register_link('nfsim', NFSimProcess)

workflow = Composite({'state': {
    'nfsim': {
        '_type': 'process',
        'address': 'local:nfsim',
        'config': {
            'model_file': 'path/to/model.bngl',
            'n_steps': 100,
        },
        'inputs': {
            'parameters': ['parameters'],
            'time': ['time'],
        },
        'outputs': {
            'observables': ['observables'],
            'time': ['time'],
        },
        'interval': 10.0,
    },
}}, core=core)

workflow.run(100)
```

## Flagella Complexation Demo

The included demo models hierarchical assembly of the *E. coli* flagellum from 30 protein monomers through 7 sequential complexation reactions, using stoichiometry from [vivarium-chemotaxis](https://github.com/vivarium-collective/vivarium-chemotaxis). A BNGL model with 237 rules is generated programmatically and simulated with NFSim.

Two experiments are provided:

1. **Standalone complexation** — fixed monomer pool assembles into ~5 flagella
2. **Production + complexation composite** — a `MonomerProduction` process continuously supplies monomers while `FlagellaComplexation` runs NFSim to assemble them, demonstrating process-bigraph composition

```bash
# generate the BNGL model from stoichiometry data
python -m vivarium_nfsim.models.generate_flagella_bngl

# run standalone complexation
python -m vivarium_nfsim.experiments.flagella_complexation

# run the composed production + complexation experiment
python -m vivarium_nfsim.experiments.flagella_with_production

# generate the full HTML report with all plots
python -m vivarium_nfsim.experiments.flagella_report
```

See the **[full report](https://vivarium-collective.github.io/vivarium-nfsim/)** for detailed results, plots, and explanation of the modeling approach.
