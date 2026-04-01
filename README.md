# vivarium-nfsim

A [process-bigraph](https://github.com/vivarium-collective/process-bigraph) wrapper for [BioNetGen](https://bionetgen.org/)/[NFSim](https://github.com/RuleWorld/nfsim).

This package provides a `Process` that runs rule-based models written in BNGL using the NFSim network-free simulator, enabling their composition with other processes in the process-bigraph framework.

## Installation

```bash
pip install vivarium-nfsim
```

## Usage

```python
from process_bigraph import Composite, process_registry
from vivarium_nfsim.process import NFSimProcess

process_registry.register('nfsim', NFSimProcess)

instance = {
    'nfsim': {
        '_type': 'process',
        'address': 'local:nfsim',
        'config': {
            'model_file': 'path/to/model.bngl',
            't_end': 100,
            'n_steps': 50,
        },
        'wires': {
            'observables': ['observables_store'],
            'parameters': ['parameters_store'],
        }
    },
}

workflow = Composite({'state': instance})
workflow.run(100)
results = workflow.gather_results()
```

## Model Files

Place your `.bngl` model files in the `vivarium_nfsim/models/` directory or provide an absolute path in the config.

## Running Tests

```bash
pytest vivarium_nfsim/
```
