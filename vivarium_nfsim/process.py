"""
NFSim Process
=============

A process-bigraph wrapper for BioNetGen/NFSim rule-based simulations.
"""
import os
import tempfile
import numpy as np
import bionetgen

from process_bigraph import Process, Composite, process_registry


class NFSimProcess(Process):
    """A process that wraps BioNetGen/NFSim network-free simulations.

    Loads a BNGL model, runs NFSim for each time step, and returns
    observable values as outputs.
    """

    config_schema = {
        'model_file': 'string',
        'n_steps': {
            '_type': 'integer',
            '_default': 10,
        },
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self.model_file = self.config['model_file']

        # load the model to extract structure
        model = bionetgen.bngmodel(self.model_file)
        self.parameter_names = [p.name for p in model.parameters]
        self.observable_names = [o.name for o in model.observables]
        self.species_names = [str(s.name) for s in model.species]

        # store initial parameter values
        self.default_parameters = {
            p.name: float(p.value) for p in model.parameters
        }

    def initial_state(self):
        return {
            'parameters': dict(self.default_parameters),
            'observables': {name: 0.0 for name in self.observable_names},
            'time': 0.0,
        }

    def inputs(self):
        return {
            'parameters': {
                name: 'float' for name in self.parameter_names
            },
            'time': 'float',
        }

    def outputs(self):
        return {
            'observables': {
                name: 'float' for name in self.observable_names
            },
            'time': 'float',
        }

    def update(self, state, interval):
        # reload the model so we can modify parameters
        model = bionetgen.bngmodel(self.model_file)

        # set parameters from input state
        for name, value in state['parameters'].items():
            if hasattr(model.parameters, name):
                setattr(model.parameters, name, value)

        # clear existing actions and add NFSim simulation
        model.actions.clear_actions()
        model.add_action(
            'simulate',
            {
                'method': '"nf"',
                't_start': state['time'],
                't_end': state['time'] + interval,
                'n_steps': self.config['n_steps'],
            }
        )

        # write modified model to temp file and run
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_bngl = os.path.join(tmpdir, 'model.bngl')
            with open(tmp_bngl, 'w') as f:
                f.write(str(model))

            result = bionetgen.run(tmp_bngl, out=tmpdir)

            # extract observable values from the last time point
            gdat = result['model']
            observables = {}
            for name in self.observable_names:
                if name in gdat.dtype.names:
                    observables[name] = float(gdat[name][-1])
                else:
                    observables[name] = 0.0

        return {
            'observables': observables,
            'time': interval,
        }


process_registry.register('nfsim', NFSimProcess)


def test_nfsim_process():
    """Test NFSim process with a simple BNGL model."""
    model_path = os.path.join(
        os.path.dirname(__file__), 'models', 'toggle_switch.bngl')

    if not os.path.exists(model_path):
        print(f'Model file not found: {model_path}')
        return

    instance = {
        'nfsim': {
            '_type': 'process',
            'address': 'local:nfsim',
            'config': {
                'model_file': model_path,
                'n_steps': 10,
            },
            'wires': {
                'observables': ['observables_store'],
                'parameters': ['parameters_store'],
                'time': ['time_store'],
            }
        },
        'emitter': {
            '_type': 'step',
            'address': 'local:ram-emitter',
            'config': {
                'ports': {
                    'inputs': {
                        'observables': 'tree[float]',
                    }
                }
            },
            'wires': {
                'inputs': {
                    'observables': ['observables_store'],
                }
            }
        }
    }

    workflow = Composite({'state': instance})
    workflow.run(10)
    results = workflow.gather_results()
    print(f'RESULTS: {results}')


if __name__ == '__main__':
    test_nfsim_process()
