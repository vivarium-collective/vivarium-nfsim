"""
NFSim Processes
===============

Process-bigraph wrappers for BioNetGen/NFSim rule-based simulations.

NFSimProcess: runs a BNGL model with NFSim, returning delta changes
in observable values. Suitable for composition with other processes.

MonomerProduction: produces flagellar protein monomers at constant
rates, mimicking gene expression.
"""
import os
import re
import tempfile

import bionetgen

from process_bigraph.composite import Process


def _parse_bngl_text(bngl_text):
    """Parse a BNGL model text to extract observables, seed species, and molecule types.

    Returns:
        observable_names: list of observable names
        obs_to_pattern: dict mapping observable name -> BNGL pattern
        seed_pattern_to_param: dict mapping seed species pattern -> parameter name
        simple_molecule_types: set of molecule type names with no internal states
    """
    observable_names = []
    obs_to_pattern = {}
    seed_pattern_to_param = {}
    simple_molecule_types = set()

    # Parse observables block
    in_observables = False
    for line in bngl_text.splitlines():
        stripped = line.strip()
        if stripped == 'begin observables':
            in_observables = True
            continue
        if stripped == 'end observables':
            break
        if in_observables and stripped and not stripped.startswith('#'):
            parts = stripped.split()
            if len(parts) >= 3:
                name = parts[1]
                pattern = parts[2]
                observable_names.append(name)
                obs_to_pattern[name] = pattern

    # Parse seed species block
    in_seeds = False
    for line in bngl_text.splitlines():
        stripped = line.strip()
        if stripped == 'begin seed species':
            in_seeds = True
            continue
        if stripped == 'end seed species':
            break
        if in_seeds and stripped and not stripped.startswith('#'):
            parts = stripped.split()
            if len(parts) >= 2:
                pattern = parts[0]
                param = parts[1]
                seed_pattern_to_param[pattern] = param

    # Parse molecule types block
    in_mol_types = False
    for line in bngl_text.splitlines():
        stripped = line.strip()
        if stripped == 'begin molecule types':
            in_mol_types = True
            continue
        if stripped == 'end molecule types':
            break
        if in_mol_types and stripped and not stripped.startswith('#'):
            mol_type = stripped
            # Simple molecule type: Name() with no internal states (no ~)
            if '~' not in mol_type and mol_type.endswith('()'):
                name = mol_type[:-2]
                simple_molecule_types.add(name)

    return observable_names, obs_to_pattern, seed_pattern_to_param, simple_molecule_types


class NFSimProcess(Process):
    """A generic process that wraps BioNetGen/NFSim network-free simulations.

    Loads a BNGL model, runs NFSim for each time step, and returns
    delta changes in observable values. Seed species counts are set
    from the current state before each run, allowing composition with
    other processes.

    Observables whose molecule types have no internal states (simple
    molecules) are "seedable" -- their counts carry over between steps.
    Growing intermediates with counter states are not seedable and are
    lost between steps.
    """

    config_schema = {
        'model_file': 'string',
        'n_steps': {
            '_type': 'integer',
            '_default': 100,
        },
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self.model_file = self.config['model_file']

        # Read and store the BNGL template
        with open(self.model_file) as f:
            self.bngl_template = f.read()

        # Parse model structure
        (self.observable_names,
         self.obs_to_pattern,
         self.seed_pattern_to_param,
         self.simple_molecule_types) = _parse_bngl_text(self.bngl_template)

        # Build seedability mapping for each observable
        # seedable_obs maps observable_name -> ('param', param_name) or ('add', pattern)
        self.seedable_obs = {}
        for name in self.observable_names:
            pattern = self.obs_to_pattern.get(name, '')
            if pattern in self.seed_pattern_to_param:
                # Has an existing seed species with a parameter
                self.seedable_obs[name] = ('param', self.seed_pattern_to_param[pattern])
            elif pattern.endswith('()'):
                # Check if molecule type is simple (no internal states)
                mol_name = pattern[:-2]
                if mol_name in self.simple_molecule_types:
                    self.seedable_obs[name] = ('add', pattern)

    def initial_state(self):
        return {
            'observables': {name: 0.0 for name in self.observable_names},
        }

    def inputs(self):
        return {
            'observables': {
                name: 'float' for name in self.observable_names
            },
        }

    def outputs(self):
        return {
            'observables': {
                name: 'float' for name in self.observable_names
            },
        }

    def update(self, state, interval):
        # Read current observable values
        current = {
            name: max(0, int(state['observables'].get(name, 0)))
            for name in self.observable_names
        }

        # Skip if nothing to simulate
        total_seedable = sum(
            current[name] for name in self.observable_names
            if name in self.seedable_obs
        )
        if total_seedable == 0:
            return {
                'observables': {name: 0.0 for name in self.observable_names},
            }

        # Build BNGL text with current state
        bngl_text = self.bngl_template

        # Update existing seed species parameters
        for name in self.observable_names:
            if name not in self.seedable_obs:
                continue
            kind, ref = self.seedable_obs[name]
            count = current[name]
            if kind == 'param':
                # Replace parameter value (match within a single line only)
                pattern = rf'([ \t]+{re.escape(ref)}[ \t]+)\S+'
                bngl_text = re.sub(pattern, rf'\g<1>{count}', bngl_text)

        # Add new seed species for seedable observables without existing params
        extra_seeds = []
        for name in self.observable_names:
            if name not in self.seedable_obs:
                continue
            kind, ref = self.seedable_obs[name]
            count = current[name]
            if kind == 'add' and count > 0:
                extra_seeds.append(f'    {ref}  {count}')

        if extra_seeds:
            bngl_text = bngl_text.replace(
                'end seed species',
                '\n'.join(extra_seeds) + '\nend seed species')

        # Append simulate action
        bngl_text += (
            f'\nsimulate({{method=>"nf",'
            f't_end=>{interval},'
            f'n_steps=>{self.config["n_steps"]}}});\n'
        )

        # Run NFSim
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_bngl = os.path.join(tmpdir, 'model.bngl')
            with open(tmp_bngl, 'w') as f:
                f.write(bngl_text)

            result = bionetgen.run(tmp_bngl, out=tmpdir)
            gdat = result['model']

        # Compute deltas
        deltas = {}
        for name in self.observable_names:
            initial = current.get(name, 0)
            final = float(gdat[name][-1]) if name in gdat.dtype.names else 0.0
            deltas[name] = final - initial

        return {
            'observables': deltas,
        }


class MonomerProduction(Process):
    """Produces flagellar protein monomers at constant rates.

    Config:
        production_rates: dict mapping monomer name -> rate (molecules/second).
            Defaults to rates that produce ~1 flagellum worth per 100 seconds.
    """

    config_schema = {
        'production_rates': {
            '_type': 'map[float]',
            '_default': {},
        },
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self.rates = self.config['production_rates']

    def inputs(self):
        return {}

    def outputs(self):
        return {
            'monomers': {
                name: 'float' for name in self.rates
            },
        }

    def update(self, state, interval):
        return {
            'monomers': {
                name: rate * interval
                for name, rate in self.rates.items()
            },
        }
