"""
Flagella Complexation Process
=============================

A process-bigraph Process that runs the flagella complexation BNGL model
with NFSim, reading monomer counts from shared state and returning
the changes (monomers consumed, complexes produced).

Each update:
  1. Reads current monomer counts from the 'monomers' port
  2. Writes a temporary BNGL model with those counts as seed species
  3. Runs NFSim for the time interval
  4. Computes deltas: how many monomers were consumed, how many
     complexes were produced
  5. Returns the deltas as the update

Growing intermediates that don't complete within a single step are lost
(their monomers are effectively consumed). This represents failed
assembly attempts when subunits are scarce.
"""
import os
import tempfile

import bionetgen

from process_bigraph.composite import Process

from vivarium_nfsim.models.generate_flagella_bngl import (
    COMPLEXATION_STOICHIOMETRY,
    generate_bngl,
    _safe_name,
    _parse_reaction,
)


# Names of completed complexes
COMPLEX_NAMES = []
MONOMER_NAMES = []

_complex_set = set()
_all_consumed = set()
for _stoich in COMPLEXATION_STOICHIOMETRY.values():
    for _sp, _c in _stoich.items():
        if _c > 0:
            _complex_set.add(_sp)
        else:
            _all_consumed.add(_sp)

MONOMER_NAMES = sorted(_all_consumed - _complex_set)
COMPLEX_NAMES = [
    _sp for _stoich in COMPLEXATION_STOICHIOMETRY.values()
    for _sp, _c in _stoich.items() if _c > 0
]
# Deduplicate while preserving order
_seen = set()
COMPLEX_NAMES = [
    x for x in COMPLEX_NAMES
    if x not in _seen and not _seen.add(x)
]

# Observable names used in the BNGL
MONOMER_OBS = [f'Free_{_safe_name(m)}' for m in MONOMER_NAMES]
COMPLEX_OBS = [_safe_name(c) for c in COMPLEX_NAMES]


class FlagellaComplexation(Process):
    """Runs flagella complexation via NFSim.

    Reads free monomer counts from shared state, runs NFSim, and returns
    the change in monomer and complex counts.

    Config:
        n_steps: number of NFSim steps per update interval (default 100)
    """

    config_schema = {
        'n_steps': {
            '_type': 'integer',
            '_default': 100,
        },
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        self.monomer_obs = list(MONOMER_OBS)
        self.complex_obs = list(COMPLEX_OBS)

        # Map from observable name -> BNGL seed species pattern
        self.obs_to_seed = {}
        for monomer in MONOMER_NAMES:
            obs_name = f'Free_{_safe_name(monomer)}'
            seed_name = f'{_safe_name(monomer)}()'
            self.obs_to_seed[obs_name] = seed_name

    def inputs(self):
        schema = {}
        schema['monomers'] = {
            name: 'float' for name in self.monomer_obs
        }
        schema['complexes'] = {
            name: 'float' for name in self.complex_obs
        }
        return schema

    def outputs(self):
        schema = {}
        schema['monomers'] = {
            name: 'float' for name in self.monomer_obs
        }
        schema['complexes'] = {
            name: 'float' for name in self.complex_obs
        }
        return schema

    def initial_state(self):
        return {
            'monomers': {name: 0.0 for name in self.monomer_obs},
            'complexes': {name: 0.0 for name in self.complex_obs},
        }

    def update(self, state, interval):
        import re

        # Read current monomer counts (floor to integers for NFSim)
        monomer_counts = {
            name: max(0, int(state['monomers'].get(name, 0)))
            for name in self.monomer_obs
        }

        # Read current complex counts
        complex_counts = {
            name: max(0, int(state['complexes'].get(name, 0)))
            for name in self.complex_obs
        }

        # Skip if nothing available
        total = sum(monomer_counts.values()) + sum(complex_counts.values())
        if total == 0:
            return {
                'monomers': {n: 0.0 for n in self.monomer_obs},
                'complexes': {n: 0.0 for n in self.complex_obs},
            }

        # Generate BNGL with current counts as initial species.
        bngl_text = generate_bngl(n_flagella=1)

        # Replace monomer seed species counts
        for obs_name, count in monomer_counts.items():
            monomer_safe = obs_name.replace('Free_', '')
            param_name = f'{monomer_safe}_0'
            pattern = rf'(\s+{re.escape(param_name)}\s+)\d+'
            bngl_text = re.sub(pattern, rf'\g<1>{count}', bngl_text)

        # Add sub-complex seed species so they can participate in
        # higher-order assembly reactions within this NFSim step
        extra_seeds = []
        for obs_name, count in complex_counts.items():
            if count > 0:
                extra_seeds.append(f'    {obs_name}()  {count}')

        if extra_seeds:
            bngl_text = bngl_text.replace(
                'end seed species',
                '\n'.join(extra_seeds) + '\nend seed species')

        # Add simulate action
        bngl_text += f'\nsimulate({{method=>"nf",t_end=>{interval},n_steps=>{self.config["n_steps"]}}});\n'

        # Run NFSim
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_bngl = os.path.join(tmpdir, 'model.bngl')
            with open(tmp_bngl, 'w') as f:
                f.write(bngl_text)

            result = bionetgen.run(tmp_bngl, out=tmpdir)
            gdat = result['model']

        # Compute deltas
        monomer_deltas = {}
        for name in self.monomer_obs:
            initial = monomer_counts.get(name, 0)
            final = float(gdat[name][-1]) if name in gdat.dtype.names else 0.0
            monomer_deltas[name] = final - initial  # negative = consumed

        complex_deltas = {}
        for name in self.complex_obs:
            initial = complex_counts.get(name, 0)
            final = float(gdat[name][-1]) if name in gdat.dtype.names else 0.0
            complex_deltas[name] = final - initial  # new complexes formed

        return {
            'monomers': monomer_deltas,
            'complexes': complex_deltas,
        }


def register(core):
    core.register_link('flagella-complexation', FlagellaComplexation)
    return core
