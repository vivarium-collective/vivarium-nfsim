"""
Monomer Production Process
==========================

A process-bigraph Process that produces flagellar protein monomers at
configured rates, mimicking gene expression. Intended to be composed
with the flagella complexation process.

Production rates are loosely based on the flagellar gene regulatory
hierarchy (Kalir & Alon, 2004):
  - Class 1 (flhD, flhC): constitutive, low rate
  - Class 2 (motor, export, hook structural genes): activated by FlhDC
  - Class 3 (filament, chemotaxis genes): activated by FliA
"""
from process_bigraph.composite import Process

from vivarium_nfsim.models.generate_flagella_bngl import (
    COMPLEXATION_STOICHIOMETRY,
)


# Default production rates (molecules per second) for each monomer.
# Scaled so that ~1 flagellum worth of monomers is produced per ~100s.
def _default_production_rates():
    """Compute default rates: produce enough monomers for 1 flagellum per 100s."""
    duration = 100.0  # seconds to produce one flagellum worth

    # Collect total demand per monomer across all reactions
    demand = {}
    complex_names = set()
    for stoich in COMPLEXATION_STOICHIOMETRY.values():
        for species, count in stoich.items():
            if count > 0:
                complex_names.add(species)
            else:
                demand[species] = demand.get(species, 0) + abs(count)

    # Only produce true monomers (not sub-complexes)
    rates = {}
    for species, count in demand.items():
        if species not in complex_names:
            safe = species.replace(' ', '_').replace('-', '_')
            rates[f'Free_{safe}'] = count / duration

    return rates


DEFAULT_RATES = _default_production_rates()


class MonomerProduction(Process):
    """Produces flagellar protein monomers at constant rates.

    Config:
        production_rates: dict mapping monomer name -> rate (molecules/second).
            Defaults to rates that produce ~1 flagellum worth per 100 seconds.
    """

    config_schema = {
        'production_rates': {
            '_type': 'map[float]',
            '_default': DEFAULT_RATES,
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


def register(core):
    core.register_link('monomer-production', MonomerProduction)
    return core
