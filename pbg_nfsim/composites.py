"""Pre-built composite document factories for NFSim simulations."""

from pbg_nfsim.models.generate_flagella_bngl import (
    get_model_path,
    default_production_rates,
)


def make_complexation_document(
    model_file=None,
    n_steps=100,
    interval=50.0,
):
    """Create a composite document for standalone NFSim complexation.

    Returns a document dict ready for use with Composite().

    Args:
        model_file: Path to BNGL model file (defaults to flagella model)
        n_steps: Number of NFSim steps per interval
        interval: Time interval between process updates (seconds)

    Returns:
        dict: Composite document with NFSim process, stores, and emitter
    """
    if model_file is None:
        model_file = get_model_path()

    return {
        'complexation': {
            '_type': 'process',
            'address': 'local:nfsim',
            'config': {
                'model_file': model_file,
                'n_steps': n_steps,
            },
            'inputs': {
                'observables': ['species'],
            },
            'outputs': {
                'observables': ['species'],
            },
            'interval': interval,
        },
        'species': {},
        'emitter': {
            '_type': 'step',
            'address': 'local:ram-emitter',
            'config': {
                'emit': {
                    'species': 'map[float]',
                    'time': 'float',
                },
            },
            'inputs': {
                'species': ['species'],
                'time': ['global_time'],
            },
        },
    }


def make_production_document(
    model_file=None,
    n_steps=100,
    complexation_interval=50.0,
    production_interval=1.0,
    production_rate_scale=1.0,
):
    """Create a composite document for production + complexation.

    Composes MonomerProduction and NFSimProcess wired to shared species store.

    Args:
        model_file: Path to BNGL model file (defaults to flagella model)
        n_steps: Number of NFSim steps per complexation interval
        complexation_interval: Time interval for NFSim process (seconds)
        production_interval: Time interval for monomer production (seconds)
        production_rate_scale: Multiplier for production rates

    Returns:
        dict: Composite document with production, complexation, stores, emitter
    """
    if model_file is None:
        model_file = get_model_path()

    rates = default_production_rates()
    scaled_rates = {
        name: rate * production_rate_scale
        for name, rate in rates.items()
    }

    return {
        'production': {
            '_type': 'process',
            'address': 'local:monomer-production',
            'config': {
                'production_rates': scaled_rates,
            },
            'outputs': {
                'monomers': ['species'],
            },
            'interval': production_interval,
        },
        'complexation': {
            '_type': 'process',
            'address': 'local:nfsim',
            'config': {
                'model_file': model_file,
                'n_steps': n_steps,
            },
            'inputs': {
                'observables': ['species'],
            },
            'outputs': {
                'observables': ['species'],
            },
            'interval': complexation_interval,
        },
        'species': {},
        'emitter': {
            '_type': 'step',
            'address': 'local:ram-emitter',
            'config': {
                'emit': {
                    'species': 'map[float]',
                    'time': 'float',
                },
            },
            'inputs': {
                'species': ['species'],
                'time': ['global_time'],
            },
        },
    }
