"""Unit tests for NFSimProcess and MonomerProduction."""

import pytest
from process_bigraph import allocate_core
from pbg_nfsim.processes import NFSimProcess, MonomerProduction
from pbg_nfsim.models.generate_flagella_bngl import get_model_path, default_production_rates


@pytest.fixture
def core():
    c = allocate_core()
    c.register_link('nfsim', NFSimProcess)
    c.register_link('monomer-production', MonomerProduction)
    return c


@pytest.fixture
def model_path():
    return get_model_path()


def test_nfsim_instantiation(core, model_path):
    proc = NFSimProcess(
        config={'model_file': model_path, 'n_steps': 10},
        core=core)
    assert proc.config['model_file'] == model_path
    assert proc.config['n_steps'] == 10
    assert len(proc.observable_names) > 0


def test_nfsim_initial_state(core, model_path):
    proc = NFSimProcess(
        config={'model_file': model_path, 'n_steps': 10},
        core=core)
    state = proc.initial_state()
    assert 'observables' in state
    # All observables should start at 0
    for name, val in state['observables'].items():
        assert val == 0.0, f'{name} should be 0.0, got {val}'


def test_nfsim_inputs_outputs(core, model_path):
    proc = NFSimProcess(
        config={'model_file': model_path, 'n_steps': 10},
        core=core)
    inputs = proc.inputs()
    outputs = proc.outputs()
    assert 'observables' in inputs
    assert 'observables' in outputs
    # Should have the same observable names
    assert set(inputs['observables'].keys()) == set(outputs['observables'].keys())
    assert len(inputs['observables']) == len(proc.observable_names)


def test_nfsim_seedability(core, model_path):
    proc = NFSimProcess(
        config={'model_file': model_path, 'n_steps': 10},
        core=core)
    # Should have seedable observables (free monomers)
    assert len(proc.seedable_obs) > 0
    # Growing intermediates should NOT be seedable
    for name in proc.observable_names:
        if name.startswith('Growing_'):
            assert name not in proc.seedable_obs


def test_nfsim_update_zero_state(core, model_path):
    """With zero state, update should return zero deltas (no simulation)."""
    proc = NFSimProcess(
        config={'model_file': model_path, 'n_steps': 10},
        core=core)
    state = proc.initial_state()
    result = proc.update(state, interval=10.0)
    assert 'observables' in result
    # All deltas should be zero since there's nothing to simulate
    for name, val in result['observables'].items():
        assert val == 0.0


def test_nfsim_config_defaults(core, model_path):
    proc = NFSimProcess(
        config={'model_file': model_path},
        core=core)
    assert proc.config['n_steps'] == 100


def test_monomer_production_instantiation(core):
    rates = {'Free_fliG': 1.3, 'Free_flgE': 6.0}
    proc = MonomerProduction(
        config={'production_rates': rates},
        core=core)
    assert proc.rates == rates


def test_monomer_production_outputs(core):
    rates = {'Free_fliG': 1.3, 'Free_flgE': 6.0}
    proc = MonomerProduction(
        config={'production_rates': rates},
        core=core)
    outputs = proc.outputs()
    assert 'monomers' in outputs
    assert 'Free_fliG' in outputs['monomers']
    assert 'Free_flgE' in outputs['monomers']


def test_monomer_production_update(core):
    rates = {'Free_fliG': 1.3, 'Free_flgE': 6.0}
    proc = MonomerProduction(
        config={'production_rates': rates},
        core=core)
    result = proc.update({}, interval=10.0)
    assert abs(result['monomers']['Free_fliG'] - 13.0) < 1e-6
    assert abs(result['monomers']['Free_flgE'] - 60.0) < 1e-6


def test_default_production_rates():
    rates = default_production_rates()
    assert len(rates) > 0
    # All rates should be positive
    for name, rate in rates.items():
        assert rate > 0, f'{name} rate should be positive'
        assert name.startswith('Free_'), f'{name} should start with Free_'
