"""Integration tests for NFSim composites."""

import pytest
from process_bigraph import Composite, allocate_core
from process_bigraph.emitter import RAMEmitter
from pbg_nfsim.processes import NFSimProcess, MonomerProduction
from pbg_nfsim.composites import make_complexation_document, make_production_document


@pytest.fixture
def core():
    c = allocate_core()
    c.register_link('nfsim', NFSimProcess)
    c.register_link('monomer-production', MonomerProduction)
    c.register_link('ram-emitter', RAMEmitter)
    return c


def test_complexation_assembly(core):
    doc = make_complexation_document(n_steps=10, interval=50.0)
    sim = Composite({'state': doc}, core=core)
    assert sim is not None


def test_production_assembly(core):
    doc = make_production_document(
        n_steps=10,
        complexation_interval=50.0,
        production_interval=1.0)
    sim = Composite({'state': doc}, core=core)
    assert sim is not None


def test_document_factory_params(core):
    doc = make_complexation_document(n_steps=50, interval=25.0)
    assert doc['complexation']['config']['n_steps'] == 50
    assert doc['complexation']['interval'] == 25.0


def test_production_document_rates(core):
    doc = make_production_document(production_rate_scale=2.0)
    rates = doc['production']['config']['production_rates']
    assert len(rates) > 0
    # All rates should be positive (scaled by 2.0)
    for name, rate in rates.items():
        assert rate > 0
