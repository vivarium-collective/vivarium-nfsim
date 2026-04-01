"""
Flagella Assembly with Continuous Monomer Production
====================================================

Composes two processes:
  1. MonomerProduction — continuously produces flagellar protein monomers
     at rates mimicking gene expression
  2. FlagellaComplexation — runs NFSim to assemble monomers into flagella

Both processes are wired to a shared 'species' store. Production adds
monomers, complexation consumes them and produces completed complexes.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# Avoid circular import from local bigraph-viz
sys.modules.setdefault('bigraph_viz', type(sys)('bigraph_viz'))

from bigraph_schema import allocate_core
from process_bigraph.composite import Composite
from process_bigraph.types.process import register_types

from vivarium_nfsim.processes.monomer_production import (
    MonomerProduction, DEFAULT_RATES,
)
from vivarium_nfsim.processes.flagella_complexation import (
    FlagellaComplexation, MONOMER_OBS, COMPLEX_OBS,
)


OUTPUT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
    'output')


def make_core():
    """Create and configure the process-bigraph core."""
    core = allocate_core()
    core = register_types(core)
    core.register_link('monomer-production', MonomerProduction)
    core.register_link('flagella-complexation', FlagellaComplexation)
    return core


def build_composite(
        core,
        production_interval=1.0,
        complexation_interval=50.0,
        production_rate_scale=1.0,
):
    """Build the composed production + complexation system."""
    scaled_rates = {
        name: rate * production_rate_scale
        for name, rate in DEFAULT_RATES.items()
    }

    state_spec = {
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
            'address': 'local:flagella-complexation',
            'config': {
                'n_steps': 100,
            },
            'inputs': {
                'monomers': ['species'],
                'complexes': ['complexes'],
            },
            'outputs': {
                'monomers': ['species'],
                'complexes': ['complexes'],
            },
            'interval': complexation_interval,
        },
    }

    workflow = Composite(
        {'state': state_spec},
        core=core,
    )

    return workflow


def run_experiment(total_time=2000, complexation_interval=50.0, production_rate_scale=1.0):
    """Run the composed experiment."""
    print(f'Building composite: production + complexation')
    print(f'  Total time: {total_time}s')
    print(f'  Complexation interval: {complexation_interval}s')
    print(f'  Production rate scale: {production_rate_scale}x')

    core = make_core()
    workflow = build_composite(
        core,
        complexation_interval=complexation_interval,
        production_rate_scale=production_rate_scale,
    )

    print('Running simulation...')
    workflow.run(total_time)

    # Read final state
    final_state = workflow.state
    return final_state, workflow


def collect_timecourse(total_time, complexation_interval=50.0, production_rate_scale=1.0):
    """Run the experiment and collect state at regular intervals.

    Since the composite doesn't have a built-in emitter wired up,
    we run in segments and record the state after each segment.
    """
    core = make_core()

    # Segment duration = complexation interval (so we see each complexation step)
    segment = complexation_interval
    n_segments = int(total_time / segment)

    timecourse = {
        'time': [],
        'species': {name: [] for name in MONOMER_OBS},
        'complexes': {name: [] for name in COMPLEX_OBS},
    }

    print(f'Running {n_segments} segments of {segment}s each...')

    for i in range(n_segments):
        workflow = build_composite(
            core,
            complexation_interval=complexation_interval,
            production_rate_scale=production_rate_scale,
        )

        # Set initial state from previous segment's final state
        if i > 0:
            if 'species' not in workflow.state:
                workflow.state['species'] = {}
            for name in MONOMER_OBS:
                workflow.state['species'][name] = timecourse['species'][name][-1]
            if 'complexes' not in workflow.state:
                workflow.state['complexes'] = {}
            for name in COMPLEX_OBS:
                workflow.state['complexes'][name] = timecourse['complexes'][name][-1]

        workflow.run(segment)

        t = (i + 1) * segment
        timecourse['time'].append(t)

        # Record state
        for name in MONOMER_OBS:
            val = workflow.state.get('species', {}).get(name, 0.0)
            timecourse['species'][name].append(float(val) if val else 0.0)

        for name in COMPLEX_OBS:
            val = workflow.state.get('complexes', {}).get(name, 0.0)
            timecourse['complexes'][name].append(float(val) if val else 0.0)

        # Progress
        if (i + 1) % 10 == 0 or i == 0:
            flagella_count = timecourse['complexes'].get('flagella', [0])[-1]
            print(f'  t={t:.0f}s  flagella={flagella_count:.0f}')

    return timecourse


def plot_results(timecourse, output_dir=OUTPUT_DIR):
    """Plot the composed simulation results."""
    os.makedirs(output_dir, exist_ok=True)

    time = timecourse['time']

    fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    fig.suptitle(
        'Flagellum Assembly with Continuous Monomer Production\n'
        '(MonomerProduction + FlagellaComplexation composed via process-bigraph)',
        fontsize=13,
    )

    # Panel A: Complex counts
    ax = axes[0]
    ax.set_title('Completed Complexes (accumulated over time)')
    complex_colors = {
        'flhDC': '#1f77b4',
        'flagellar_motor_switch': '#ff7f0e',
        'flagellar_export_apparatus_subunit': '#2ca02c',
        'flagellar_export_apparatus': '#d62728',
        'flagellar_motor': '#9467bd',
        'flagellar_hook': '#8c564b',
        'flagella': '#e377c2',
    }
    complex_labels = {
        'flhDC': 'FlhD$_4$C$_2$',
        'flagellar_motor_switch': 'Motor Switch',
        'flagellar_export_apparatus_subunit': 'Export Subunit',
        'flagellar_export_apparatus': 'Export Apparatus',
        'flagellar_motor': 'Motor',
        'flagellar_hook': 'Hook',
        'flagella': 'Complete Flagellum',
    }
    for name in COMPLEX_OBS:
        data = timecourse['complexes'].get(name, [])
        if data:
            ax.plot(
                time, data,
                label=complex_labels.get(name, name),
                color=complex_colors.get(name),
                linewidth=2,
            )
    ax.set_ylabel('Count')
    ax.legend(loc='upper left', fontsize=9)
    ax.set_ylim(bottom=-0.5)

    # Panel B: Key monomer pools
    ax = axes[1]
    ax.set_title('Free Monomer Pools (production - consumption)')
    key_monomers = [
        ('Free_fliG', 'FliG (26/switch)'),
        ('Free_fliM', 'FliM (34/switch)'),
        ('Free_flgE', 'FlgE (120/hook)'),
        ('Free_fliH', 'FliH (12/export)'),
        ('Free_fliI', 'FliI (6/export)'),
        ('Free_fliD', 'FliD (5/flagellum)'),
    ]
    for obs_name, label in key_monomers:
        data = timecourse['species'].get(obs_name, [])
        if data:
            ax.plot(time, data, label=label, linewidth=1.5)
    ax.set_ylabel('Molecule count')
    ax.legend(loc='upper left', fontsize=8)

    # Panel C: Low-copy monomers
    ax = axes[2]
    ax.set_title('Low-copy Monomers')
    low_copy = [
        ('Free_flhD', 'FlhD (4/FlhDC)'),
        ('Free_flhC', 'FlhC (2/FlhDC)'),
        ('Free_flhA', 'FlhA (1/export)'),
        ('Free_flhB', 'FlhB (1/export)'),
        ('Free_fliN', 'FliN (1/switch)'),
        ('Free_fliE', 'FliE (1/motor)'),
        ('Free_motA', 'MotA (1/motor)'),
    ]
    for obs_name, label in low_copy:
        data = timecourse['species'].get(obs_name, [])
        if data:
            ax.plot(time, data, label=label, linewidth=1.5)
    ax.set_ylabel('Molecule count')
    ax.set_xlabel('Time (s)')
    ax.legend(loc='upper left', fontsize=8, ncol=2)

    plt.tight_layout()
    fig_path = os.path.join(output_dir, 'flagella_with_production.png')
    plt.savefig(fig_path, dpi=150)
    plt.close()
    print(f'Saved plot to {fig_path}')


def main():
    print('=' * 60)
    print('Flagella Assembly with Continuous Monomer Production')
    print('=' * 60)
    print()
    print('Two composed processes:')
    print('  1. MonomerProduction: adds monomers every 1s')
    print('  2. FlagellaComplexation: runs NFSim every 50s')
    print()

    timecourse = collect_timecourse(
        total_time=2000,
        complexation_interval=50.0,
        production_rate_scale=1.0,
    )

    plot_results(timecourse)

    # Print final counts
    print('\nFinal complex counts:')
    for name in COMPLEX_OBS:
        data = timecourse['complexes'].get(name, [0])
        print(f'  {name}: {data[-1]:.0f}')

    print('\nDone.')


if __name__ == '__main__':
    main()
