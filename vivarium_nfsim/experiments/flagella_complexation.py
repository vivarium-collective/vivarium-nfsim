"""
Flagella Complexation Experiment
================================

Runs the E. coli flagella structural assembly model using NFSim.
Models the hierarchical assembly of a complete flagellum from protein
monomers, following the complexation pathway defined in vivarium-chemotaxis.

Assembly hierarchy:
    FlhDC (master regulator complex)
    Motor Switch (FliG/M/N C-ring)
    Export Apparatus (type III secretion system)
    Flagellar Motor (basal body + stator)
    Hook (FlgE polymer)
    Complete Flagellum (all subcomplexes + filament)
"""
import os
import tempfile

import bionetgen
import numpy as np
import matplotlib.pyplot as plt


MODEL_FILE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    'models', 'flagella_complexation.bngl')

OUTPUT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
    'output')


def run_nfsim(model_file=MODEL_FILE, t_end=500, n_steps=500):
    """Run the flagella complexation model with NFSim."""
    model = bionetgen.bngmodel(model_file)

    model.actions.clear_actions()
    model.add_action('simulate', {
        'method': '"nf"',
        't_end': t_end,
        'n_steps': n_steps,
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_bngl = os.path.join(tmpdir, 'model.bngl')
        with open(tmp_bngl, 'w') as f:
            f.write(str(model))

        result = bionetgen.run(tmp_bngl, out=tmpdir)
        gdat = result['model']

    return gdat


def plot_results(gdat, output_dir=OUTPUT_DIR):
    """Generate plots showing flagella assembly dynamics."""
    os.makedirs(output_dir, exist_ok=True)
    time = gdat['time']
    names = [n for n in gdat.dtype.names if n != 'time']

    # Color scheme
    colors = {
        'flhDC': '#1f77b4',
        'flagellar_motor_switch': '#ff7f0e',
        'flagellar_export_apparatus_subunit': '#2ca02c',
        'flagellar_export_apparatus': '#d62728',
        'flagellar_motor': '#9467bd',
        'flagellar_hook': '#8c564b',
        'flagella': '#e377c2',
    }

    # ---- Figure 1: Complex Assembly Hierarchy ----
    fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    fig.suptitle(
        'E. coli Flagellum Assembly from Protein Monomers\n'
        '(stoichiometry from vivarium-chemotaxis)',
        fontsize=14)

    # Panel A: Sub-complex formation
    ax = axes[0]
    ax.set_title('Sub-complex Formation')
    subcomplexes = [
        ('flhDC', 'FlhD$_4$C$_2$'),
        ('flagellar_motor_switch', 'Motor Switch (26 FliG + 34 FliM + FliN)'),
        ('flagellar_export_apparatus_subunit', 'Export Apparatus Subunit'),
        ('flagellar_export_apparatus', 'Export Apparatus (+ 12 FliH)'),
    ]
    for obs, label in subcomplexes:
        if obs in names:
            ax.plot(time, gdat[obs], label=label, color=colors.get(obs), linewidth=2)
    ax.set_ylabel('Complex count')
    ax.legend(loc='right', fontsize=9)
    ax.set_ylim(bottom=-0.5)

    # Panel B: Higher-order assembly
    ax = axes[1]
    ax.set_title('Higher-order Assembly')
    higher = [
        ('flagellar_motor', 'Flagellar Motor'),
        ('flagellar_hook', 'Hook (120 FlgE)'),
        ('flagella', 'Complete Flagellum'),
    ]
    for obs, label in higher:
        if obs in names:
            ax.plot(time, gdat[obs], label=label, color=colors.get(obs), linewidth=2)
    ax.set_ylabel('Complex count')
    ax.legend(loc='right', fontsize=9)
    ax.set_ylim(bottom=-0.5)

    # Panel C: Key monomer depletion
    ax = axes[2]
    ax.set_title('Rate-limiting Monomer Depletion')
    monomers = [
        ('Free_fliG', 'Free FliG (need 26/switch)', '-'),
        ('Free_fliM', 'Free FliM (need 34/switch)', '-'),
        ('Free_flgE', 'Free FlgE (need 120/hook)', '-'),
        ('Free_fliH', 'Free FliH (need 12/export)', '--'),
        ('Free_fliI', 'Free FliI (need 6/export)', '--'),
        ('Free_fliD', 'Free FliD (need 5/flagellum)', ':'),
    ]
    for obs, label, ls in monomers:
        if obs in names:
            ax.plot(time, gdat[obs], label=label, linestyle=ls, linewidth=1.5)
    ax.set_ylabel('Molecule count')
    ax.set_xlabel('Time (s)')
    ax.legend(loc='right', fontsize=8)

    plt.tight_layout()
    fig_path = os.path.join(output_dir, 'flagella_complexation.png')
    plt.savefig(fig_path, dpi=150)
    plt.close()
    print(f'Saved assembly plot to {fig_path}')

    # ---- Figure 2: Growing intermediates ----
    growing_obs = [n for n in names if n.startswith('Growing_')]
    if growing_obs:
        fig, ax = plt.subplots(figsize=(12, 5))
        ax.set_title('Growing Intermediates (partially assembled complexes)')
        for obs in growing_obs:
            label = obs.replace('Growing_', '').replace('_total', '')
            ax.plot(time, gdat[obs], label=label, linewidth=1.5)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Count')
        ax.legend(fontsize=8)
        plt.tight_layout()
        fig_path = os.path.join(output_dir, 'flagella_intermediates.png')
        plt.savefig(fig_path, dpi=150)
        plt.close()
        print(f'Saved intermediates plot to {fig_path}')

    # ---- Figure 3: All observables grid ----
    n_obs = len(names)
    n_cols = 5
    n_rows = (n_obs + n_cols - 1) // n_cols
    fig, axes_grid = plt.subplots(n_rows, n_cols, figsize=(20, 3 * n_rows))
    fig.suptitle('All Observables — Flagella Assembly', fontsize=14)
    axes_flat = axes_grid.flatten()

    for i, obs in enumerate(names):
        ax = axes_flat[i]
        ax.plot(time, gdat[obs], color='steelblue', linewidth=1)
        ax.set_title(obs, fontsize=8)
        ax.tick_params(labelsize=7)

    for i in range(n_obs, len(axes_flat)):
        axes_flat[i].set_visible(False)

    plt.tight_layout()
    fig_path = os.path.join(output_dir, 'flagella_all_observables.png')
    plt.savefig(fig_path, dpi=150)
    plt.close()
    print(f'Saved all observables plot to {fig_path}')


def main():
    print('Running E. coli flagella complexation model with NFSim...')
    print(f'Model: {MODEL_FILE}')
    print('This models the hierarchical assembly of flagella from monomers:')
    print('  FlhDC -> Motor Switch -> Export Apparatus -> Motor -> Hook -> Flagellum')

    gdat = run_nfsim(t_end=500, n_steps=500)

    obs_names = [n for n in gdat.dtype.names if n != 'time']
    print(f'\nSimulation complete. {len(obs_names)} observables, {len(gdat["time"])} time points.')

    # Print final counts of complexes
    complex_obs = [n for n in obs_names if not n.startswith('Free_') and not n.startswith('Growing_')]
    print('\nFinal complex counts:')
    for obs in complex_obs:
        print(f'  {obs}: {gdat[obs][-1]:.0f}')

    plot_results(gdat)
    print('\nDone.')


if __name__ == '__main__':
    main()
