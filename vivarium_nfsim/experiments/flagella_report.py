"""
Flagella Complexation Report Generator
======================================

Generates an HTML report with detailed plots for both experiments:
  1. Standalone complexation (fixed monomer pool)
  2. Composed: continuous production + complexation

Saves the report to doc/flagella_report.html
"""
import os
import sys
import base64
import io

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Avoid circular import from local bigraph-viz
sys.modules.setdefault('bigraph_viz', type(sys)('bigraph_viz'))

from vivarium_nfsim.experiments.flagella_complexation import (
    run_nfsim, MODEL_FILE,
)
from vivarium_nfsim.experiments.flagella_with_production import (
    collect_timecourse, MONOMER_OBS, COMPLEX_OBS,
)
from vivarium_nfsim.models.generate_flagella_bngl import (
    COMPLEXATION_STOICHIOMETRY, _safe_name,
)

DOC_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
    'docs')


def fig_to_base64(fig):
    """Convert a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=140, bbox_inches='tight')
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return b64


# ── Color palette ───────────────────────────────────────────────────
COMPLEX_COLORS = {
    'flhDC': '#1f77b4',
    'flagellar_motor_switch': '#ff7f0e',
    'flagellar_export_apparatus_subunit': '#2ca02c',
    'flagellar_export_apparatus': '#d62728',
    'flagellar_motor': '#9467bd',
    'flagellar_hook': '#8c564b',
    'flagella': '#e377c2',
}
COMPLEX_LABELS = {
    'flhDC': 'FlhD\u2084C\u2082',
    'flagellar_motor_switch': 'Motor Switch',
    'flagellar_export_apparatus_subunit': 'Export Subunit',
    'flagellar_export_apparatus': 'Export Apparatus',
    'flagellar_motor': 'Motor',
    'flagellar_hook': 'Hook',
    'flagella': 'Complete Flagellum',
}


# ── Experiment 1: standalone complexation ───────────────────────────

def plot_standalone_assembly(gdat):
    """Three-panel figure for standalone assembly."""
    time = gdat['time']
    names = [n for n in gdat.dtype.names if n != 'time']

    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)

    # Panel A: completed complexes
    ax = axes[0]
    ax.set_title('Completed Complexes', fontweight='bold')
    for obs in COMPLEX_OBS:
        if obs in names:
            ax.step(time, gdat[obs], label=COMPLEX_LABELS.get(obs, obs),
                    color=COMPLEX_COLORS.get(obs), linewidth=2, where='post')
    ax.set_ylabel('Count')
    ax.legend(fontsize=9, loc='upper left')
    ax.set_ylim(bottom=-0.3)
    ax.grid(alpha=0.3)

    # Panel B: monomer depletion
    ax = axes[1]
    ax.set_title('Monomer Depletion', fontweight='bold')
    highlight = [
        ('Free_fliG', 'FliG (26/switch)'),
        ('Free_fliM', 'FliM (34/switch)'),
        ('Free_flgE', 'FlgE (120/hook)'),
        ('Free_fliH', 'FliH (12/export)'),
        ('Free_fliI', 'FliI (6/export)'),
        ('Free_fliD', 'FliD (5/flagellum)'),
    ]
    for obs, label in highlight:
        if obs in names:
            ax.plot(time, gdat[obs], label=label, linewidth=1.5)
    ax.set_ylabel('Molecule count')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(alpha=0.3)

    # Panel C: growing intermediates
    ax = axes[2]
    ax.set_title('Growing Intermediates (partially assembled scaffolds)', fontweight='bold')
    growing = [n for n in names if n.startswith('Growing_')]
    for obs in growing:
        label = obs.replace('Growing_', '').replace('_total', '')
        ax.plot(time, gdat[obs], label=label, linewidth=1.2)
    ax.set_ylabel('Count')
    ax.set_xlabel('Time (s)')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(alpha=0.3)

    fig.tight_layout()
    return fig


def plot_standalone_grid(gdat):
    """Grid of all observables."""
    time = gdat['time']
    names = [n for n in gdat.dtype.names if n != 'time']
    n = len(names)
    cols = 5
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(18, 2.8 * rows))
    for i, obs in enumerate(names):
        ax = axes.flat[i]
        ax.plot(time, gdat[obs], color='steelblue', linewidth=1)
        ax.set_title(obs, fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.2)
    for i in range(n, len(axes.flat)):
        axes.flat[i].set_visible(False)
    fig.tight_layout()
    return fig


# ── Experiment 2: composed production + complexation ────────────────

def plot_composed_overview(tc):
    """Three-panel overview for the composed experiment."""
    time = tc['time']

    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)

    # Panel A: accumulated complexes
    ax = axes[0]
    ax.set_title('Accumulated Complexes', fontweight='bold')
    for obs in COMPLEX_OBS:
        data = tc['complexes'].get(obs, [])
        if data:
            ax.step(time, data, label=COMPLEX_LABELS.get(obs, obs),
                    color=COMPLEX_COLORS.get(obs), linewidth=2, where='post')
    ax.set_ylabel('Count')
    ax.legend(fontsize=9, loc='upper left')
    ax.set_ylim(bottom=-0.3)
    ax.grid(alpha=0.3)

    # Panel B: high-copy monomers (sawtooth dynamics)
    ax = axes[1]
    ax.set_title('High-copy Monomer Pools (production \u2192 consumption cycles)', fontweight='bold')
    for obs, label in [
        ('Free_fliG', 'FliG'),
        ('Free_fliM', 'FliM'),
        ('Free_flgE', 'FlgE'),
        ('Free_fliH', 'FliH'),
        ('Free_fliI', 'FliI'),
    ]:
        data = tc['species'].get(obs, [])
        if data:
            ax.plot(time, data, label=label, linewidth=1.2)
    ax.set_ylabel('Molecule count')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.3)

    # Panel C: low-copy monomers
    ax = axes[2]
    ax.set_title('Low-copy Monomer Pools', fontweight='bold')
    for obs, label in [
        ('Free_flhD', 'FlhD'),
        ('Free_flhC', 'FlhC'),
        ('Free_fliN', 'FliN'),
        ('Free_flhA', 'FlhA'),
        ('Free_fliE', 'FliE'),
        ('Free_motA', 'MotA'),
        ('Free_fliD', 'FliD'),
    ]:
        data = tc['species'].get(obs, [])
        if data:
            ax.plot(time, data, label=label, linewidth=1.2)
    ax.set_ylabel('Molecule count')
    ax.set_xlabel('Time (s)')
    ax.legend(fontsize=8, loc='upper left', ncol=2)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    return fig


def plot_composed_detail(tc):
    """Detail: flagella production rate and monomer budget."""
    time = np.array(tc['time'])
    flagella = np.array(tc['complexes'].get('flagella', []))

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

    # Panel A: cumulative flagella with rate annotation
    ax = axes[0]
    ax.set_title('Flagellum Production Over Time', fontweight='bold')
    ax.step(time, flagella, color=COMPLEX_COLORS['flagella'], linewidth=2.5, where='post')
    ax.fill_between(time, flagella, step='post', alpha=0.15, color=COMPLEX_COLORS['flagella'])
    if len(time) > 1 and flagella[-1] > 0:
        rate = flagella[-1] / time[-1]
        ax.annotate(
            f'{flagella[-1]:.0f} flagella in {time[-1]:.0f}s\n'
            f'({rate * 60:.1f} per minute)',
            xy=(time[-1], flagella[-1]),
            xytext=(time[-1] * 0.55, flagella[-1] * 0.7),
            fontsize=11,
            arrowprops=dict(arrowstyle='->', color='gray'),
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow'),
        )
    ax.set_ylabel('Complete flagella')
    ax.grid(alpha=0.3)

    # Panel B: sub-complex inventory (stacked)
    ax = axes[1]
    ax.set_title('Sub-complex Inventory', fontweight='bold')
    sub_obs = [o for o in COMPLEX_OBS if o != 'flagella']
    for obs in sub_obs:
        data = tc['complexes'].get(obs, [])
        if data:
            ax.plot(time, data, label=COMPLEX_LABELS.get(obs, obs),
                    color=COMPLEX_COLORS.get(obs), linewidth=1.5)
    ax.set_ylabel('Count')
    ax.set_xlabel('Time (s)')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.3)

    fig.tight_layout()
    return fig


# ── Stoichiometry table ─────────────────────────────────────────────

def stoichiometry_table_html():
    rows = []
    rxn_order = [
        'flhDC',
        'flagellar motor switch reaction',
        'flagellar export apparatus reaction 1',
        'flagellar export apparatus reaction 2',
        'flagellar motor reaction',
        'flagellar hook reaction',
        'flagellum reaction',
    ]
    for i, rxn in enumerate(rxn_order, 1):
        stoich = COMPLEXATION_STOICHIOMETRY[rxn]
        product = None
        subs = []
        for sp, cnt in stoich.items():
            if cnt > 0:
                product = sp
            else:
                c = int(abs(cnt))
                subs.append(f'{c}&times; {sp}' if c > 1 else sp)
        rows.append(f'<tr><td>{i}</td><td><b>{product}</b></td><td>{" + ".join(subs)}</td></tr>')
    return '\n'.join(rows)


# ── HTML template ───────────────────────────────────────────────────

def build_html(images):
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Flagella Complexation Report &mdash; vivarium-nfsim</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
         max-width: 960px; margin: 2em auto; padding: 0 1em; color: #24292e; line-height: 1.6; }}
  h1 {{ border-bottom: 1px solid #e1e4e8; padding-bottom: 0.3em; }}
  h2 {{ border-bottom: 1px solid #eaecef; padding-bottom: 0.2em; margin-top: 2em; }}
  img.plot {{ width: 100%; border: 1px solid #e1e4e8; border-radius: 6px; margin: 1em 0; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
  th, td {{ border: 1px solid #dfe2e5; padding: 6px 12px; text-align: left; }}
  th {{ background: #f6f8fa; }}
  code {{ background: #f6f8fa; padding: 0.2em 0.4em; border-radius: 3px; font-size: 0.9em; }}
  pre {{ background: #f6f8fa; padding: 1em; border-radius: 6px; overflow-x: auto; }}
  .caption {{ color: #586069; font-size: 0.9em; margin-top: -0.5em; }}
  details {{ margin: 0.5em 0; }}
  summary {{ cursor: pointer; font-weight: 600; }}
</style>
</head>
<body>

<h1>E. coli Flagella Complexation &mdash; vivarium-nfsim</h1>

<p>
This report documents two experiments modeling the hierarchical assembly of the
<em>E. coli</em> flagellum from protein monomers, using
<a href="https://bionetgen.org/">BioNetGen</a>/<a href="https://github.com/RuleWorld/nfsim">NFSim</a>
wrapped as a <a href="https://github.com/vivarium-collective/process-bigraph">process-bigraph</a> process.
</p>

<h2>Biological Background</h2>

<p>
The bacterial flagellum is assembled from ~30 different proteins through seven
sequential complexation reactions. The stoichiometry is taken from
<a href="https://github.com/vivarium-collective/vivarium-chemotaxis/blob/master/chemotaxis/data/chromosomes/flagella_chromosome.py">vivarium-chemotaxis</a>:
</p>

<table>
<tr><th>#</th><th>Product</th><th>Subunits</th></tr>
{stoichiometry_table_html()}
</table>

<h2>Modeling Approach</h2>

<p>
BioNetGen only supports bimolecular reactions. Multi-subunit assembly
(up to 120 copies of FlgE for the hook) is modeled using a
<b>scaffold-with-counters</b> strategy:
</p>

<ol>
<li><b>Nucleation</b> &mdash; two monomers (the scarcest chosen first) react to create a
<code>Growing_*</code> scaffold molecule with internal state counters.</li>
<li><b>Sequential growth</b> &mdash; monomers bind one at a time, incrementing the
appropriate counter.</li>
<li><b>Completion</b> &mdash; when all counters reach their targets, the scaffold converts
to the finished complex.</li>
</ol>

<p>This generates <b>237 reaction rules</b> from 7 biological reactions. The BNGL model is
generated programmatically by
<a href="https://github.com/vivarium-collective/vivarium-nfsim/blob/main/vivarium_nfsim/models/generate_flagella_bngl.py"><code>generate_flagella_bngl.py</code></a>.
NFSim simulates these efficiently without enumerating the species space.</p>

<hr>

<h2>Experiment 1: Standalone Complexation</h2>

<p>
Starting with enough monomers for 5 flagella, NFSim runs for 500s.
Assembly proceeds hierarchically: small complexes form first, then combine
into motors, hooks, and finally complete flagella.
</p>

<img class="plot" src="data:image/png;base64,{images['standalone_assembly']}"
     alt="Standalone assembly dynamics">
<p class="caption">
Sub-complexes form within seconds (FlhDC) to minutes (motor switch, export apparatus).
Complete flagella appear once all sub-complexes and the 120-subunit hook are ready.
Stochastic timing from NFSim produces realistic single-molecule kinetics.
</p>

<details>
<summary>All 44 observables (click to expand)</summary>
<img class="plot" src="data:image/png;base64,{images['standalone_grid']}"
     alt="All observables grid">
</details>

<hr>

<h2>Experiment 2: Continuous Production + Complexation</h2>

<p>
Two processes composed via process-bigraph and wired to shared state:
</p>
<ul>
<li><b>MonomerProduction</b> &mdash; adds monomers every 1s at gene-expression rates
(~1 flagellum worth per 100s)</li>
<li><b>FlagellaComplexation</b> &mdash; runs NFSim every 50s, consuming accumulated
monomers and building complexes</li>
</ul>

<p>
Between NFSim steps, sub-complex counts carry over so that higher-order assembly
(motor + hook + export &rarr; flagellum) can proceed across steps.
</p>

<img class="plot" src="data:image/png;base64,{images['composed_overview']}"
     alt="Composed experiment overview">
<p class="caption">
Top: flagella accumulate roughly linearly as production feeds complexation.
Middle: high-copy monomer pools show sawtooth dynamics &mdash; production accumulates
them between complexation steps, which periodically consume them.
Bottom: low-copy monomers (1&ndash;5 per flagellum) hover near zero, indicating they
are rate-limiting.
</p>

<img class="plot" src="data:image/png;base64,{images['composed_detail']}"
     alt="Composed experiment detail">
<p class="caption">
Top: cumulative flagella production with annotated rate.
Bottom: sub-complex inventory shows transient accumulation of intermediates
(motor switches, export apparatus) before they are consumed by higher-order assembly.
</p>

<hr>

<h2>Running the Experiments</h2>

<pre><code># Generate the BNGL model
python -m vivarium_nfsim.models.generate_flagella_bngl

# Experiment 1: standalone complexation
python -m vivarium_nfsim.experiments.flagella_complexation

# Experiment 2: production + complexation composite
python -m vivarium_nfsim.experiments.flagella_with_production

# Regenerate this report
python -m vivarium_nfsim.experiments.flagella_report
</code></pre>

<hr>
<p style="color:#586069; font-size:0.85em;">
Generated by <a href="https://github.com/vivarium-collective/vivarium-nfsim">vivarium-nfsim</a>.
Source: <code>vivarium_nfsim/experiments/flagella_report.py</code>
</p>

</body>
</html>
"""


# ── Main ────────────────────────────────────────────────────────────

def main():
    os.makedirs(DOC_DIR, exist_ok=True)
    images = {}

    # ── Experiment 1 ──
    print('Running Experiment 1: standalone complexation...')
    gdat = run_nfsim(t_end=500, n_steps=500)
    print(f'  {len(gdat["time"])} time points')

    fig = plot_standalone_assembly(gdat)
    images['standalone_assembly'] = fig_to_base64(fig)

    fig = plot_standalone_grid(gdat)
    images['standalone_grid'] = fig_to_base64(fig)

    # ── Experiment 2 ──
    print('Running Experiment 2: production + complexation composite...')
    tc = collect_timecourse(
        total_time=2000,
        complexation_interval=50.0,
        production_rate_scale=1.0,
    )
    print(f'  {len(tc["time"])} segments, final flagella: '
          f'{tc["complexes"]["flagella"][-1]:.0f}')

    fig = plot_composed_overview(tc)
    images['composed_overview'] = fig_to_base64(fig)

    fig = plot_composed_detail(tc)
    images['composed_detail'] = fig_to_base64(fig)

    # ── Write HTML ──
    html = build_html(images)
    report_path = os.path.join(DOC_DIR, 'flagella_report.html')
    with open(report_path, 'w') as f:
        f.write(html)
    print(f'\nReport saved to {report_path}')


if __name__ == '__main__':
    main()
