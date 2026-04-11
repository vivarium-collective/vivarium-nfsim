"""Demo: NFSim multi-configuration flagella assembly report.

Runs three distinct flagella assembly simulations and generates an
interactive HTML report with Plotly charts, bigraph-viz architecture
diagrams, and navigatable PBG document trees.

Configurations:
  1. Standalone Assembly — fixed monomer pool, shows hierarchical assembly
  2. Composed: Steady Production — monomer production + complexation
  3. Composed: Fast Production — higher production rate, rapid assembly
"""

import json
import os
import sys
import io
import base64
import tempfile
import time as _time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from process_bigraph import allocate_core
from process_bigraph.composite import Composite
from process_bigraph.emitter import RAMEmitter

from pbg_nfsim.processes import NFSimProcess, MonomerProduction, _parse_bngl_text
from pbg_nfsim.composites import make_complexation_document, make_production_document
from pbg_nfsim.models.generate_flagella_bngl import (
    get_model_path, COMPLEXATION_STOICHIOMETRY, default_production_rates,
)


# ── Observable classification ──────────────────────────────────────

def _classify_observables(model_path):
    """Classify observables into monomers, complexes, and growing."""
    with open(model_path) as f:
        bngl_text = f.read()
    obs_names, _, _, _ = _parse_bngl_text(bngl_text)
    monomer_obs = [n for n in obs_names if n.startswith('Free_')]
    complex_obs = [
        n for n in obs_names
        if not n.startswith('Free_') and not n.startswith('Growing_')
    ]
    growing_obs = [n for n in obs_names if n.startswith('Growing_')]
    return obs_names, monomer_obs, complex_obs, growing_obs


# ── Color palette ──────────────────────────────────────────────────

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
    'flhDC': 'FlhD4C2',
    'flagellar_motor_switch': 'Motor Switch',
    'flagellar_export_apparatus_subunit': 'Export Subunit',
    'flagellar_export_apparatus': 'Export Apparatus',
    'flagellar_motor': 'Motor',
    'flagellar_hook': 'Hook',
    'flagella': 'Complete Flagellum',
}

COLOR_SCHEMES = {
    'indigo': {'primary': '#6366f1', 'light': '#e0e7ff', 'dark': '#4338ca',
               'bg': '#eef2ff', 'accent': '#818cf8', 'text': '#312e81'},
    'emerald': {'primary': '#10b981', 'light': '#d1fae5', 'dark': '#059669',
                'bg': '#ecfdf5', 'accent': '#34d399', 'text': '#064e3b'},
    'rose': {'primary': '#f43f5e', 'light': '#ffe4e6', 'dark': '#e11d48',
             'bg': '#fff1f2', 'accent': '#fb7185', 'text': '#881337'},
}


# ── Simulation configs ─────────────────────────────────────────────

CONFIGS = [
    {
        'id': 'standalone',
        'title': 'Standalone Assembly',
        'subtitle': 'Fixed monomer pool, hierarchical flagella assembly',
        'description': (
            'Starting with enough monomers for 5 flagella (N=5), NFSim '
            'simulates 500 seconds of stochastic rule-based assembly. '
            'Complexes form hierarchically: FlhDC first (smallest, 6 subunits), '
            'then motor switch (61 subunits), export apparatus, motor, the '
            '120-subunit hook, and finally complete flagella. This demonstrates '
            'the NFSimProcess bridge pattern running standalone.'
        ),
        'mode': 'standalone',
        'n_flagella': 5,
        'total_time': 500.0,
        'n_steps_nfsim': 500,
        'color_scheme': 'indigo',
    },
    {
        'id': 'production_normal',
        'title': 'Composed: Steady Production',
        'subtitle': 'Continuous monomer production + complexation via process-bigraph',
        'description': (
            'Two processes composed via process-bigraph wired to shared state: '
            'MonomerProduction adds monomers every 1s at gene-expression rates '
            '(~1 flagellum worth per 100s), while NFSimProcess runs NFSim every '
            '50s to assemble accumulated monomers. Monomer pools show sawtooth '
            'dynamics as production accumulates them between complexation steps. '
            'This is the core workflow demonstrating process composition.'
        ),
        'mode': 'composed',
        'total_time': 2000.0,
        'complexation_interval': 50.0,
        'production_rate_scale': 1.0,
        'color_scheme': 'emerald',
    },
    {
        'id': 'production_fast',
        'title': 'Composed: Fast Production',
        'subtitle': 'Accelerated monomer production for rapid assembly',
        'description': (
            'Same composed workflow but with 3x monomer production rate, '
            'simulating overexpression of flagellar genes. The higher supply '
            'rate leads to faster assembly kinetics and larger monomer pools '
            'between complexation steps. Compare with the steady production '
            'configuration to see how production rate affects assembly yield.'
        ),
        'mode': 'composed',
        'total_time': 1500.0,
        'complexation_interval': 50.0,
        'production_rate_scale': 3.0,
        'color_scheme': 'rose',
    },
]


# ── Simulation runners ─────────────────────────────────────────────

def run_standalone(cfg, model_path, obs_names, monomer_obs, complex_obs, growing_obs):
    """Run standalone NFSim directly (not via Composite) and collect timecourse."""
    import bionetgen

    t0 = _time.perf_counter()

    model = bionetgen.bngmodel(model_path)
    model.actions.clear_actions()
    model.add_action('simulate', {
        'method': '"nf"',
        't_end': cfg['total_time'],
        'n_steps': cfg['n_steps_nfsim'],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_bngl = os.path.join(tmpdir, 'model.bngl')
        with open(tmp_bngl, 'w') as f:
            f.write(str(model))
        result = bionetgen.run(tmp_bngl, out=tmpdir)
        gdat = result['model']

    runtime = _time.perf_counter() - t0

    # Convert to timecourse dict
    time_arr = list(gdat['time'].astype(float))
    timecourse = {
        'time': time_arr,
        'complexes': {},
        'monomers': {},
        'growing': {},
    }
    for obs in complex_obs:
        if obs in gdat.dtype.names:
            timecourse['complexes'][obs] = list(gdat[obs].astype(float))
    for obs in monomer_obs:
        if obs in gdat.dtype.names:
            timecourse['monomers'][obs] = list(gdat[obs].astype(float))
    for obs in growing_obs:
        if obs in gdat.dtype.names:
            timecourse['growing'][obs] = list(gdat[obs].astype(float))

    return timecourse, runtime


def run_composed(cfg, model_path, obs_names, monomer_obs, complex_obs, growing_obs):
    """Run composed production+complexation and collect timecourse."""
    t0 = _time.perf_counter()

    core = allocate_core()
    core.register_link('nfsim', NFSimProcess)
    core.register_link('monomer-production', MonomerProduction)
    core.register_link('ram-emitter', RAMEmitter)

    segment = cfg['complexation_interval']
    total_time = cfg['total_time']
    n_segments = int(total_time / segment)
    production_rate_scale = cfg['production_rate_scale']

    timecourse = {
        'time': [],
        'complexes': {name: [] for name in complex_obs},
        'monomers': {name: [] for name in monomer_obs},
        'growing': {},
    }

    for i in range(n_segments):
        doc = make_production_document(
            model_file=model_path,
            n_steps=100,
            complexation_interval=segment,
            production_interval=1.0,
            production_rate_scale=production_rate_scale,
        )
        workflow = Composite({'state': doc}, core=core)

        # Set initial state from previous segment
        if i > 0:
            if 'species' not in workflow.state:
                workflow.state['species'] = {}
            for name in monomer_obs:
                workflow.state['species'][name] = timecourse['monomers'][name][-1]
            for name in complex_obs:
                workflow.state['species'][name] = timecourse['complexes'][name][-1]

        workflow.run(segment)

        t = (i + 1) * segment
        timecourse['time'].append(t)

        for name in monomer_obs:
            val = workflow.state.get('species', {}).get(name, 0.0)
            timecourse['monomers'][name].append(float(val) if val else 0.0)

        for name in complex_obs:
            val = workflow.state.get('species', {}).get(name, 0.0)
            timecourse['complexes'][name].append(float(val) if val else 0.0)

        if (i + 1) % 10 == 0 or i == 0:
            flagella_count = timecourse['complexes'].get('flagella', [0])[-1]
            print(f'    t={t:.0f}s  flagella={flagella_count:.0f}')

    runtime = _time.perf_counter() - t0
    return timecourse, runtime


# ── Bigraph diagram ────────────────────────────────────────────────

def generate_bigraph_image(cfg):
    """Generate a colored bigraph-viz PNG for the composite architecture."""
    from bigraph_viz import plot_bigraph

    if cfg['mode'] == 'standalone':
        doc = {
            'complexation': {
                '_type': 'process',
                'address': 'local:nfsim',
                'inputs': {'observables': ['species']},
                'outputs': {'observables': ['species']},
            },
            'species': {},
            'emitter': {
                '_type': 'step',
                'address': 'local:ram-emitter',
                'inputs': {
                    'species': ['species'],
                    'time': ['global_time'],
                },
            },
        }
        node_colors = {
            ('complexation',): '#6366f1',
            ('emitter',): '#8b5cf6',
            ('species',): '#e0e7ff',
        }
    else:
        doc = {
            'production': {
                '_type': 'process',
                'address': 'local:monomer-production',
                'outputs': {'monomers': ['species']},
            },
            'complexation': {
                '_type': 'process',
                'address': 'local:nfsim',
                'inputs': {'observables': ['species']},
                'outputs': {'observables': ['species']},
            },
            'species': {},
            'emitter': {
                '_type': 'step',
                'address': 'local:ram-emitter',
                'inputs': {
                    'species': ['species'],
                    'time': ['global_time'],
                },
            },
        }
        node_colors = {
            ('production',): '#10b981',
            ('complexation',): '#6366f1',
            ('emitter',): '#8b5cf6',
            ('species',): '#e0e7ff',
        }

    outdir = tempfile.mkdtemp()
    plot_bigraph(
        state=doc,
        out_dir=outdir,
        filename='bigraph',
        file_format='png',
        remove_process_place_edges=True,
        rankdir='LR',
        node_fill_colors=node_colors,
        node_label_size='16pt',
        port_labels=False,
        dpi='150',
    )
    png_path = os.path.join(outdir, 'bigraph.png')
    with open(png_path, 'rb') as f:
        b64 = base64.b64encode(f.read()).decode()
    return f'data:image/png;base64,{b64}'


def build_pbg_document(cfg):
    """Build the PBG composite document dict for display."""
    if cfg['mode'] == 'standalone':
        return make_complexation_document(
            n_steps=cfg.get('n_steps_nfsim', 100),
            interval=cfg['total_time'],
        )
    else:
        return make_production_document(
            n_steps=100,
            complexation_interval=cfg['complexation_interval'],
            production_rate_scale=cfg['production_rate_scale'],
        )


# ── Stoichiometry table ───────────────────────────────────────────

def stoichiometry_table_rows():
    rxn_order = [
        'flhDC',
        'flagellar motor switch reaction',
        'flagellar export apparatus reaction 1',
        'flagellar export apparatus reaction 2',
        'flagellar motor reaction',
        'flagellar hook reaction',
        'flagellum reaction',
    ]
    rows = []
    for i, rxn in enumerate(rxn_order, 1):
        stoich = COMPLEXATION_STOICHIOMETRY[rxn]
        product = None
        subs = []
        total = 0
        for sp, cnt in stoich.items():
            if cnt > 0:
                product = sp
            else:
                c = int(abs(cnt))
                total += c
                subs.append(f'{c}&times; {sp}' if c > 1 else sp)
        rows.append(
            f'<tr><td>{i}</td><td><b>{product}</b></td>'
            f'<td>{" + ".join(subs)}</td><td>{total}</td></tr>'
        )
    return '\n'.join(rows)


# ── HTML report generation ─────────────────────────────────────────

def generate_html(sim_results, output_path):
    """Generate comprehensive HTML report."""

    sections_html = []
    all_chart_data = {}

    for idx, (cfg, timecourse, runtime) in enumerate(sim_results):
        sid = cfg['id']
        cs = COLOR_SCHEMES[cfg['color_scheme']]

        times = timecourse['time']
        complexes = timecourse['complexes']
        monomers = timecourse['monomers']
        growing = timecourse.get('growing', {})

        # Metrics
        n_timepoints = len(times)
        final_flagella = complexes.get('flagella', [0])[-1] if complexes.get('flagella') else 0
        total_complexes = sum(
            data[-1] for data in complexes.values() if data
        )

        # Chart data for Plotly
        all_chart_data[sid] = {
            'times': times,
            'complexes': {k: v for k, v in complexes.items() if v},
            'monomers': {k: v for k, v in monomers.items() if v},
            'growing': {k: v for k, v in growing.items() if v},
        }

        # Bigraph diagram
        print(f'  Generating bigraph diagram for {sid}...')
        bigraph_img = generate_bigraph_image(cfg)

        # PBG document
        pbg_doc = build_pbg_document(cfg)

        # Key monomer highlights
        key_monomers_js = json.dumps([
            ('Free_fliG', 'FliG (26/switch)', '#1f77b4'),
            ('Free_fliM', 'FliM (34/switch)', '#ff7f0e'),
            ('Free_flgE', 'FlgE (120/hook)', '#2ca02c'),
            ('Free_fliH', 'FliH (12/export)', '#d62728'),
            ('Free_fliI', 'FliI (6/export)', '#9467bd'),
            ('Free_fliD', 'FliD (5/flagellum)', '#8c564b'),
        ])

        # Mode-specific description
        mode_label = 'Standalone NFSim' if cfg['mode'] == 'standalone' else 'Composed Workflow'

        section = f"""
    <div class="sim-section" id="sim-{sid}">
      <div class="sim-header" style="border-left: 4px solid {cs['primary']};">
        <div class="sim-number" style="background:{cs['light']}; color:{cs['dark']};">{idx+1}</div>
        <div>
          <h2 class="sim-title">{cfg['title']}</h2>
          <p class="sim-subtitle">{cfg['subtitle']}</p>
        </div>
      </div>
      <p class="sim-description">{cfg['description']}</p>

      <div class="metrics-row">
        <div class="metric"><span class="metric-label">Mode</span><span class="metric-value" style="font-size:1rem;">{mode_label}</span></div>
        <div class="metric"><span class="metric-label">Flagella</span><span class="metric-value">{final_flagella:.0f}</span></div>
        <div class="metric"><span class="metric-label">Total Complexes</span><span class="metric-value">{total_complexes:.0f}</span></div>
        <div class="metric"><span class="metric-label">Time Points</span><span class="metric-value">{n_timepoints}</span></div>
        <div class="metric"><span class="metric-label">Sim Duration</span><span class="metric-value">{cfg['total_time']:.0f}s</span></div>
        <div class="metric"><span class="metric-label">Runtime</span><span class="metric-value">{runtime:.1f}s</span></div>
      </div>

      <h3 class="subsection-title">Assembly Dynamics</h3>
      <div class="charts-row">
        <div class="chart-box"><div id="chart-complexes-{sid}" class="chart"></div></div>
        <div class="chart-box"><div id="chart-flagella-{sid}" class="chart"></div></div>
        <div class="chart-box"><div id="chart-monomers-{sid}" class="chart"></div></div>
        <div class="chart-box"><div id="chart-monomers2-{sid}" class="chart"></div></div>
      </div>

      <div class="pbg-row">
        <div class="pbg-col">
          <h3 class="subsection-title">Bigraph Architecture</h3>
          <div class="bigraph-img-wrap">
            <img src="{bigraph_img}" alt="Bigraph architecture diagram">
          </div>
        </div>
        <div class="pbg-col">
          <h3 class="subsection-title">Composite Document</h3>
          <div class="json-tree" id="json-{sid}"></div>
        </div>
      </div>
    </div>
"""
        sections_html.append(section)

    # Navigation
    nav_items = ''.join(
        f'<a href="#sim-{c["id"]}" class="nav-link" '
        f'style="border-color:{COLOR_SCHEMES[c["color_scheme"]]["primary"]};">'
        f'{c["title"]}</a>'
        for c in [r[0] for r in sim_results])

    # PBG docs for JSON viewer
    pbg_docs = {r[0]['id']: build_pbg_document(r[0]) for r in sim_results}

    # Complex color/label data for JS
    complex_colors_js = json.dumps(COMPLEX_COLORS)
    complex_labels_js = json.dumps(COMPLEX_LABELS)

    # Stoichiometry table
    stoich_rows = stoichiometry_table_rows()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Flagella Assembly Report &mdash; pbg-nfsim</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;
       background:#fff; color:#1e293b; line-height:1.6; }}
.page-header {{
  background:linear-gradient(135deg,#f8fafc 0%,#eef2ff 50%,#ecfdf5 100%);
  border-bottom:1px solid #e2e8f0; padding:3rem;
}}
.page-header h1 {{ font-size:2.2rem; font-weight:800; color:#0f172a; margin-bottom:.3rem; }}
.page-header p {{ color:#64748b; font-size:.95rem; max-width:700px; }}
.nav {{ display:flex; gap:.8rem; padding:1rem 3rem; background:#f8fafc;
        border-bottom:1px solid #e2e8f0; position:sticky; top:0; z-index:100;
        flex-wrap:wrap; }}
.nav-link {{ padding:.4rem 1rem; border-radius:8px; border:1.5px solid;
             text-decoration:none; font-size:.85rem; font-weight:600;
             transition:all .15s; color:#334155; }}
.nav-link:hover {{ transform:translateY(-1px); box-shadow:0 2px 8px rgba(0,0,0,.08); }}

/* Biology section */
.bio-section {{ padding:2rem 3rem; background:#f8fafc; border-bottom:1px solid #e2e8f0; }}
.bio-section h2 {{ font-size:1.3rem; font-weight:700; color:#0f172a; margin-bottom:1rem; }}
.bio-section p {{ color:#475569; font-size:.9rem; max-width:800px; margin-bottom:.8rem; }}
.bio-section ol {{ color:#475569; font-size:.88rem; padding-left:1.5rem; margin-bottom:1rem; }}
.bio-section li {{ margin-bottom:.3rem; }}
table.stoich {{ border-collapse:collapse; width:100%; max-width:700px; margin:1rem 0; font-size:.85rem; }}
table.stoich th, table.stoich td {{ border:1px solid #e2e8f0; padding:6px 12px; text-align:left; }}
table.stoich th {{ background:#eef2ff; font-weight:600; color:#334155; }}
table.stoich td {{ color:#475569; }}

.sim-section {{ padding:2.5rem 3rem; border-bottom:1px solid #e2e8f0; }}
.sim-header {{ display:flex; align-items:center; gap:1rem; margin-bottom:.8rem;
               padding-left:1rem; }}
.sim-number {{ width:36px; height:36px; border-radius:10px; display:flex;
               align-items:center; justify-content:center; font-weight:800; font-size:1.1rem; }}
.sim-title {{ font-size:1.5rem; font-weight:700; color:#0f172a; }}
.sim-subtitle {{ font-size:.9rem; color:#64748b; }}
.sim-description {{ color:#475569; font-size:.9rem; margin-bottom:1.5rem; max-width:800px; }}
.subsection-title {{ font-size:1.05rem; font-weight:600; color:#334155;
                     margin:1.5rem 0 .8rem; }}
.metrics-row {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr));
                gap:.8rem; margin-bottom:1.5rem; }}
.metric {{ background:#f8fafc; border:1px solid #e2e8f0; border-radius:10px;
           padding:.8rem; text-align:center; }}
.metric-label {{ display:block; font-size:.7rem; text-transform:uppercase;
                 letter-spacing:.06em; color:#94a3b8; margin-bottom:.2rem; }}
.metric-value {{ display:block; font-size:1.3rem; font-weight:700; color:#1e293b; }}
.metric-sub {{ display:block; font-size:.7rem; color:#94a3b8; }}
.charts-row {{ display:grid; grid-template-columns:1fr 1fr; gap:1rem; margin-bottom:1rem; }}
.chart-box {{ background:#f8fafc; border:1px solid #e2e8f0; border-radius:10px; overflow:hidden; }}
.chart {{ height:320px; }}
.pbg-row {{ display:grid; grid-template-columns:1fr 1fr; gap:1.5rem; margin-top:1rem; }}
.pbg-col {{ min-width:0; }}
.bigraph-img-wrap {{ background:#fafafa; border:1px solid #e2e8f0; border-radius:10px;
                     padding:1.5rem; text-align:center; }}
.bigraph-img-wrap img {{ max-width:100%; height:auto; }}
.json-tree {{ background:#f8fafc; border:1px solid #e2e8f0; border-radius:10px;
              padding:1rem; max-height:500px; overflow-y:auto; font-family:'SF Mono',
              Menlo,Monaco,'Courier New',monospace; font-size:.78rem; line-height:1.5; }}
.jt-key {{ color:#7c3aed; font-weight:600; }}
.jt-str {{ color:#059669; }}
.jt-num {{ color:#2563eb; }}
.jt-bool {{ color:#d97706; }}
.jt-null {{ color:#94a3b8; }}
.jt-toggle {{ cursor:pointer; user-select:none; color:#94a3b8; margin-right:.3rem; }}
.jt-toggle:hover {{ color:#1e293b; }}
.jt-collapsed {{ display:none; }}
.jt-bracket {{ color:#64748b; }}
.footer {{ text-align:center; padding:2rem; color:#94a3b8; font-size:.8rem;
           border-top:1px solid #e2e8f0; }}
@media(max-width:900px) {{
  .charts-row,.pbg-row {{ grid-template-columns:1fr; }}
  .sim-section,.page-header,.bio-section {{ padding:1.5rem; }}
}}
</style>
</head>
<body>

<div class="page-header">
  <h1>E. coli Flagella Assembly Report</h1>
  <p>Three flagella assembly simulations wrapped as <strong>process-bigraph</strong>
  Processes using <a href="https://bionetgen.org/" style="color:#6366f1;">BioNetGen</a>/<a href="https://github.com/RuleWorld/nfsim" style="color:#6366f1;">NFSim</a>
  network-free stochastic simulation. Each configuration demonstrates a distinct
  assembly scenario with interactive visualization.</p>
</div>

<div class="nav">{nav_items}</div>

<div class="bio-section">
  <h2>Biological Background</h2>
  <p>
  The bacterial flagellum is assembled from ~30 different proteins through seven
  sequential complexation reactions. Multi-subunit assembly (up to 120 copies of
  FlgE for the hook) is modeled using a <b>scaffold-with-counters</b> strategy
  that generates <b>237 reaction rules</b> from 7 biological reactions.
  </p>

  <table class="stoich">
  <tr><th>#</th><th>Product</th><th>Subunits</th><th>Total</th></tr>
  {stoich_rows}
  </table>

  <p>
  NFSim simulates these rules efficiently without enumerating the species space,
  using a network-free stochastic approach that scales well for large rule sets.
  </p>

  <h2>Modeling Approach</h2>
  <ol>
  <li><b>Nucleation</b> &mdash; two monomers react to create a Growing scaffold with internal counters</li>
  <li><b>Sequential growth</b> &mdash; monomers bind one at a time, incrementing counters</li>
  <li><b>Completion</b> &mdash; when all counters reach targets, the scaffold converts to the finished complex</li>
  </ol>
</div>

{''.join(sections_html)}

<div class="footer">
  Generated by <strong>pbg-nfsim</strong> &mdash;
  BioNetGen/NFSim + process-bigraph &mdash;
  Network-Free Stochastic Rule-Based Simulation
</div>

<script>
const DATA = {json.dumps(all_chart_data)};
const DOCS = {json.dumps(pbg_docs, indent=2, default=str)};
const COMPLEX_COLORS = {complex_colors_js};
const COMPLEX_LABELS = {complex_labels_js};

// ─── JSON Tree Viewer ───
function renderJson(obj, depth) {{
  if (depth === undefined) depth = 0;
  if (obj === null) return '<span class="jt-null">null</span>';
  if (typeof obj === 'boolean') return '<span class="jt-bool">' + obj + '</span>';
  if (typeof obj === 'number') return '<span class="jt-num">' + obj + '</span>';
  if (typeof obj === 'string') return '<span class="jt-str">"' + obj.replace(/</g,'&lt;') + '"</span>';
  if (Array.isArray(obj)) {{
    if (obj.length === 0) return '<span class="jt-bracket">[]</span>';
    if (obj.length <= 5 && obj.every(x => typeof x !== 'object' || x === null)) {{
      const items = obj.map(x => renderJson(x, depth+1)).join(', ');
      return '<span class="jt-bracket">[</span>' + items + '<span class="jt-bracket">]</span>';
    }}
    const id = 'jt' + Math.random().toString(36).slice(2,9);
    let html = '<span class="jt-toggle" onclick="toggleJt(\\'' + id + '\\')">&blacktriangledown;</span>';
    html += '<span class="jt-bracket">[</span> <span style="color:#94a3b8;font-size:.7rem;">' + obj.length + ' items</span>';
    html += '<div id="' + id + '" style="margin-left:1.2rem;">';
    obj.forEach((v, i) => {{ html += '<div>' + renderJson(v, depth+1) + (i < obj.length-1 ? ',' : '') + '</div>'; }});
    html += '</div><span class="jt-bracket">]</span>';
    return html;
  }}
  if (typeof obj === 'object') {{
    const keys = Object.keys(obj);
    if (keys.length === 0) return '<span class="jt-bracket">{{}}</span>';
    const id = 'jt' + Math.random().toString(36).slice(2,9);
    const collapsed = depth >= 2;
    let html = '<span class="jt-toggle" onclick="toggleJt(\\'' + id + '\\')">' +
               (collapsed ? '&blacktriangleright;' : '&blacktriangledown;') + '</span>';
    html += '<span class="jt-bracket">{{</span>';
    html += '<div id="' + id + '"' + (collapsed ? ' class="jt-collapsed"' : '') + ' style="margin-left:1.2rem;">';
    keys.forEach((k, i) => {{
      html += '<div><span class="jt-key">' + k + '</span>: ' +
              renderJson(obj[k], depth+1) + (i < keys.length-1 ? ',' : '') + '</div>';
    }});
    html += '</div><span class="jt-bracket">}}</span>';
    return html;
  }}
  return String(obj);
}}
function toggleJt(id) {{
  const el = document.getElementById(id);
  if (el.classList.contains('jt-collapsed')) {{
    el.classList.remove('jt-collapsed');
    const prev = el.previousElementSibling;
    if (prev && prev.previousElementSibling && prev.previousElementSibling.classList.contains('jt-toggle'))
      prev.previousElementSibling.innerHTML = '&blacktriangledown;';
  }} else {{
    el.classList.add('jt-collapsed');
    const prev = el.previousElementSibling;
    if (prev && prev.previousElementSibling && prev.previousElementSibling.classList.contains('jt-toggle'))
      prev.previousElementSibling.innerHTML = '&blacktriangleright;';
  }}
}}
// Render JSON trees
Object.keys(DOCS).forEach(sid => {{
  const el = document.getElementById('json-' + sid);
  if (el) el.innerHTML = renderJson(DOCS[sid], 0);
}});

// ─── Plotly Charts ───
const pLayout = {{
  paper_bgcolor:'#f8fafc', plot_bgcolor:'#f8fafc',
  font:{{ color:'#64748b', family:'-apple-system,sans-serif', size:11 }},
  margin:{{ l:55, r:15, t:40, b:45 }},
  xaxis:{{ gridcolor:'#e2e8f0', zerolinecolor:'#e2e8f0',
           title:{{ text:'Time (s)', font:{{ size:10 }} }} }},
  yaxis:{{ gridcolor:'#e2e8f0', zerolinecolor:'#e2e8f0' }},
}};
const pCfg = {{ responsive:true, displayModeBar:false }};

// Key monomers to highlight
const KEY_MONOMERS = [
  ['Free_fliG', 'FliG (26/switch)', '#1f77b4'],
  ['Free_fliM', 'FliM (34/switch)', '#ff7f0e'],
  ['Free_flgE', 'FlgE (120/hook)', '#2ca02c'],
  ['Free_fliH', 'FliH (12/export)', '#d62728'],
  ['Free_fliI', 'FliI (6/export)', '#9467bd'],
  ['Free_fliD', 'FliD (5/flagellum)', '#8c564b'],
];

const LOW_COPY = [
  ['Free_flhD', 'FlhD', '#1f77b4'],
  ['Free_flhC', 'FlhC', '#ff7f0e'],
  ['Free_fliN', 'FliN', '#2ca02c'],
  ['Free_flhA', 'FlhA', '#d62728'],
  ['Free_fliE', 'FliE', '#9467bd'],
  ['Free_motA', 'MotA', '#8c564b'],
  ['Free_fliD', 'FliD', '#e377c2'],
];

Object.keys(DATA).forEach(sid => {{
  const d = DATA[sid];

  // Chart 1: Complex assembly
  const complexTraces = [];
  Object.keys(d.complexes).forEach(obs => {{
    complexTraces.push({{
      x: d.times, y: d.complexes[obs], type:'scatter', mode:'lines',
      line:{{ color: COMPLEX_COLORS[obs] || '#94a3b8', width:2 }},
      name: COMPLEX_LABELS[obs] || obs,
    }});
  }});
  Plotly.newPlot('chart-complexes-'+sid, complexTraces, {{
    ...pLayout,
    title:{{ text:'Completed Complexes', font:{{ size:12, color:'#334155' }} }},
    yaxis:{{...pLayout.yaxis, title:{{ text:'Count', font:{{ size:10 }} }} }},
    legend:{{ font:{{ size:8 }}, bgcolor:'rgba(248,250,252,0.9)',
              bordercolor:'#e2e8f0', borderwidth:1 }},
    showlegend:true,
  }}, pCfg);

  // Chart 2: Flagella production detail
  const flagData = d.complexes['flagella'] || [];
  Plotly.newPlot('chart-flagella-'+sid, [{{
    x: d.times, y: flagData, type:'scatter', mode:'lines+markers',
    line:{{ color:'#e377c2', width:2.5 }}, marker:{{ size:4 }},
    fill:'tozeroy', fillcolor:'rgba(227,119,194,0.08)',
  }}], {{
    ...pLayout,
    title:{{ text:'Flagella Production', font:{{ size:12, color:'#334155' }} }},
    yaxis:{{...pLayout.yaxis, title:{{ text:'Complete Flagella', font:{{ size:10 }} }} }},
    showlegend:false,
  }}, pCfg);

  // Chart 3: Key monomers (high copy)
  const monoTraces = [];
  KEY_MONOMERS.forEach(([obs, label, color]) => {{
    if (d.monomers[obs]) {{
      monoTraces.push({{
        x: d.times, y: d.monomers[obs], type:'scatter', mode:'lines',
        line:{{ color, width:1.5 }}, name: label,
      }});
    }}
  }});
  Plotly.newPlot('chart-monomers-'+sid, monoTraces, {{
    ...pLayout,
    title:{{ text:'High-Copy Monomer Pools', font:{{ size:12, color:'#334155' }} }},
    yaxis:{{...pLayout.yaxis, title:{{ text:'Molecule Count', font:{{ size:10 }} }} }},
    legend:{{ font:{{ size:8 }}, bgcolor:'rgba(248,250,252,0.9)',
              bordercolor:'#e2e8f0', borderwidth:1 }},
    showlegend:true,
  }}, pCfg);

  // Chart 4: Low-copy monomers
  const lowTraces = [];
  LOW_COPY.forEach(([obs, label, color]) => {{
    if (d.monomers[obs]) {{
      lowTraces.push({{
        x: d.times, y: d.monomers[obs], type:'scatter', mode:'lines',
        line:{{ color, width:1.5 }}, name: label,
      }});
    }}
  }});
  Plotly.newPlot('chart-monomers2-'+sid, lowTraces, {{
    ...pLayout,
    title:{{ text:'Low-Copy Monomer Pools', font:{{ size:12, color:'#334155' }} }},
    yaxis:{{...pLayout.yaxis, title:{{ text:'Molecule Count', font:{{ size:10 }} }} }},
    legend:{{ font:{{ size:8 }}, bgcolor:'rgba(248,250,252,0.9)',
              bordercolor:'#e2e8f0', borderwidth:1 }},
    showlegend:true,
  }}, pCfg);
}});

</script>
</body>
</html>"""

    with open(output_path, 'w') as f:
        f.write(html)
    print(f'Report saved to {output_path}')


# ── Main ───────────────────────────────────────────────────────────

def run_demo():
    import subprocess

    model_path = get_model_path()
    obs_names, monomer_obs, complex_obs, growing_obs = _classify_observables(model_path)

    demo_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(demo_dir, 'report.html')

    sim_results = []
    for cfg in CONFIGS:
        print(f'Running: {cfg["title"]}...')

        if cfg['mode'] == 'standalone':
            tc, runtime = run_standalone(
                cfg, model_path, obs_names, monomer_obs, complex_obs, growing_obs)
        else:
            tc, runtime = run_composed(
                cfg, model_path, obs_names, monomer_obs, complex_obs, growing_obs)

        sim_results.append((cfg, tc, runtime))
        print(f'  Runtime: {runtime:.2f}s')
        print(f'  {len(tc["time"])} time points')
        flagella_count = tc['complexes'].get('flagella', [0])[-1] if tc['complexes'].get('flagella') else 0
        print(f'  Final flagella: {flagella_count:.0f}')

    print('Generating HTML report...')
    generate_html(sim_results, output_path)

    # Auto-open in Safari
    subprocess.run(['open', '-a', 'Safari', output_path])
    print('Report opened in Safari.')


if __name__ == '__main__':
    run_demo()
