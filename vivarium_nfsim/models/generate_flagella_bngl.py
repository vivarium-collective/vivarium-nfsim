"""
Generate a BioNetGen model for E. coli flagella complexation.

Reads the complexation stoichiometry from vivarium-chemotaxis and generates
a BNGL file with sequential bimolecular binding rules for each assembly step.

Each multi-subunit complex is modeled as a scaffold molecule with counter
states that track how many of each subunit have been incorporated. Monomers
bind one at a time via bimolecular reactions. When all subunits are bound,
a final rule converts the scaffold into the completed complex.
"""
import os

# ---------------------------------------------------------------------------
# Stoichiometry from vivarium-chemotaxis flagella_chromosome.py
# Negative values = consumed, positive = produced
# ---------------------------------------------------------------------------
COMPLEXATION_STOICHIOMETRY = {
    'flhDC': {
        'flhD': -4.0,
        'flhC': -2.0,
        'flhDC': 1.0,
    },
    'flagellar motor switch reaction': {
        'flagellar motor switch': 1.0,
        'fliG': -26.0,
        'fliM': -34.0,
        'fliN': -1.0,
    },
    'flagellar export apparatus reaction 1': {
        'flagellar export apparatus subunit': 1.0,
        'flhA': -1.0,
        'flhB': -1.0,
        'fliO': -1.0,
        'fliP': -1.0,
        'fliQ': -1.0,
        'fliR': -1.0,
        'fliJ': -1.0,
        'fliI': -6.0,
    },
    'flagellar export apparatus reaction 2': {
        'flagellar export apparatus': 1.0,
        'flagellar export apparatus subunit': -1.0,
        'fliH': -12.0,
    },
    'flagellar motor reaction': {
        'flagellar motor': 1.0,
        'flagellar motor switch': -1.0,
        'fliL': -2.0,
        'flgH': -1.0,
        'motA': -1.0,
        'motB': -1.0,
        'flgB': -1.0,
        'flgC': -1.0,
        'flgF': -1.0,
        'flgG': -1.0,
        'flgI': -1.0,
        'fliF': -1.0,
        'fliE': -1.0,
    },
    'flagellar hook reaction': {
        'flagellar hook': 1,
        'flgE': -120.0,
    },
    'flagellum reaction': {
        'flagella': 1.0,
        'flagellar export apparatus': -1.0,
        'flagellar motor': -1.0,
        'fliC': -1.0,
        'flgL': -1.0,
        'flgK': -1.0,
        'fliD': -5.0,
        'flagellar hook': -1,
    },
}

# Rate constants
# Note: vivarium-chemotaxis uses k=1e-4 in a deterministic ODE context.
# For NFSim stochastic simulation with small molecule counts, we need
# higher rates to get assembly on a reasonable timescale.
K_BIND = 5e-1  # bimolecular binding rate for growth steps (1/molecule/s)
K_NUCLEATION = 5e-2  # slower nucleation rate to limit intermediate count
K_COMPLETION = 10.0  # unimolecular completion rate (1/s)

# Number of flagella worth of monomers to provide
N_FLAGELLA = 5


def _safe_name(name):
    """Convert a name to a valid BNG identifier."""
    return name.replace(' ', '_').replace('-', '_')


def _parse_reaction(rxn_name, stoich):
    """Parse a reaction into consumed monomers and produced complex."""
    consumed = {}
    product = None
    for species, count in stoich.items():
        if count < 0:
            consumed[species] = int(abs(count))
        elif count > 0:
            product = species
    return consumed, product


def generate_bngl(n_flagella=N_FLAGELLA, k_bind=K_BIND, k_nucleation=K_NUCLEATION, k_completion=K_COMPLETION):
    """Generate the complete BNGL model string."""

    # Ordered reactions (assembly hierarchy)
    reaction_order = [
        'flhDC',
        'flagellar motor switch reaction',
        'flagellar export apparatus reaction 1',
        'flagellar export apparatus reaction 2',
        'flagellar motor reaction',
        'flagellar hook reaction',
        'flagellum reaction',
    ]

    # Parse all reactions
    reactions = {}
    for rxn_name in reaction_order:
        consumed, product = _parse_reaction(rxn_name, COMPLEXATION_STOICHIOMETRY[rxn_name])
        reactions[rxn_name] = {
            'consumed': consumed,
            'product': product,
        }

    # Collect all monomer species (those that are never products of a reaction)
    complex_names = set()
    for rxn in reactions.values():
        complex_names.add(rxn['product'])

    all_consumed = set()
    for rxn in reactions.values():
        all_consumed.update(rxn['consumed'].keys())

    monomer_names = sorted(all_consumed - complex_names)
    complex_names_ordered = [reactions[r]['product'] for r in reaction_order]

    # Calculate initial monomer counts (enough for n_flagella flagella)
    monomer_counts = {}
    for rxn in reactions.values():
        for species, count in rxn['consumed'].items():
            if species in monomer_names:
                needed = count * n_flagella
                monomer_counts[species] = max(monomer_counts.get(species, 0), needed)

    # For sub-complexes consumed by later reactions, they start at 0
    # (they'll be produced during simulation)

    # ---- Build BNGL ----
    lines = []
    lines.append('begin model')
    lines.append('')

    # -- Parameters --
    lines.append('begin parameters')
    lines.append(f'    n_flagella  {n_flagella}')
    lines.append(f'    k_bind      {k_bind}')
    lines.append(f'    k_nucleation {k_nucleation}')
    lines.append(f'    k_completion {k_completion}')
    lines.append('')
    # Per-complex nucleation rates: computed to yield ~n_flagella nucleation
    # events in the first ~50s. Rate is set based on initial monomer counts
    # so that propensity = target_propensity.
    # For A+A: propensity = k * N*(N-1)/2
    # For A+B: propensity = k * NA * NB
    target_propensity = n_flagella / 50.0  # nucleate over ~50 seconds
    for rxn_name in reaction_order:
        rxn = reactions[rxn_name]
        consumed = rxn['consumed']
        product = rxn['product']
        total = sum(consumed.values())
        if total > 2:
            safe_product = _safe_name(product)
            sorted_sp = sorted(consumed.keys(), key=lambda s: consumed[s])
            nuc1 = sorted_sp[0]
            if consumed[nuc1] >= 2:
                nuc2 = nuc1
            else:
                nuc2 = sorted_sp[1]

            # Get initial counts
            n1 = monomer_counts.get(nuc1, n_flagella)
            n2 = monomer_counts.get(nuc2, n_flagella)

            if nuc1 == nuc2:
                combinatorial = n1 * (n1 - 1) / 2.0
            else:
                combinatorial = n1 * n2

            if combinatorial > 0:
                nuc_rate = target_propensity / combinatorial
            else:
                nuc_rate = k_nucleation

            lines.append(f'    k_nuc_{safe_product}  {nuc_rate:.6e}')
    lines.append('')

    # Initial counts for each monomer
    for monomer in sorted(monomer_counts.keys()):
        safe = _safe_name(monomer)
        lines.append(f'    {safe}_0  {monomer_counts[monomer]}')
    lines.append('end parameters')
    lines.append('')

    # -- Molecule Types --
    lines.append('begin molecule types')

    # Monomers - simple molecules
    for monomer in monomer_names:
        lines.append(f'    {_safe_name(monomer)}()')

    # Sub-complexes that are both products and reactants
    intermediate_complexes = sorted(complex_names & all_consumed)
    for cx in intermediate_complexes:
        lines.append(f'    {_safe_name(cx)}()')

    # For each reaction, create a "Growing" scaffold with counter states
    for rxn_name in reaction_order:
        rxn = reactions[rxn_name]
        consumed = rxn['consumed']
        product = rxn['product']
        safe_product = _safe_name(product)

        total_subunits = sum(consumed.values())

        if total_subunits <= 2:
            # Simple bimolecular - no scaffold needed
            continue

        # Scaffold molecule with one counter state per subunit type
        state_parts = []
        for species in sorted(consumed.keys()):
            count = consumed[species]
            safe_species = _safe_name(species)
            states = '~'.join(str(i) for i in range(count + 1))
            state_parts.append(f'{safe_species}~{states}')

        scaffold_name = f'Growing_{safe_product}'
        lines.append(f'    {scaffold_name}({",".join(state_parts)})')

    # Final complexes (products that are never consumed)
    final_complexes = sorted(complex_names - all_consumed)
    for cx in final_complexes:
        lines.append(f'    {_safe_name(cx)}()')

    lines.append('end molecule types')
    lines.append('')

    # -- Seed Species --
    lines.append('begin seed species')
    for monomer in sorted(monomer_counts.keys()):
        safe = _safe_name(monomer)
        lines.append(f'    {safe}()  {safe}_0')
    lines.append('end seed species')
    lines.append('')

    # -- Observables --
    lines.append('begin observables')

    # Free monomers
    for monomer in monomer_names:
        safe = _safe_name(monomer)
        lines.append(f'    Molecules  Free_{safe}  {safe}()')

    # Completed complexes
    for cx_name in complex_names_ordered:
        safe = _safe_name(cx_name)
        lines.append(f'    Molecules  {safe}  {safe}()')

    # Growing intermediates (total count)
    for rxn_name in reaction_order:
        rxn = reactions[rxn_name]
        consumed = rxn['consumed']
        product = rxn['product']
        safe_product = _safe_name(product)
        total_subunits = sum(consumed.values())
        if total_subunits > 2:
            scaffold_name = f'Growing_{safe_product}'
            lines.append(f'    Molecules  {scaffold_name}_total  {scaffold_name}()')

    lines.append('end observables')
    lines.append('')

    # -- Reaction Rules --
    lines.append('begin reaction rules')

    for rxn_name in reaction_order:
        rxn = reactions[rxn_name]
        consumed = rxn['consumed']
        product = rxn['product']
        safe_product = _safe_name(product)
        total_subunits = sum(consumed.values())

        lines.append(f'')
        lines.append(f'    # === {rxn_name} ===')
        lines.append(f'    # Product: {product}')
        lines.append(f'    # Subunits: {", ".join(f"{c}x {s}" for s, c in consumed.items())}')

        if total_subunits == 1:
            # Unimolecular conversion (e.g., a single subcomplex -> new complex)
            species = list(consumed.keys())[0]
            safe_species = _safe_name(species)
            lines.append(f'    {safe_species}() -> {safe_product}()  k_bind')

        elif total_subunits == 2:
            # Simple bimolecular
            species_list = []
            for species, count in consumed.items():
                for _ in range(count):
                    species_list.append(species)

            if len(species_list) == 2:
                s1, s2 = species_list
                lines.append(f'    {_safe_name(s1)}() + {_safe_name(s2)}() -> {safe_product}()  k_bind')

        else:
            # Multi-subunit: use Growing scaffold
            scaffold_name = f'Growing_{safe_product}'
            sorted_species = sorted(consumed.keys())

            # Pick nucleation pair: use the SCARCEST monomer to avoid
            # creating too many intermediates that can never complete.
            # Among species with count=1 (scarcest), pick first two distinct.
            # If only one species, use the one with lowest required count.
            species_by_count = sorted(consumed.keys(), key=lambda s: consumed[s])
            nuc_species_1 = species_by_count[0]
            if consumed[nuc_species_1] >= 2:
                # Can self-nucleate with the scarcest
                nuc_species_2 = nuc_species_1
            else:
                # Use two different scarce species
                nuc_species_2 = species_by_count[1]

            # Initial state after nucleation
            init_states = []
            nuc_counts = {}
            for species in sorted_species:
                if species == nuc_species_1:
                    nuc_counts[species] = nuc_counts.get(species, 0) + 1
                if species == nuc_species_2:
                    nuc_counts[species] = nuc_counts.get(species, 0) + 1

            for species in sorted_species:
                safe_sp = _safe_name(species)
                c = nuc_counts.get(species, 0)
                init_states.append(f'{safe_sp}~{c}')

            safe_nuc1 = _safe_name(nuc_species_1)
            safe_nuc2 = _safe_name(nuc_species_2)

            # Scale nucleation rate inversely with complex size to prevent
            # over-nucleation of large complexes (e.g., 120-subunit hook)
            nuc_rate_name = f'k_nuc_{safe_product}'
            lines.append(f'    # Nucleation (rate scaled by 1/total_subunits)')
            lines.append(f'    {safe_nuc1}() + {safe_nuc2}() -> '
                         f'{scaffold_name}({",".join(init_states)})  {nuc_rate_name}')

            # Sequential binding rules for each species
            for species in sorted_species:
                safe_sp = _safe_name(species)
                count = consumed[species]

                # Start from count already incorporated during nucleation
                start = nuc_counts.get(species, 0)

                for i in range(start, count):
                    lines.append(
                        f'    {scaffold_name}({safe_sp}~{i}) + {_safe_name(species)}() -> '
                        f'{scaffold_name}({safe_sp}~{i + 1})  k_bind')

            # Completion: when all counters are at max, convert to product
            complete_states = []
            for species in sorted_species:
                safe_sp = _safe_name(species)
                complete_states.append(f'{safe_sp}~{consumed[species]}')

            lines.append(f'    # Completion')
            lines.append(f'    {scaffold_name}({",".join(complete_states)}) -> '
                         f'{safe_product}()  k_completion')

    lines.append('')
    lines.append('end reaction rules')
    lines.append('')
    lines.append('end model')

    return '\n'.join(lines)


def write_bngl(output_path=None, **kwargs):
    """Generate and write the BNGL model file."""
    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(__file__), 'flagella_complexation.bngl')

    bngl_text = generate_bngl(**kwargs)

    with open(output_path, 'w') as f:
        f.write(bngl_text)

    print(f'Wrote BNGL model to {output_path}')
    return output_path


if __name__ == '__main__':
    path = write_bngl()
    print(f'\nModel written to: {path}')

    # Print summary
    with open(path) as f:
        text = f.read()
    n_rules = text.count(' -> ')
    n_mol_types = text.count('    ') // 2  # rough
    print(f'Reaction rules: {n_rules}')
