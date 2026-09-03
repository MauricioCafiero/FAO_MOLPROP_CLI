#!/usr/bin/env python3
"""
run_ga_baseline.py - Non-LLM baseline: a small population-based genetic
algorithm over the same base-ring vocabulary already defined in MolPropOp.py
(the vocabulary code/adversarial_set.md was combinatorially enumerated from),
paired with one of two substituent menus selected via --pool (see below).

Purpose: every other baseline in this study (zero/few/frag-shot, the agentic
loop) is an LLM proposing molecules. This one isn't -- it's a classic
evolutionary search using ONLY real GFN2-xTB HOMO-LUMO gap as fitness, no
chemical "reasoning" at all. It answers "how much of the agentic loop's
advantage is really about search-with-feedback (dock, keep what's good,
discard what's bad, repeat), versus something only an LLM proposer supplies?"
A GA that matches or beats the agentic loop's HL gaps at a comparable
evaluation budget would say the loop's gain is mostly generic feedback-driven
search; a GA that falls well short would say the LLM's chemical judgment is
doing real work beyond just having a scored feedback loop.

--pool controls which chemical space the GA can search, so the comparison
stays apples-to-apples with whichever condition it's being set against:
'frag10' (default) restricts it to the exact 10 fragments frag-shot's prompt
showed the LLM, so a GA-vs-frag-shot gap reflects search strategy, not a
bigger available chemical space. 'full' opens up the entire combined
e_withdraw/e_donate(+linker) pool (~390 items), for comparison against
zero-shot, which was never shown a fragment menu at all -- there, frag10
would be the artificial handicap instead.

Genome: (ring_index, {position: substituent_smiles}) -- one of the 9
MolPropOp.base_rings plus 0..min(3, len(clean positions)) substituents drawn
from the selected --pool, placed at MolPropOp's pre-defined "clean"
(symmetrically unique) ring positions. Built into a SMILES string the same
way MolPropOp.sub_cycle does (splice substituents into the ring SMILES at
each clean position, highest position index first so earlier indices stay
valid), then validated with RDKit.

Fitness: the real GFN2-xTB HOMO-LUMO gap via hl_gap_module.scoring_function
-- the identical scorer every other HL batch uses -- with results cached by
canonical SMILES so a genotype revisited by mutation/crossover (or kept via
elitism) is never recomputed. Lower is better, as in docking; the 100.0
failure sentinel makes an unscorable genotype lose every tournament. There is
no pocket here, so the docking version's in-pocket columns are dropped.
QED/aLogP come from RDKit's own
QED.properties (not parsed from any model text, since there is no model here);
SAS/NP come from sascorer/npscorer directly.

Budget: --pop x --gens is the total unique-genotype evaluation budget per
replicate (duplicates from caching/elitism are free). --preset 5x4 sets
pop=5, gens=4 (20-eval budget, matching the agentic study's 5x4 naming and
its typical ~5-compounds-delivered-per-replicate scale); --preset 10x8 sets
pop=8, gens=10 (80-eval budget). Each replicate's *final population* is
written out as that replicate's "delivered compounds" -- the GA's analogue of
a proposer's final turn -- in the same compounds_<batch>.csv column schema
analyze_replicates.py produces, so this drops straight into the same
downstream analysis (analyze_shot_vs_agentic_stats.py, SUMMARY_TABLES.md
conventions) as every other condition in this study.

Checkpointing: every dock result is appended to the output CSV immediately
(not buffered to the end). Gap evaluation is far cheaper than docking
(~0.02-0.1 s/molecule vs minutes), so a full GA run is minutes not hours, but
the checkpointing is kept so a reaped process still loses nothing.

Usage:
  fao-env/bin/python code/run_ga_baseline.py --preset 5x4 --replicates 5                     # vs. frag-shot (default pool=frag10)
  fao-env/bin/python code/run_ga_baseline.py --preset 10x8 --replicates 3
  fao-env/bin/python code/run_ga_baseline.py --preset 5x4 --replicates 5 --pool full          # vs. zero-shot
  fao-env/bin/python code/run_ga_baseline.py --preset 5x4 --replicates 1 --skip-docking       # smoke test, no Vina calls
"""
import argparse
import csv
import os
import random
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root

from rdkit import Chem
from rdkit.Chem import QED, RDConfig

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
sys.path.append(os.path.join(RDConfig.RDContribDir, 'NP_Score'))
import sascorer  # noqa: E402
import npscorer  # noqa: E402

from hl_gap_module import scoring_function, scoring_args  # noqa: E402
from MolPropOp import (  # noqa: E402
    base_rings, clean_ring_locations, e_withdraw, e_donate,
    withdraw_with_linkers, donate_with_linkers,
)

# Two substituent menus, selected via --pool:
#
# 'frag10' (default) -- the exact 10-item functional-group menu shown to the
# LLM in run_zero_few_shot.py's FRAG_SHOT_SYSTEM prompt (itself a fixed sample
# drawn from the 'full' pool below). Use this to compare against frag-shot:
# holding rings AND fragments identical isolates search-with-feedback vs.
# static LLM reasoning as the only difference. Using the full pool there would
# let a GA-vs-frag-shot score gap just reflect a larger accessible substituent
# space instead.
#
# 'full' -- the entire combined e_withdraw/e_donate(+linker) pool (~390
# items). Use this to compare against zero-shot, which was never shown any
# fragment menu at all -- there, constraining the GA to 10 items would be the
# artificial handicap.
# NOTE: these are the exact 10 in code_hl/adversarial_set.md and in frag-shot's HL
# prompt -- all linker+group combinations. They differ from the docking menu, which
# uses bare 'I' and 'C#N' where the HL seed uses 'N(I)' (-NH- + iodo) and 'O(C#N)'
# (-O- + nitrile). Matching the HL seed is what keeps GA-vs-frag-shot pool-matched.
FRAG10_SUBSTITUENTS = [
    'N(I)',
    'O(C#N)',
    'C(=O)O(C(C)C)',
    'C#C(SC)',
    'C(C(=O)[O-])',
    'C(C)',
    'C=C([N+](=O)[O-])',
    'C(N)',
    'C([O-])',
    'CC(N(C)C)',
]
POOLS = {
    'frag10': FRAG10_SUBSTITUENTS,
    'full': e_withdraw + e_donate + withdraw_with_linkers + donate_with_linkers,
}
SUBSTITUENT_POOL = POOLS['frag10']  # overwritten in main() per --pool
_FAIL_SENTINEL = 100.0  # hl_gap_module's failure value
MAX_SUBS = 3  # cap substituents per genome, regardless of how many clean positions a ring has

PRESETS = {
    '5x4':  {'pop': 5, 'gens': 4},
    '10x8': {'pop': 8, 'gens': 10},
}


# --- genome <-> SMILES -------------------------------------------------------

def build_smiles(ring_idx, subs):
    """subs: dict {position: substituent_smiles}. Splice highest position first
    so earlier indices (defined against the original ring string) stay valid."""
    s = base_rings[ring_idx]
    for loc in sorted(subs.keys(), reverse=True):
        e = subs[loc]
        s = f'{s[:loc + 1]}({e}){s[loc + 1:]}'
    return s


def _validate(ring_idx, subs):
    smi = build_smiles(ring_idx, subs)
    mol = Chem.MolFromSmiles(smi)
    return (smi, mol) if mol is not None else (None, None)


def random_genome(rng):
    """Retry until a valid (ring, substituent-set) combination is found."""
    for _ in range(200):
        ring_idx = rng.randrange(len(base_rings))
        locs = clean_ring_locations[ring_idx]
        k = rng.randint(0, min(MAX_SUBS, len(locs)))
        chosen_locs = rng.sample(locs, k) if k else []
        subs = {loc: rng.choice(SUBSTITUENT_POOL) for loc in chosen_locs}
        smi, mol = _validate(ring_idx, subs)
        if smi is not None:
            return {'ring_idx': ring_idx, 'subs': subs, 'smiles': smi}
    raise RuntimeError("could not build a valid random genome after 200 tries")


def mutate(genome, rng):
    ring_idx, subs = genome['ring_idx'], dict(genome['subs'])
    op = rng.choice(['swap_ring', 'add_sub', 'remove_sub', 'replace_sub'])
    for _ in range(50):
        new_ring_idx, new_subs = ring_idx, dict(subs)
        if op == 'swap_ring':
            new_ring_idx = rng.randrange(len(base_rings))
            valid_locs = set(clean_ring_locations[new_ring_idx])
            new_subs = {loc: e for loc, e in subs.items() if loc in valid_locs}
            if not new_subs and valid_locs:
                loc = rng.choice(list(valid_locs))
                new_subs = {loc: rng.choice(SUBSTITUENT_POOL)}
        elif op == 'add_sub':
            free = [l for l in clean_ring_locations[ring_idx] if l not in subs]
            if free and len(subs) < MAX_SUBS:
                new_subs[rng.choice(free)] = rng.choice(SUBSTITUENT_POOL)
            else:
                op = 'replace_sub'
                continue
        elif op == 'remove_sub':
            if subs:
                del new_subs[rng.choice(list(subs.keys()))]
            else:
                op = 'add_sub'
                continue
        elif op == 'replace_sub':
            if subs:
                loc = rng.choice(list(subs.keys()))
                new_subs[loc] = rng.choice(SUBSTITUENT_POOL)
            else:
                op = 'add_sub'
                continue
        smi, mol = _validate(new_ring_idx, new_subs)
        if smi is not None:
            return {'ring_idx': new_ring_idx, 'subs': new_subs, 'smiles': smi}
    return dict(genome)  # gave up mutating cleanly; keep parent unchanged


def crossover(parent_a, parent_b, rng):
    if parent_a['ring_idx'] == parent_b['ring_idx']:
        ring_idx = parent_a['ring_idx']
        subs = {}
        for loc in clean_ring_locations[ring_idx]:
            src = None
            if loc in parent_a['subs'] and loc in parent_b['subs']:
                src = rng.choice([parent_a, parent_b])
            elif loc in parent_a['subs']:
                src = parent_a if rng.random() < 0.5 else None
            elif loc in parent_b['subs']:
                src = parent_b if rng.random() < 0.5 else None
            if src is not None:
                subs[loc] = src['subs'][loc]
        if len(subs) > MAX_SUBS:
            keep = rng.sample(list(subs.keys()), MAX_SUBS)
            subs = {k: subs[k] for k in keep}
        smi, mol = _validate(ring_idx, subs)
        if smi is not None:
            return {'ring_idx': ring_idx, 'subs': subs, 'smiles': smi}
    # different rings (or same-ring crossover produced nothing valid):
    # clone one parent's ring+subs, then try folding in one substituent value
    # from the other parent at a random empty position.
    base, other = rng.sample([parent_a, parent_b], 2)
    ring_idx, subs = base['ring_idx'], dict(base['subs'])
    free = [l for l in clean_ring_locations[ring_idx] if l not in subs]
    if free and other['subs'] and len(subs) < MAX_SUBS:
        subs[rng.choice(free)] = rng.choice(list(other['subs'].values()))
        smi, mol = _validate(ring_idx, subs)
        if smi is not None:
            return {'ring_idx': ring_idx, 'subs': subs, 'smiles': smi}
    return dict(base)


# --- fitness (real docking, cached by canonical SMILES) ---------------------

class Evaluator:
    """Dock-and-score genomes, caching by canonical SMILES so a genotype
    revisited across generations (elitism, convergence, crossover) is never
    re-docked. Does NOT write to the output CSV itself -- only the caller
    decides which evaluated genotypes are worth persisting (the final
    population, see run_replicate/main), matching the agentic study's
    compounds_<batch>.csv convention of logging final delivered compounds
    only, not every intermediate turn."""

    def __init__(self, method, skip_gaps, np_model=None, log_writer=None, log_meta=None):
        self.method = method
        self.skip_gaps = skip_gaps
        self.cache = {}  # canonical_smiles -> result dict
        self.n_gap_calls = 0
        self._np_model = np_model
        self._log_writer = log_writer  # every evaluated genotype, checkpointed immediately
        self._log_meta = log_meta or {}

    def evaluate(self, genome, generation):
        mol = Chem.MolFromSmiles(genome['smiles'])
        canon = Chem.MolToSmiles(mol)
        if canon in self.cache:
            return self.cache[canon]

        gap = None
        if not self.skip_gaps:
            scoring_args[0] = os.cpu_count()
            scoring_args[1] = self.method
            try:
                score, _aux = scoring_function(genome['smiles'])
                self.n_gap_calls += 1
                # hl_gap_module returns 100.0 on failure (bad SMILES, embed failure,
                # non-convergence); record None rather than a fake enormous gap.
                gap = None if (score is None or score >= _FAIL_SENTINEL) else float(score)
            except Exception:
                gap = None

        try:
            qed = QED.default(mol)
            alogp = QED.properties(mol)[1]
        except Exception:
            qed = alogp = None
        try:
            sas = sascorer.calculateScore(mol)
            np_score = npscorer.scoreMol(mol, self._np_model) if self._np_model else None
        except Exception:
            sas = np_score = None

        # unscored genotype: the 100.0 failure sentinel, so it loses every
        # tournament against any real gap (all real gaps are single-digit eV).
        fitness = gap if gap is not None else _FAIL_SENTINEL
        result = {
            'canonical_smiles': canon, 'original_smiles': genome['smiles'],
            'inchikey': Chem.MolToInchiKey(mol), 'gap': gap,
            'qed': qed, 'alogp': alogp, 'sas': sas, 'np': np_score,
            'fitness': fitness,
        }
        result['generation'] = generation
        self.cache[canon] = result
        if self._log_writer is not None:
            row = dict(self._log_meta)
            row.update({k: v for k, v in result.items() if k != 'fitness'})
            row['method'] = self.method
            self._log_writer.writerow(row)
        return result


# --- GA loop ------------------------------------------------------------------

def run_replicate(rng, pop_size, gens, evaluator, tournament_k=2, elitism=1, mutation_rate=0.7):
    population = [random_genome(rng) for _ in range(pop_size)]
    fitness = {id(g): evaluator.evaluate(g, generation=0) for g in population}

    for gen in range(1, gens + 1):
        ranked = sorted(population, key=lambda g: fitness[id(g)]['fitness'])
        next_pop = ranked[:elitism]  # elites carry over unchanged (already evaluated)

        def tournament():
            contenders = rng.sample(population, min(tournament_k, len(population)))
            return min(contenders, key=lambda g: fitness[id(g)]['fitness'])

        while len(next_pop) < pop_size:
            p1, p2 = tournament(), tournament()
            child = crossover(p1, p2, rng)
            if rng.random() < mutation_rate:
                child = mutate(child, rng)
            next_pop.append(child)

        population = next_pop
        fitness = {id(g): evaluator.evaluate(g, generation=gen) for g in population}

    ranked = sorted(population, key=lambda g: fitness[id(g)]['fitness'])
    return ranked, fitness


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--preset', choices=list(PRESETS), default='5x4')
    p.add_argument('--pop', type=int, default=None, help='override preset population size')
    p.add_argument('--gens', type=int, default=None, help='override preset generation count')
    p.add_argument('--replicates', type=int, default=5)
    p.add_argument('--method', default='GFN2-xTB')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--mutation-rate', type=float, default=0.7)
    p.add_argument('--tournament-k', type=int, default=2)
    p.add_argument('--elitism', type=int, default=1)
    p.add_argument('--pool', choices=list(POOLS), default='frag10',
                    help="'frag10' (default): the exact 10-fragment menu shown to the LLM in "
                         "frag-shot -- use for the GA-vs-frag-shot comparison. 'full': the entire "
                         "~390-item combined substituent pool -- use for the GA-vs-zero-shot "
                         "comparison, where no fragment menu was shown to the LLM either.")
    p.add_argument('--skip-gaps', action='store_true',
                    help='CPU-light wiring check: build/mutate genomes but skip the GFN2-xTB call')
    p.add_argument('--out-dir', default=None,
                    help='default: results/batches/hl_batches/ga_baseline/<preset>[_<pool>]/analysis')
    args = p.parse_args(argv)

    global SUBSTITUENT_POOL
    SUBSTITUENT_POOL = POOLS[args.pool]

    pop = args.pop or PRESETS[args.preset]['pop']
    gens = args.gens or PRESETS[args.preset]['gens']
    # frag10 keeps the original ga_<preset> naming (existing analysis scripts
    # reference this path directly); full gets a distinct suffix so it never
    # collides with or overwrites the frag10 run.
    preset_dir = args.preset if args.pool == 'frag10' else f'{args.preset}_{args.pool}'
    set_label = f'ga_{preset_dir}'

    out_dir = args.out_dir or os.path.join(_ROOT, 'results', 'batches', 'hl_batches',
                                          'ga_baseline', preset_dir, 'analysis')
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, f'compounds_{set_label}.csv')
    log_csv = os.path.join(out_dir, f'ga_eval_log_{set_label}.csv')

    fieldnames = ['set_label', 'replicate', 'proposer_provider', 'proposer_model',
                  'adversary_provider', 'adversary_model', 'method',
                  'original_smiles', 'canonical_smiles', 'inchikey', 'gap',
                  'qed', 'alogp', 'sas', 'np',
                  'source_turn']
    log_fieldnames = fieldnames + ['generation']

    np_model = npscorer.readNPModel()  # cheap (~1s); loaded regardless of --skip-docking

    write_out_header = not os.path.exists(out_csv)
    out_f = open(out_csv, 'a', newline='')
    out_writer = csv.DictWriter(out_f, fieldnames=fieldnames)
    if write_out_header:
        out_writer.writeheader()
    out_f.flush()

    # Checkpointed as each genotype is docked -- the CPU-heavy step -- so a
    # killed/reaped run loses at most the in-flight dock call, not the whole
    # replicate. compounds_<set_label>.csv (above) is derived from the final
    # population only and written once per completed replicate.
    write_log_header = not os.path.exists(log_csv)
    log_f = open(log_csv, 'a', newline='')
    log_writer = csv.DictWriter(log_f, fieldnames=log_fieldnames)
    if write_log_header:
        log_writer.writeheader()
    log_f.flush()

    print(f"GA baseline: preset={args.preset} pool={args.pool} ({len(SUBSTITUENT_POOL)} substituents) "
          f"pop={pop} gens={gens} replicates={args.replicates} "
          f"method={args.method} skip_gaps={args.skip_gaps}")
    print(f"Delivered-compounds CSV: {out_csv}")
    print(f"Full evaluation log (checkpointed per gap call): {log_csv}")

    t0 = time.time()
    for rep in range(1, args.replicates + 1):
        rng = random.Random(args.seed + rep)
        meta = {
            'set_label': set_label, 'replicate': rep,
            'proposer_provider': None, 'proposer_model': 'genetic_algorithm',
            'adversary_provider': None, 'adversary_model': None,
        }
        evaluator = Evaluator(args.method, args.skip_gaps, np_model=np_model,
                               log_writer=log_writer, log_meta=meta)
        ranked, fitness = run_replicate(rng, pop, gens, evaluator,
                                          tournament_k=args.tournament_k,
                                          elitism=args.elitism,
                                          mutation_rate=args.mutation_rate)
        log_f.flush()

        # final population = this replicate's "delivered compounds" (dedup by
        # canonical SMILES, same way a proposer's turn can repeat a molecule).
        seen = set()
        for g in ranked:
            r = fitness[id(g)]
            if r['canonical_smiles'] in seen:
                continue
            seen.add(r['canonical_smiles'])
            row = dict(meta)
            row.update({k: v for k, v in r.items() if k not in ('fitness', 'generation')})
            row['method'] = args.method
            row['source_turn'] = gens
            out_writer.writerow(row)
        out_f.flush()

        best = fitness[id(ranked[0])]
        print(f"  rep {rep}/{args.replicates}: {evaluator.n_gap_calls} gap calls, "
              f"{len(seen)} delivered compounds, best={best['gap']} "
              f"({best['canonical_smiles']}) [{time.time() - t0:.0f}s elapsed]")

    out_f.close()
    log_f.close()
    print(f"Done. {time.time() - t0:.0f}s total.")


if __name__ == '__main__':
    main()
