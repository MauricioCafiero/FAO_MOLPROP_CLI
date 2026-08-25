#!/usr/bin/env python3
"""
run_ga_baseline.py - Non-LLM baseline: a small population-based genetic
algorithm over the same fragment vocabulary (base rings + electron-
withdrawing/donating substituents, with/without linkers) already defined in
MolPropOp.py -- the vocabulary code/adversarial_set.md was combinatorially
enumerated from and the frag-shot baseline (ZERO_FEW_SHOT_BASELINE.md) draws
its ring/functional-group menu from.

Purpose: every other baseline in this study (zero/few/frag-shot, the agentic
loop) is an LLM proposing molecules. This one isn't -- it's a classic
evolutionary search using ONLY real Vina docking score as fitness, no
chemical "reasoning" at all. It answers "how much of the agentic loop's
advantage is really about search-with-feedback (dock, keep what's good,
discard what's bad, repeat), versus something only an LLM proposer supplies?"
A GA that matches or beats the agentic loop's docking scores at a comparable
evaluation budget would say the loop's gain is mostly generic feedback-driven
search; a GA that falls well short would say the LLM's chemical judgment is
doing real work beyond just having a scored feedback loop.

Genome: (ring_index, {position: substituent_smiles}) -- one of the 9
MolPropOp.base_rings plus 0..min(3, len(clean positions)) substituents drawn
from the combined e_withdraw/e_donate(+linker) pool, placed at MolPropOp's
pre-defined "clean" (symmetrically unique) ring positions. Built into a SMILES
string the same way MolPropOp.sub_cycle does (splice substituents into the
ring SMILES at each clean position, highest position index first so earlier
indices stay valid), then validated with RDKit.

Fitness: real dockstring/Vina docking score via docking_module.scoring_function
-- the identical scorer every other batch in this study uses -- with results
cached by canonical SMILES so a genotype revisited by mutation/crossover
(or kept via elitism) is never re-docked. In-pocket status uses
docking_module.contacted_residues/target_residues, matching
analyze_replicates.py's definition exactly. QED/aLogP come from RDKit's own
QED.properties (not parsed from any model text, since there is no model here);
SAS/NP come from docking_module.calculate_SAS_and_NP's underlying scorers.

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
(not buffered to the end), since docking is the CPU-heavy, potentially
multi-hour step in this pipeline and a killed/reaped process should not lose
completed work (see project memory on this exact failure mode).

Usage:
  fao-env/bin/python code/run_ga_baseline.py --preset 5x4 --replicates 5
  fao-env/bin/python code/run_ga_baseline.py --preset 10x8 --replicates 3
  fao-env/bin/python code/run_ga_baseline.py --preset 5x4 --replicates 1 --skip-docking   # smoke test, no Vina calls
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

from docking_module import scoring_function, scoring_args, contacted_residues, target_residues  # noqa: E402
from MolPropOp import (  # noqa: E402
    base_rings, clean_ring_locations, e_withdraw, e_donate,
    withdraw_with_linkers, donate_with_linkers,
)

SUBSTITUENT_POOL = e_withdraw + e_donate + withdraw_with_linkers + donate_with_linkers
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

    def __init__(self, protein, skip_docking, np_model=None, log_writer=None, log_meta=None):
        self.protein = protein
        self.skip_docking = skip_docking
        self.cache = {}  # canonical_smiles -> result dict
        self.n_dock_calls = 0
        self._np_model = np_model
        self._log_writer = log_writer  # every evaluated genotype, checkpointed immediately
        self._log_meta = log_meta or {}

    def evaluate(self, genome, generation):
        mol = Chem.MolFromSmiles(genome['smiles'])
        canon = Chem.MolToSmiles(mol)
        if canon in self.cache:
            return self.cache[canon]

        docking = docked_in_pocket = n_target_contacts = None
        if not self.skip_docking:
            scoring_args[0] = os.cpu_count()
            scoring_args[1] = self.protein
            try:
                score, aux = scoring_function(genome['smiles'])
                self.n_dock_calls += 1
                docking = score if aux is not None else None
                if aux is None:
                    docked_in_pocket = False
                else:
                    tgt_res = target_residues()
                    contacts = contacted_residues(aux)
                    if contacts is None:
                        docked_in_pocket = True
                    elif tgt_res:
                        n_target_contacts = len(contacts & tgt_res)
                        docked_in_pocket = n_target_contacts > 0
                    else:
                        docked_in_pocket = True
            except Exception:
                docking = None
                docked_in_pocket = False

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

        # unscored genotype: 0.0 (dockstring's own failure sentinel) so it
        # loses every tournament against any real (negative) docking score.
        fitness = docking if docking is not None else 0.0
        result = {
            'canonical_smiles': canon, 'original_smiles': genome['smiles'],
            'inchikey': Chem.MolToInchiKey(mol), 'docking': docking,
            'docked_in_pocket': docked_in_pocket, 'n_target_contacts': n_target_contacts,
            'qed': qed, 'alogp': alogp, 'sas': sas, 'np': np_score,
            'fitness': fitness,
        }
        result['generation'] = generation
        self.cache[canon] = result
        if self._log_writer is not None:
            row = dict(self._log_meta)
            row.update({k: v for k, v in result.items() if k != 'fitness'})
            row['protein'] = self.protein
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
    p.add_argument('--protein', default='HMGCR')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--mutation-rate', type=float, default=0.7)
    p.add_argument('--tournament-k', type=int, default=2)
    p.add_argument('--elitism', type=int, default=1)
    p.add_argument('--skip-docking', action='store_true',
                    help='CPU-light wiring check: build/mutate genomes but skip the real Vina call')
    p.add_argument('--out-dir', default=None,
                    help='default: results/batches/ga_baseline/<preset>/analysis')
    args = p.parse_args(argv)

    pop = args.pop or PRESETS[args.preset]['pop']
    gens = args.gens or PRESETS[args.preset]['gens']
    set_label = f'ga_{args.preset}'

    out_dir = args.out_dir or os.path.join(_ROOT, 'results', 'batches', 'ga_baseline', args.preset, 'analysis')
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, f'compounds_{set_label}.csv')
    log_csv = os.path.join(out_dir, f'ga_eval_log_{set_label}.csv')

    fieldnames = ['set_label', 'replicate', 'proposer_provider', 'proposer_model',
                  'adversary_provider', 'adversary_model', 'protein',
                  'original_smiles', 'canonical_smiles', 'inchikey', 'docking',
                  'docked_in_pocket', 'n_target_contacts', 'qed', 'alogp', 'sas', 'np',
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

    print(f"GA baseline: preset={args.preset} pop={pop} gens={gens} replicates={args.replicates} "
          f"protein={args.protein} skip_docking={args.skip_docking}")
    print(f"Delivered-compounds CSV: {out_csv}")
    print(f"Full evaluation log (checkpointed per dock call): {log_csv}")

    t0 = time.time()
    for rep in range(1, args.replicates + 1):
        rng = random.Random(args.seed + rep)
        meta = {
            'set_label': set_label, 'replicate': rep,
            'proposer_provider': None, 'proposer_model': 'genetic_algorithm',
            'adversary_provider': None, 'adversary_model': None,
        }
        evaluator = Evaluator(args.protein, args.skip_docking, np_model=np_model,
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
            row['protein'] = args.protein
            row['source_turn'] = gens
            out_writer.writerow(row)
        out_f.flush()

        best = fitness[id(ranked[0])]
        print(f"  rep {rep}/{args.replicates}: {evaluator.n_dock_calls} dock calls, "
              f"{len(seen)} delivered compounds, best={best['docking']} "
              f"({best['canonical_smiles']}) [{time.time() - t0:.0f}s elapsed]")

    out_f.close()
    log_f.close()
    print(f"Done. {time.time() - t0:.0f}s total.")


if __name__ == '__main__':
    main()
