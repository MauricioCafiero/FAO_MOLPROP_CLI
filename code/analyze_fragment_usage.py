#!/usr/bin/env python3
"""
analyze_fragment_usage.py - For the frag-shot baseline (run_zero_few_shot.py
--shot frag), measure whether proposed compounds actually build on the
suggested ring skeletons / functional groups, not just their docking score.

The frag-shot prompt lists 9 base rings and 10 functional groups the model
"should use as the base of your molecules." The point of that prompt design
is that those fragments represent what a bench scientist can/wants to
synthesize -- so whether the model actually used them is the headline
question for this baseline, docking score is secondary.

Rings are checked as an RDKit substructure match (query mol, not SMARTS) --
this matches a ring system regardless of what's substituted onto it, which
is exactly "did the model build on this skeleton." A molecule built on
naphthalene also substructure-matches plain benzene (it contains one), so
each compound is also assigned a single "primary skeleton" = the largest
(most-atom) matching ring, to avoid double-counting.

Functional groups are checked the same way, but flagged by distinctiveness:
'I', 'C#N', 'C#C(SC)', 'C=C([N+](=O)[O-])', 'CC(N(C)C)' are structurally
specific tags; 'C(C)' and 'C(N)' are too generic (a bare ethyl/methylamine
substructure) to mean anything as a "used this functional group" signal --
reported separately, not folded into the headline rate.

Usage:
  fao-env/bin/python code/analyze_fragment_usage.py --compounds results/batches/frag_shot/analysis/compounds_frag_shot.csv
"""
import argparse
import csv
import os
from collections import defaultdict

from rdkit import Chem

# The 3 empirically best base rings by real Vina docking of every
# adversarial_set.md entry (flavone mean -7.33, anthracene mean -6.96,
# naphthalene mean -6.21 kcal/mol; the other 6 rings cluster at -4.55 to
# -5.20). "top3-ring rate" = the fraction of a proposer's compounds whose
# primary skeleton is one of these -- a ground-truthed test of whether a
# no-score mode (zero-shot, frag-shot) picks good fragments by intuition
# alone, vs. a with-scores mode (few-shot) that's handed the answer.
TOP3_RINGS = {'flavone', 'anthracene', 'naphthalene'}

RINGS = {
    'benzene': 'c1ccccc1',
    'pyridine': 'n1ccccc1',
    'furan': 'o1cccc1',
    'thiophene': 's1cccc1',
    'pyrrole': '[nH]1cccc1',
    'imidazole': 'n1c[nH]cc1',
    'naphthalene': 'c1ccc2ccccc2c1',
    'anthracene': 'c1ccc2cc3ccccc3cc2c1',
    'flavone': 'O=c1cc(-c2ccccc2)oc2ccccc12',
}

# distinctive (structurally specific) suggested functional groups
FGROUPS_SPECIFIC = {
    'iodo': 'I',
    'nitrile': 'C#N',
    'thioalkyne': 'C#C(SC)',
    'nitrovinyl': 'C=C([N+](=O)[O-])',
    'dimethylamino': 'CC(N(C)C)',
    'carboxylate_branch': 'C(C(=O)[O-])',
}
# generic (too unspecific to count as real signal -- reported separately)
FGROUPS_GENERIC = {
    'ethyl_branch': 'C(C)',
    'amino_branch': 'C(N)',
    'alkoxide_branch': 'C([O-])',
    'isopropyl_ester': 'C(=O)O(C(C)C)',
}


def _mol(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return m


RING_MOLS = {name: (_mol(smi), _mol(smi).GetNumAtoms()) for name, smi in RINGS.items()}
FG_SPECIFIC_MOLS = {name: _mol(smi) for name, smi in FGROUPS_SPECIFIC.items()}
FG_GENERIC_MOLS = {name: _mol(smi) for name, smi in FGROUPS_GENERIC.items()}


def ring_matches(mol):
    hits = []
    for name, (qmol, _n) in RING_MOLS.items():
        if qmol is not None and mol.HasSubstructMatch(qmol):
            hits.append(name)
    return hits


def primary_skeleton(hits):
    if not hits:
        return None
    return max(hits, key=lambda name: RING_MOLS[name][1])


def fg_matches(mol, fg_mols):
    hits = []
    for name, qmol in fg_mols.items():
        if qmol is not None and mol.HasSubstructMatch(qmol):
            hits.append(name)
    return hits


def analyze(compounds_csv):
    rows = list(csv.DictReader(open(compounds_csv)))
    per_set = defaultdict(lambda: {
        'n': 0, 'n_ring': 0, 'n_specific_fg': 0, 'n_top3': 0,
        'ring_counts': defaultdict(int), 'fg_counts': defaultdict(int),
        'primary_counts': defaultdict(int),
    })
    per_compound = []

    for r in rows:
        smi = r.get('canonical_smiles') or r.get('original_smiles')
        mol = _mol(smi) if smi else None
        label = r['set_label']
        d = per_set[label]
        d['n'] += 1
        if mol is None:
            per_compound.append({**r, 'ring_hits': '', 'primary_skeleton': '',
                                  'specific_fg_hits': '', 'generic_fg_hits': ''})
            continue
        rhits = ring_matches(mol)
        prim = primary_skeleton(rhits)
        sfg = fg_matches(mol, FG_SPECIFIC_MOLS)
        gfg = fg_matches(mol, FG_GENERIC_MOLS)
        if rhits:
            d['n_ring'] += 1
            for h in rhits:
                d['ring_counts'][h] += 1
            d['primary_counts'][prim] += 1
            if prim in TOP3_RINGS:
                d['n_top3'] += 1
        if sfg:
            d['n_specific_fg'] += 1
            for h in sfg:
                d['fg_counts'][h] += 1
        per_compound.append({
            **r,
            'ring_hits': ';'.join(rhits),
            'primary_skeleton': prim or '',
            'specific_fg_hits': ';'.join(sfg),
            'generic_fg_hits': ';'.join(gfg),
        })

    return per_set, per_compound


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--compounds', required=True,
                    help='Path to a compounds_<batch>.csv from analyze_replicates.py')
    p.add_argument('--out-csv', default=None,
                    help='Where to write the per-compound fragment-usage CSV '
                         '(default: alongside --compounds, suffixed _fragusage.csv)')
    args = p.parse_args()

    per_set, per_compound = analyze(args.compounds)

    out_csv = args.out_csv or os.path.splitext(args.compounds)[0] + '_fragusage.csv'
    fieldnames = list(per_compound[0].keys())
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(per_compound)

    print(f"{'set':18} {'n':4} {'skeleton%':10} {'top3-ring%':11} {'specific_fg%':13} top skeletons")
    for label, d in sorted(per_set.items()):
        skel_pct = 100.0 * d['n_ring'] / d['n'] if d['n'] else 0.0
        top3_pct = 100.0 * d['n_top3'] / d['n'] if d['n'] else 0.0
        fg_pct = 100.0 * d['n_specific_fg'] / d['n'] if d['n'] else 0.0
        top = sorted(d['primary_counts'].items(), key=lambda kv: -kv[1])[:3]
        top_str = ', '.join(f'{k}:{v}' for k, v in top)
        print(f"{label:18} {d['n']:4} {skel_pct:9.1f}% {top3_pct:10.1f}% {fg_pct:12.1f}% {top_str}")

    print(f"\nWrote per-compound fragment-usage CSV: {out_csv}")


if __name__ == '__main__':
    main()
