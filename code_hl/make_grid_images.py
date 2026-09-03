#!/usr/bin/env python3
"""make_grid_images.py - Grid images for the HL-gap study, from analyze_replicates.py CSVs.

The HL counterpart of code/make_grid_images.py. Unlike that one (which hardcodes the five
docking proposer paths), this discovers batches automatically: any
<results-root>/<batch>/analysis/ holding compounds_*.csv is picked up, so agentic sets and
the zero/frag/few-shot baselines all work without editing this file.

Writes into <results-root>/images/:
  top_gap_per_batch.png      - the single smallest-gap compound from each batch
  accumulated_<batch>.png    - every compound a batch produced, one grid per batch

Pure RDKit 2D depiction from canonical_smiles in the CSVs; no gap recomputation.

Usage:
  python code_hl/make_grid_images.py
  python code_hl/make_grid_images.py --results-root results/batches/hl_batches --only kimi
"""
import os
import csv
import glob
import argparse
from rdkit import Chem
from rdkit.Chem import Draw

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_DEFAULT_ROOT = os.path.join(_ROOT, 'results', 'batches', 'hl_batches')


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def _num(row, key):
    """float(row[key]) or None -- gap is blank when the calculation failed."""
    try:
        return float(row[key])
    except (TypeError, ValueError, KeyError):
        return None


def find_batches(root, only=None):
    """[(batch_name, compounds_csv, best_csv|None)] for every analysed batch."""
    out = []
    for comp in sorted(glob.glob(os.path.join(root, '*', 'analysis', 'compounds_*.csv'))):
        name = os.path.basename(os.path.dirname(os.path.dirname(comp)))
        if only and only not in name:
            continue
        best = comp.replace('compounds_', 'best_per_replicate_')
        out.append((name, comp, best if os.path.isfile(best) else None))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description='Grid images for HL-gap batches.')
    p.add_argument('--results-root', default=_DEFAULT_ROOT)
    p.add_argument('--only', default=None, help='Substring filter on batch name.')
    p.add_argument('--out-dir', default=None, help='Default: <results-root>/images')
    args = p.parse_args(argv)

    batches = find_batches(args.results_root, args.only)
    if not batches:
        print(f'No analysed batches under {args.results_root} '
              f'(run analyze_replicates.py first).')
        return 1
    out_dir = args.out_dir or os.path.join(args.results_root, 'images')
    os.makedirs(out_dir, exist_ok=True)

    # --- One image: smallest-gap compound from each batch ------------------
    top_mols, top_legends = [], []
    for name, comp_csv, best_csv in batches:
        rows = [r for r in read_csv(best_csv or comp_csv) if _num(r, 'gap') is not None]
        if not rows:
            print(f'  [skip {name}: no compound with a gap]')
            continue
        best = min(rows, key=lambda r: _num(r, 'gap'))
        mol = Chem.MolFromSmiles(best['canonical_smiles'])
        if mol is None:
            continue
        qed = _num(best, 'qed')
        top_mols.append(mol)
        top_legends.append(f"{name}\ngap {_num(best,'gap'):.3f} eV"
                           + (f", QED {qed:.2f}" if qed is not None else ''))
    if top_mols:
        img = Draw.MolsToGridImage(top_mols, molsPerRow=3, subImgSize=(320, 320),
                                   legends=top_legends, useSVG=False)
        path = os.path.join(out_dir, 'top_gap_per_batch.png')
        img.save(path)
        print(f'wrote {path} ({len(top_mols)} mols)')

    # --- One image per batch: every compound it produced --------------------
    for name, comp_csv, _ in batches:
        mols, legends = [], []
        for i, row in enumerate(read_csv(comp_csv), start=1):
            mol = Chem.MolFromSmiles(row['canonical_smiles'])
            if mol is None:
                continue
            gap = _num(row, 'gap')
            mols.append(mol)
            legends.append(f"{i}: {gap:.3f} eV" if gap is not None else f'{i}: gap failed')
        if not mols:
            continue
        img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(220, 220),
                                   legends=legends, useSVG=False)
        path = os.path.join(out_dir, f'accumulated_{name}.png')
        img.save(path)
        print(f'wrote {path} ({len(mols)} mols)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
