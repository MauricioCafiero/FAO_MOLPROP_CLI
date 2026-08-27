#!/usr/bin/env python3
"""make_ga_grid_images.py - Build grid images for the GA baseline
(run_ga_baseline.py), paralleling make_baseline_grid_images.py: one
best-per-replicate image per pool + one accumulated-compounds grid per pool
(4 images total). ga_5x4 = restricted (frag10) pool, ga_5x4_full =
unrestricted (full) pool.

There is no best_per_replicate CSV for the GA baseline, so the best-docking
compound per replicate is picked directly from each pool's compounds CSV.
One pool = one "set_label" = one accumulated image.

Writes into results/batches/images/ga_baseline/ (private repo).
"""
import csv
import os

from rdkit import Chem
from rdkit.Chem import Draw

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root
BATCHES = os.path.join(_ROOT, "results", "batches")

SETS = ["ga_5x4", "ga_5x4_full"]
DISPLAY = {
    "ga_5x4": "GA (restricted)",
    "ga_5x4_full": "GA (unrestricted)",
}


def main():
    out_dir = os.path.join(BATCHES, "images", "ga_baseline")
    os.makedirs(out_dir, exist_ok=True)

    all_rows = {}
    for label in SETS:
        path = os.path.join(BATCHES, "ga_baseline", label.replace("ga_", ""),
                            "analysis", f"compounds_{label}.csv")
        with open(path) as f:
            all_rows[label] = [r for r in csv.DictReader(f) if r["docking"]]

    # --- One image per pool: single best binder per replicate --------------
    for label in SETS:
        by_rep = {}
        for row in all_rows[label]:
            rep = row["replicate"]
            if rep not in by_rep or float(row["docking"]) < float(by_rep[rep]["docking"]):
                by_rep[rep] = row
        mols, legends = [], []
        for rep in sorted(by_rep, key=int):
            row = by_rep[rep]
            mol = Chem.MolFromSmiles(row["canonical_smiles"])
            if mol is None:
                continue
            mols.append(mol)
            legends.append(f"{DISPLAY[label]} r{rep}\n"
                           f"dock {row['docking']}, QED {float(row['qed']):.2f}")
        if not mols:
            print(f"  (no scored compounds for {label}, skipping best-per-replicate grid)")
            continue
        img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(320, 320),
                                   legends=legends, useSVG=False)
        slug = "restricted" if label == "ga_5x4" else "unrestricted"
        out_path = os.path.join(out_dir, f"top_binders_ga_{slug}.png")
        img.save(out_path)
        print(f"wrote {out_path} ({len(mols)} mols, one per replicate)")

    # --- One image per pool: all accumulated compounds ---------------------
    for label in SETS:
        mols, legends = [], []
        for i, row in enumerate(all_rows[label], start=1):
            mol = Chem.MolFromSmiles(row["canonical_smiles"])
            if mol is None:
                continue
            mols.append(mol)
            legends.append(f"{DISPLAY[label]} r{row['replicate']}\n{i}")
        if not mols:
            print(f"  (no compounds for {label}, skipping accumulated grid)")
            continue
        img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(220, 220),
                                   legends=legends, useSVG=False)
        slug = "restricted" if label == "ga_5x4" else "unrestricted"
        out_path = os.path.join(out_dir, f"accumulated_ga_{slug}.png")
        img.save(out_path)
        print(f"wrote {out_path} ({len(mols)} mols)")


if __name__ == "__main__":
    main()