#!/usr/bin/env python3
"""make_first_response_grid_images.py - Grid images for the first-response
baseline, paralleling make_baseline_grid_images.py (docking arm) and
make_grid_images.py (HL arm): one top-compound-per-model grid + one
accumulated-compounds grid per model, from the first_response_5x4 CSVs.

Both arms share the same baseline-style CSV schema (set_label / replicate /
canonical_smiles + metric + descriptors); the arm only changes the metric
column and legend fields:
  dock -> best = most negative docking, legend shows dock + QED
  hl   -> best = smallest gap,         legend shows gap + SAS

Writes into:
  results/batches/images/first_response_5x4/          (--study dock)
  results/batches/hl_batches/images/first_response_5x4/ (--study hl)

Usage:
  python code/make_first_response_grid_images.py --study dock
  python code/make_first_response_grid_images.py --study hl
"""
import argparse
import os

from rdkit import Chem
from rdkit.Chem import Draw

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root
BATCHES = os.path.join(_ROOT, "results", "batches")

LABELS = ["openai", "anthropic", "gemini", "kimi", "deepseek"]
DISPLAY = {
    "openai": "OpenAI gpt-5.2",
    "anthropic": "Anthropic haiku-4.5",
    "gemini": "Gemini 3-flash",
    "kimi": "kimi k2.6",
    "deepseek": "deepseek v4-pro",
}

ARMS = {
    "dock": dict(root=os.path.join(BATCHES, "first_response_5x4"),
                 out=os.path.join(BATCHES, "images", "first_response_5x4"),
                 metric="docking", metric_name="dock", fmt="{:.1f}",
                 legend="{score}, QED {aux:.2f}"),
    "hl": dict(root=os.path.join(BATCHES, "hl_batches", "first_response_5x4"),
               out=os.path.join(BATCHES, "hl_batches", "images", "first_response_5x4"),
               metric="gap", metric_name="gap", fmt="{:.2f}",
               legend="{score} eV, SAS {aux:.2f}"),
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--study", choices=["dock", "hl"], required=True)
    args = p.parse_args()
    cfg = ARMS[args.study]
    analysis = os.path.join(cfg["root"], "analysis")
    best_csv = os.path.join(analysis, f"best_per_replicate_first_response_5x4.csv")
    comp_csv = os.path.join(analysis, f"compounds_first_response_5x4.csv")
    out_dir = cfg["out"]
    os.makedirs(out_dir, exist_ok=True)

    import csv
    best_rows = list(csv.DictReader(open(best_csv)))
    comp_rows = list(csv.DictReader(open(comp_csv)))
    metric = cfg["metric"]

    def scored(row):
        return row[metric] not in ("", None)

    # --- One image: single top compound from each model --------------------
    top_mols, top_legends = [], []
    for label in LABELS:
        rows = [r for r in best_rows if r["set_label"] == label and scored(r)]
        if not rows:
            print(f"  (no scored compounds for {label}, skipping)")
            continue
        best_row = min(rows, key=lambda r: float(r[metric]))
        mol = Chem.MolFromSmiles(best_row["canonical_smiles"])
        if mol is None:
            continue
        aux = "qed" if metric == "docking" else "sas"
        top_mols.append(mol)
        top_legends.append(
            f"{DISPLAY[label]}\n"
            + cfg["legend"].format(score=cfg["fmt"].format(float(best_row[metric])),
                                   aux=float(best_row[aux])))

    top_img = Draw.MolsToGridImage(top_mols, molsPerRow=3, subImgSize=(320, 320),
                                   legends=top_legends, useSVG=False)
    stem = "top_binder_per_model" if metric == "docking" else "top_gap_per_model"
    top_path = os.path.join(out_dir, f"{stem}_first_response_5x4.png")
    top_img.save(top_path)
    print(f"wrote {top_path} ({len(top_mols)} mols)")

    # --- One image per model: all accumulated compounds --------------------
    for label in LABELS:
        rows = [r for r in comp_rows if r["set_label"] == label]
        mols, legends = [], []
        for i, row in enumerate(rows, start=1):
            mol = Chem.MolFromSmiles(row["canonical_smiles"])
            if mol is None:
                continue
            mols.append(mol)
            legends.append(f"{DISPLAY[label]}\n{i}")
        if not mols:
            print(f"  (no compounds for {label}, skipping accumulated grid)")
            continue
        img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(220, 220),
                                   legends=legends, useSVG=False)
        out_path = os.path.join(out_dir, f"accumulated_{label}_first_response_5x4.png")
        img.save(out_path)
        print(f"wrote {out_path} ({len(mols)} mols) -- {label}")


if __name__ == "__main__":
    main()