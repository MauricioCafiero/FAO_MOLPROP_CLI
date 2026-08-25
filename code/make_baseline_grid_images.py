#!/usr/bin/env python3
"""make_baseline_grid_images.py - Build grid images for the non-agentic
zero-shot / few-shot baseline (run_zero_few_shot.py), paralleling
make_grid_images.py's output for the agentic 5x4 study: one top-binder-per-
model grid + one accumulated-compounds grid per model (6 images per shot type).

Unlike the agentic 5x4 studies (one proposer per manifest/folder),
run_zero_few_shot.py writes all 5 models into a single batch_dir/manifest, so
analyze_replicates.py's output is one pair of CSVs per shot type with a
set_label column distinguishing models -- this script groups on that column
instead of reading 5 separate per-proposer CSVs.

Writes into results/batches/images/<zero_shot|few_shot>/ (private repo).

Usage:
  python3 code/make_baseline_grid_images.py --shot zero
  python3 code/make_baseline_grid_images.py --shot few
"""
import argparse
import csv
import os

from rdkit import Chem
from rdkit.Chem import Draw

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root
BATCHES = os.path.join(_ROOT, "results", "batches")

LABELS = ["openai", "anthropic", "gemini", "kimi-k2.6", "deepseek-v4-pro",
          "gemma4", "glm-5.2", "nemotron-3-ultra", "nemotron-3-super",
          "nemotron-3-nano", "gpt-oss-20b", "gpt-oss-120b"]
DISPLAY = {
    "openai": "OpenAI gpt-5.2",
    "anthropic": "Anthropic haiku-4.5",
    "gemini": "Gemini 3-flash",
    "kimi-k2.6": "kimi k2.6",
    "deepseek-v4-pro": "deepseek v4-pro",
    "gemma4": "gemma4",
    "glm-5.2": "glm-5.2",
    "nemotron-3-ultra": "nemotron-3-ultra",
    "nemotron-3-super": "nemotron-3-super",
    "nemotron-3-nano": "nemotron-3-nano",
    "gpt-oss-20b": "gpt-oss:20b",
    "gpt-oss-120b": "gpt-oss:120b",
}


def slug(label):
    return label.replace(" ", "_").replace(".", "")


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shot", choices=["zero", "few", "frag"], required=True)
    args = p.parse_args()

    shot_dir = f"{args.shot}_shot"
    analysis_dir = os.path.join(BATCHES, shot_dir, "analysis")
    best_csv = os.path.join(analysis_dir, f"best_per_replicate_{shot_dir}.csv")
    comp_csv = os.path.join(analysis_dir, f"compounds_{shot_dir}.csv")
    out_dir = os.path.join(BATCHES, "images", shot_dir)
    os.makedirs(out_dir, exist_ok=True)

    best_rows = read_csv(best_csv)
    comp_rows = read_csv(comp_csv)

    # --- One image: single top binder from each model present -------------
    top_mols, top_legends = [], []
    for label in LABELS:
        rows = [r for r in best_rows if r["set_label"] == label and r["docking"]]
        if not rows:
            print(f"  (no scored compounds for {label}, skipping from top-binder grid)")
            continue
        best_row = min(rows, key=lambda r: float(r["docking"]))
        mol = Chem.MolFromSmiles(best_row["canonical_smiles"])
        top_mols.append(mol)
        top_legends.append(
            f"{DISPLAY[label]}\ndock {best_row['docking']}, QED {float(best_row['qed']):.2f}")

    top_img = Draw.MolsToGridImage(top_mols, molsPerRow=3, subImgSize=(320, 320),
                                    legends=top_legends, useSVG=False)
    top_path = os.path.join(out_dir, f"top_binder_per_model_{shot_dir}.png")
    top_img.save(top_path)
    print(f"wrote {top_path} ({len(top_mols)} mols)")

    # --- Five images: one per model, all accumulated compounds -------------
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
        out_path = os.path.join(out_dir, f"accumulated_{slug(label)}_{shot_dir}.png")
        img.save(out_path)
        print(f"wrote {out_path} ({len(mols)} mols) -- {label}")


if __name__ == "__main__":
    main()
