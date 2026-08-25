#!/usr/bin/env python3
"""make_ablation_grid_images.py - Build grid images for the self-critique
ablation (SELF_CRITIQUE_ABLATION.md) and the 10x8-protocol sets
(PROTOCOL_DIFF_5x4_vs_10x8.md), paralleling make_grid_images.py's output for
the main 5x4 study: one top-binder-per-model grid + one accumulated-compounds
grid per proposer.

Both studies use one proposer per batch_dir/manifest (like the main 5x4
study), so this reads batch_dir/analysis/{best_per_replicate,compounds}_
<batch_dir>.csv directly -- no set_label grouping needed.

Writes into results/batches/images/<self_critique|10x8>/ (private repo).

Usage:
  python3 code/make_ablation_grid_images.py --study self
  python3 code/make_ablation_grid_images.py --study 10x8
"""
import argparse
import csv
import os

from rdkit import Chem
from rdkit.Chem import Draw

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root
BATCHES = os.path.join(_ROOT, "results", "batches")

# (display label, batch_dir) -- csv paths are derived as
# batch_dir/analysis/{best_per_replicate,compounds}_<batch_dir>.csv
STUDIES = {
    "self": [
        ("Anthropic haiku-4.5 (self)", "anthropic-haiku-4-5_vs_anthropic_5x4"),
        ("Gemini 3-flash (self)", "gemini-3-flash-preview_vs_gemini_5x4"),
        ("kimi k2.6 (self)", "ollama_kimi-k2.6_vs_ollama_5x4"),
        ("OpenAI gpt-5.2 (self)", "openai_gpt-5.2_vs_openai_5x4"),
    ],
    "10x8": [
        ("OpenAI gpt-5.2 / Anthropic haiku-4.5", "openai_gpt-5.2_vs_anthropic_10x8"),
        ("Anthropic haiku-4.5 / OpenAI gpt-5.2", "anthropic-haiku-4-5_vs_openai_10x8"),
        ("Gemini 3-flash / OpenAI gpt-5.2", "gemini-3-flash-preview_vs_openai_10x8"),
        ("kimi k2.6 / OpenAI gpt-5.2", "ollama_kimi-k2.6_vs_openai_10x8"),
        ("deepseek v4-pro / OpenAI gpt-5.2", "ollama_deepseek-v4-pro_vs_openai_10x8"),
    ],
}


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def slug(label):
    name = label.split(" / ")[0].replace(" (self)", "")
    return name.replace(" ", "_").replace(".", "")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--study", choices=["self", "10x8"], required=True)
    args = p.parse_args()

    out_dir_name = "self_critique" if args.study == "self" else "10x8"
    out_dir = os.path.join(BATCHES, "images", out_dir_name)
    os.makedirs(out_dir, exist_ok=True)

    proposers = STUDIES[args.study]

    # --- One image: single top binder from each proposer -------------------
    top_mols, top_legends = [], []
    for label, batch_dir in proposers:
        best_csv = os.path.join(BATCHES, batch_dir, "analysis", f"best_per_replicate_{batch_dir}.csv")
        rows = [r for r in read_csv(best_csv) if r["docking"]]
        if not rows:
            print(f"  (no scored compounds for {label}, skipping from top-binder grid)")
            continue
        best_row = min(rows, key=lambda r: float(r["docking"]))
        mol = Chem.MolFromSmiles(best_row["canonical_smiles"])
        top_mols.append(mol)
        top_legends.append(
            f"{label}\ndock {best_row['docking']}, QED {float(best_row['qed']):.2f}")

    top_img = Draw.MolsToGridImage(top_mols, molsPerRow=3, subImgSize=(320, 320),
                                    legends=top_legends, useSVG=False)
    top_path = os.path.join(out_dir, f"top_binder_per_model_{out_dir_name}.png")
    top_img.save(top_path)
    print(f"wrote {top_path} ({len(top_mols)} mols)")

    # --- One image per proposer: all accumulated compounds ------------------
    for label, batch_dir in proposers:
        comp_csv = os.path.join(BATCHES, batch_dir, "analysis", f"compounds_{batch_dir}.csv")
        rows = read_csv(comp_csv)
        mols, legends = [], []
        for i, row in enumerate(rows, start=1):
            mol = Chem.MolFromSmiles(row["canonical_smiles"])
            if mol is None:
                continue
            mols.append(mol)
            legends.append(f"{label}\n{i}")
        if not mols:
            print(f"  (no compounds for {label}, skipping accumulated grid)")
            continue
        img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(220, 220),
                                    legends=legends, useSVG=False)
        out_path = os.path.join(out_dir, f"accumulated_{slug(label)}_{out_dir_name}.png")
        img.save(out_path)
        print(f"wrote {out_path} ({len(mols)} mols) -- {label}")


if __name__ == "__main__":
    main()
