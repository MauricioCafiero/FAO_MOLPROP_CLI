#!/usr/bin/env python3
"""make_grid_images.py - Build grid images for the 5-proposer / 5x4 study
(PROPOSER_COMPARISON_5.md) from analyze_replicates.py's output CSVs.

Writes, into results/batches/images/ (private repo -- images only, no code):
  1. top_binder_per_model.png  - one grid, the single best-scoring compound from
     each of the 5 proposers (legend: proposer/adversary + docking score + QED).
  2. accumulated_<proposer>.png (x5) - one grid per proposer of every compound it
     produced across all 5 replicates (legend: index number only).

No docking -- pure RDKit 2D depiction from the canonical_smiles already in the
analysis CSVs (compounds_*.csv / best_per_replicate_*.csv).
"""
import csv
import os
from rdkit import Chem
from rdkit.Chem import Draw

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root
BATCHES = os.path.join(_ROOT, "results", "batches")
OUT = os.path.join(BATCHES, "images")

# (display label, best_per_replicate csv, compounds csv) -- paths relative to BATCHES
PROPOSERS = [
    ("OpenAI gpt-5.2 / Anthropic haiku-4.5",
     "openai_gpt-5.2_vs_anthropic_5x4/analysis/best_per_replicate_e2e_test.csv",
     "openai_gpt-5.2_vs_anthropic_5x4/analysis/compounds_e2e_test.csv"),
    ("Anthropic haiku-4.5 / OpenAI gpt-5.2",
     "anthropic-haiku-4-5_vs_openai_5x4/analysis/best_per_replicate_e2e_test.csv",
     "anthropic-haiku-4-5_vs_openai_5x4/analysis/compounds_e2e_test.csv"),
    ("Gemini 3-flash / OpenAI gpt-5.2",
     "gemini-3-flash-preview_vs_openai_5x4/analysis_full5/best_per_replicate_gemini-3-flash-preview_vs_openai.csv",
     "gemini-3-flash-preview_vs_openai_5x4/analysis_full5/compounds_gemini-3-flash-preview_vs_openai.csv"),
    ("kimi k2.6 / OpenAI gpt-5.2",
     "ollama_kimi-k2.6_vs_openai_5x4/analysis_full5/best_per_replicate_ollama_kimi-k2.6_vs_openai.csv",
     "ollama_kimi-k2.6_vs_openai_5x4/analysis_full5/compounds_ollama_kimi-k2.6_vs_openai.csv"),
    ("deepseek v4-pro / OpenAI gpt-5.2",
     "ollama_deepseek-v4-pro_vs_openai_5x4/analysis_full5/best_per_replicate_ollama_deepseek-v4-pro_vs_openai.csv",
     "ollama_deepseek-v4-pro_vs_openai_5x4/analysis_full5/compounds_ollama_deepseek-v4-pro_vs_openai.csv"),
]


def read_csv(path):
    with open(os.path.join(BATCHES, path)) as f:
        return list(csv.DictReader(f))


def slug(label):
    return label.split(" / ")[0].replace(" ", "_").replace(".", "")


def main():
    os.makedirs(OUT, exist_ok=True)

    # --- One image: single top binder from each of the 5 proposers ---------
    top_mols, top_legends = [], []
    for label, best_csv, _ in PROPOSERS:
        rows = read_csv(best_csv)
        best_row = min(rows, key=lambda r: float(r["docking"]))
        mol = Chem.MolFromSmiles(best_row["canonical_smiles"])
        top_mols.append(mol)
        top_legends.append(
            f"{label}\ndock {best_row['docking']}, QED {float(best_row['qed']):.2f}")

    top_img = Draw.MolsToGridImage(top_mols, molsPerRow=3, subImgSize=(320, 320),
                                    legends=top_legends, useSVG=False)
    top_path = os.path.join(OUT, "top_binder_per_model.png")
    top_img.save(top_path)
    print(f"wrote {top_path} ({len(top_mols)} mols)")

    # --- Five images: one per proposer, all accumulated compounds ----------
    for label, _, comp_csv in PROPOSERS:
        rows = read_csv(comp_csv)
        mols, legends = [], []
        for i, row in enumerate(rows, start=1):
            mol = Chem.MolFromSmiles(row["canonical_smiles"])
            if mol is None:
                continue
            mols.append(mol)
            legends.append(f"{label}\n{i}")
        img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(220, 220),
                                    legends=legends, useSVG=False)
        out_path = os.path.join(OUT, f"accumulated_{slug(label)}.png")
        img.save(out_path)
        print(f"wrote {out_path} ({len(mols)} mols) -- {label}")


if __name__ == '__main__':
    main()
