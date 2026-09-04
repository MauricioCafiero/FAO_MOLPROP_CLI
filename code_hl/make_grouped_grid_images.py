#!/usr/bin/env python3
"""make_grouped_grid_images.py - HL-gap grid images, grouped the same way as the
docking study's results/batches/images/ tree.

The HL counterpart of code/make_grid_images.py (agentic) +
code/make_baseline_grid_images.py (zero/few/frag) + code/make_ga_grid_images.py
(GA). The earlier code_hl/make_grid_images.py dumped every batch flat into one
folder as accumulated_<batch>.png; this script reproduces the docking grouping:

  images/                          agentic 5x4 self-critique (top level, like docking)
    top_gap_per_model.png            one smallest-gap compound per proposer
    accumulated_<Proposer>.png       every compound a proposer produced
  images/<zero|few|frag>_shot/
    top_gap_per_model_<shot>.png     one smallest-gap compound per model (12 models)
    accumulated_<model>_<shot>.png   every compound a model produced
  images/ga_baseline/
    top_binders_ga_restricted.png    one best compound per replicate, restricted pool
    top_binders_ga_unrestricted.png  same, unrestricted pool
    accumulated_ga_restricted.png    every compound the restricted-pool GA produced
    accumulated_ga_unrestricted.png  same, unrestricted pool

The gemma4 self-critique batch (hl_gemma4-31b_vs_gemma4_5x4) is a test run and is
excluded, matching its exclusion from the statistics. Pure RDKit 2D depiction from
the canonical_smiles already in the analysis CSVs; no gap recomputation. Legends
report gap (eV) and SAS -- not QED, which is not an outcome metric in this study.

Usage:
  fao-env/bin/python code_hl/make_grouped_grid_images.py
"""
import csv
import os

from rdkit import Chem
from rdkit.Chem import Draw

_HERE = os.path.dirname(os.path.abspath(__file__))  # code_hl/
_ROOT = os.path.dirname(_HERE)  # repo root
HL = os.path.join(_ROOT, "results", "batches", "hl_batches")
OUT = os.path.join(HL, "images")

# Agentic 5x4 self-critique: (batch dir, display label, file slug).
# File slugs match the docking study's images/accumulated_*.png names.
AGENTIC = [
    ("hl_gpt-5.2_vs_gpt-5.2_5x4", "OpenAI gpt-5.2 (self-critic)", "OpenAI_gpt-52"),
    ("hl_claude-haiku-4-5_vs_claude-haiku-4-5_5x4", "Anthropic haiku-4.5 (self-critic)", "Anthropic_haiku-45"),
    ("hl_gemini-3-flash-preview_vs_gemini_5x4", "Gemini 3-flash (self-critic)", "Gemini_3-flash"),
    ("hl_kimi-k2.6_vs_kimi-k2.6_5x4", "kimi k2.6 (self-critic)", "kimi_k26"),
    ("hl_deepseek-v4-pro_vs_deepseek-v4-pro_5x4", "deepseek v4-pro (self-critic)", "deepseek_v4-pro"),
]

# zero/frag/few-shot: all 12 models live in one compounds_<shot>.csv, grouped
# on set_label (same layout the docking baseline script groups on).
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

# GA baselines: (compounds csv rel to HL, display label, pool slug).
GA = [
    ("ga_baseline/5x4/analysis/compounds_ga_5x4.csv", "GA (restricted)", "restricted"),
    ("ga_baseline/5x4_full/analysis/compounds_ga_5x4_full.csv", "GA (unrestricted)", "unrestricted"),
]


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def num(row, key):
    """float(row[key]) or None -- gap is blank when the calculation failed."""
    try:
        return float(row[key])
    except (TypeError, ValueError, KeyError):
        return None


def save_grid(mols, legends, path, per_row, size):
    img = Draw.MolsToGridImage(mols, molsPerRow=per_row, subImgSize=size,
                               legends=legends, useSVG=False)
    img.save(path)
    print(f"wrote {path} ({len(mols)} mols)")


def legend(row, prefix):
    gap = num(row, "gap")
    g = f"{gap:.3f} eV" if gap is not None else "gap failed"
    sas = num(row, "sas")
    s = f", SAS {sas:.2f}" if sas is not None else ""
    return f"{prefix}\n{g}{s}" if prefix else f"{g}{s}"


def main():
    os.makedirs(OUT, exist_ok=True)

    # --- Agentic 5x4 self-critique: top level, like docking's images/ ------
    top_mols, top_legends = [], []
    for batch, display, slug in AGENTIC:
        comp_csv = os.path.join(HL, batch, "analysis", f"compounds_{batch}.csv")
        rows = read_csv(comp_csv)
        scored = [r for r in rows if num(r, "gap") is not None]
        best = min(scored, key=lambda r: num(r, "gap"))
        mol = Chem.MolFromSmiles(best["canonical_smiles"])
        if mol is not None:
            top_mols.append(mol)
            top_legends.append(legend(best, display))
        mols, legends = [], []
        for i, row in enumerate(rows, start=1):
            mol = Chem.MolFromSmiles(row["canonical_smiles"])
            if mol is None:
                continue
            mols.append(mol)
            legends.append(legend(row, f"{i}"))
        if mols:
            save_grid(mols, legends, os.path.join(OUT, f"accumulated_{slug}.png"),
                      5, (220, 220))
    if top_mols:
        save_grid(top_mols, top_legends, os.path.join(OUT, "top_gap_per_model.png"),
                  3, (320, 320))

    # --- zero/few/frag-shot: one subfolder per shot type -------------------
    for shot in ["zero", "few", "frag"]:
        shot_dir = f"{shot}_shot"
        comp_csv = os.path.join(HL, shot_dir, "analysis", f"compounds_{shot_dir}.csv")
        out_dir = os.path.join(OUT, shot_dir)
        os.makedirs(out_dir, exist_ok=True)
        rows = read_csv(comp_csv)
        top_mols, top_legends = [], []
        for label in LABELS:
            sub = [r for r in rows if r["set_label"] == label]
            if not sub:
                print(f"  (no compounds for {label} in {shot_dir}, skipping)")
                continue
            # one accumulated grid per model
            mols, legends = [], []
            for i, row in enumerate(sub, start=1):
                mol = Chem.MolFromSmiles(row["canonical_smiles"])
                if mol is None:
                    continue
                mols.append(mol)
                legends.append(legend(row, f"{i}"))
            if mols:
                save_grid(mols, legends,
                          os.path.join(out_dir, f"accumulated_{label.replace('.', '')}_{shot_dir}.png"),
                          5, (220, 220))
            # this model's entry for the combined top-compound grid
            scored = [r for r in sub if num(r, "gap") is not None]
            best = min(scored, key=lambda r: num(r, "gap"))
            mol = Chem.MolFromSmiles(best["canonical_smiles"])
            if mol is not None:
                top_mols.append(mol)
                top_legends.append(legend(best, DISPLAY[label]))
        if top_mols:
            save_grid(top_mols, top_legends,
                      os.path.join(out_dir, f"top_gap_per_model_{shot_dir}.png"),
                      3, (320, 320))

    # --- GA baselines: one subfolder, one best-per-replicate + one
    #     accumulated grid per pool (same as docking's ga_baseline/) --------
    out_dir = os.path.join(OUT, "ga_baseline")
    os.makedirs(out_dir, exist_ok=True)
    for rel, display, pool in GA:
        rows = read_csv(os.path.join(HL, rel))
        scored = [r for r in rows if num(r, "gap") is not None]
        # one best compound per replicate
        by_rep = {}
        for row in scored:
            rep = row["replicate"]
            if rep not in by_rep or num(row, "gap") < num(by_rep[rep], "gap"):
                by_rep[rep] = row
        mols, legends = [], []
        for rep in sorted(by_rep, key=int):
            mol = Chem.MolFromSmiles(by_rep[rep]["canonical_smiles"])
            if mol is None:
                continue
            mols.append(mol)
            legends.append(legend(by_rep[rep], f"{display} r{rep}"))
        if mols:
            save_grid(mols, legends, os.path.join(out_dir, f"top_binders_ga_{pool}.png"),
                      3, (320, 320))
        # all accumulated compounds
        mols, legends = [], []
        for i, row in enumerate(rows, start=1):
            mol = Chem.MolFromSmiles(row["canonical_smiles"])
            if mol is None:
                continue
            mols.append(mol)
            legends.append(legend(row, f"{display} r{row['replicate']} {i}"))
        if mols:
            save_grid(mols, legends, os.path.join(out_dir, f"accumulated_ga_{pool}.png"),
                      5, (220, 220))


if __name__ == "__main__":
    main()