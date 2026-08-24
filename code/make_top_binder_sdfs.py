#!/usr/bin/env python3
"""make_top_binder_sdfs.py - Real-dock the single top binder from each of the 5
proposers (the 5-proposer / 5x4 study, PROPOSER_COMPARISON_5.md) and write each
one's actual Vina docking pose to its own SDF.

Writes into results/batches/sdf/ (private repo -- data only, no code), one file
per proposer: <slug>_top_binder.sdf. Each SDF's single record carries the docked
3D pose plus SDF tags for proposer/adversary, rep, docking score, and QED.

Re-runs Vina once per compound (5 total) -- do NOT run this while any other
docking job (a batch run, another analyze_replicates.py pass) is active; Vina
uses all available CPUs and two concurrent jobs will contend for them.
"""
import csv
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root
BATCHES = os.path.join(_ROOT, "results", "batches")
OUT = os.path.join(BATCHES, "sdf")

import docking_module as dm
from rdkit import Chem

dm.scoring_args[1] = 'HMGCR'  # default target is DRD2; this study docks to HMGCR

# (display label, best_per_replicate csv)
PROPOSERS = [
    ("OpenAI gpt-5.2 / Anthropic haiku-4.5",
     "openai_gpt-5.2_vs_anthropic_5x4/analysis/best_per_replicate_e2e_test.csv"),
    ("Anthropic haiku-4.5 / OpenAI gpt-5.2",
     "anthropic-haiku-4-5_vs_openai_5x4/analysis/best_per_replicate_e2e_test.csv"),
    ("Gemini 3-flash / OpenAI gpt-5.2",
     "gemini-3-flash-preview_vs_openai_5x4/analysis_full5/best_per_replicate_gemini-3-flash-preview_vs_openai.csv"),
    ("kimi k2.6 / OpenAI gpt-5.2",
     "ollama_kimi-k2.6_vs_openai_5x4/analysis_full5/best_per_replicate_ollama_kimi-k2.6_vs_openai.csv"),
    ("deepseek v4-pro / OpenAI gpt-5.2",
     "ollama_deepseek-v4-pro_vs_openai_5x4/analysis_full5/best_per_replicate_ollama_deepseek-v4-pro_vs_openai.csv"),
]


def read_csv(path):
    with open(os.path.join(BATCHES, path)) as f:
        return list(csv.DictReader(f))


def slug(label):
    return label.split(" / ")[0].replace(" ", "_").replace(".", "")


def main():
    os.makedirs(OUT, exist_ok=True)
    for label, best_csv in PROPOSERS:
        rows = read_csv(best_csv)
        best_row = min(rows, key=lambda r: float(r["docking"]))
        smiles = best_row["canonical_smiles"]

        print(f"Docking {label}: {smiles} (prior score {best_row['docking']}, rep {best_row['replicate']}) ...")
        score, aux = dm.scoring_function(smiles)
        if aux is None:
            print(f"  DOCKING FAILED for {label}, skipping.")
            continue

        pose_mol = aux['ligand']
        pose_mol.SetProp('_Name', slug(label))
        pose_mol.SetProp('proposer_adversary', label)
        pose_mol.SetProp('replicate', best_row['replicate'])
        pose_mol.SetProp('docking_score_redocked', str(score))
        pose_mol.SetProp('docking_score_analysis_csv', best_row['docking'])
        pose_mol.SetProp('QED', best_row['qed'])
        pose_mol.SetProp('canonical_smiles', smiles)

        out_path = os.path.join(OUT, f"{slug(label)}_top_binder.sdf")
        w = Chem.SDWriter(out_path)
        w.write(pose_mol)
        w.close()
        print(f"  wrote {out_path} (re-docked score {score}, prior {best_row['docking']})")


if __name__ == '__main__':
    main()
