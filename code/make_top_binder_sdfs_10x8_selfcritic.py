#!/usr/bin/env python3
"""make_top_binder_sdfs_10x8_selfcritic.py - Real-dock the single top binder from
each proposer for the two condition sets that lack pose SDFs:

  1. the 10x8 runs (10 molecules x 8 generations, 3 replicates, all vs OpenAI
     gpt-5.2 adversary) -> results/batches/sdf/10x8/
  2. the self-critic 5x4 runs (adversary = proposer's own model family)
     -> results/batches/sdf/self_critique/

Same recipe as make_top_binder_sdfs.py (the 5x4 vs-OpenAI set): each SDF's single
record carries the docked 3D pose plus SDF tags for condition, proposer/adversary,
rep, docking score, and QED.

Re-runs Vina once per compound (9 total) -- do NOT run this while any other
docking job (a batch run, another analyze_replicates.py pass) is active; Vina
uses all available CPUs and two concurrent jobs will contend for them.
"""
import csv
import os

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root
BATCHES = os.path.join(_ROOT, "results", "batches")
OUT = os.path.join(BATCHES, "sdf")

import docking_module as dm
from rdkit import Chem

dm.scoring_args[1] = 'HMGCR'  # default target is DRD2; this study docks to HMGCR

# condition label -> (output subdir, [(display label, best_per_replicate csv)])
CONDITIONS = {
    "10x8_vs_openai": ("10x8", [
        ("OpenAI gpt-5.2 / Anthropic haiku-4.5",
         "openai_gpt-5.2_vs_anthropic_10x8/analysis/best_per_replicate_openai_gpt-5.2_vs_anthropic_10x8.csv"),
        ("Anthropic haiku-4.5 / OpenAI gpt-5.2",
         "anthropic-haiku-4-5_vs_openai_10x8/analysis/best_per_replicate_anthropic-haiku-4-5_vs_openai_10x8.csv"),
        ("Gemini 3-flash / OpenAI gpt-5.2",
         "gemini-3-flash-preview_vs_openai_10x8/analysis/best_per_replicate_gemini-3-flash-preview_vs_openai_10x8.csv"),
        ("kimi k2.6 / OpenAI gpt-5.2",
         "ollama_kimi-k2.6_vs_openai_10x8/analysis/best_per_replicate_ollama_kimi-k2.6_vs_openai_10x8.csv"),
        ("deepseek v4-pro / OpenAI gpt-5.2",
         "ollama_deepseek-v4-pro_vs_openai_10x8/analysis/best_per_replicate_ollama_deepseek-v4-pro_vs_openai_10x8.csv"),
    ]),
    "self_critic_5x4": ("self_critique", [
        ("OpenAI gpt-5.2 / OpenAI gpt-5.2",
         "openai_gpt-5.2_vs_openai_5x4/analysis/best_per_replicate_openai_gpt-5.2_vs_openai_5x4.csv"),
        ("Anthropic haiku-4.5 / Anthropic haiku-4.5",
         "anthropic-haiku-4-5_vs_anthropic_5x4/analysis/best_per_replicate_anthropic-haiku-4-5_vs_anthropic_5x4.csv"),
        ("Gemini 3-flash / Gemini 3-flash",
         "gemini-3-flash-preview_vs_gemini_5x4/analysis/best_per_replicate_gemini-3-flash-preview_vs_gemini_5x4.csv"),
        ("kimi k2.6 / ollama",
         "ollama_kimi-k2.6_vs_ollama_5x4/analysis/best_per_replicate_ollama_kimi-k2.6_vs_ollama_5x4.csv"),
    ]),
}


def read_csv(path):
    with open(os.path.join(BATCHES, path)) as f:
        return list(csv.DictReader(f))


def slug(label):
    return label.split(" / ")[0].replace(" ", "_").replace(".", "")


def main():
    for condition, (subdir, proposers) in CONDITIONS.items():
        out_dir = os.path.join(OUT, subdir)
        os.makedirs(out_dir, exist_ok=True)
        for label, best_csv in proposers:
            rows = read_csv(best_csv)
            best_row = min(rows, key=lambda r: float(r["docking"]))
            smiles = best_row["canonical_smiles"]

            print(f"[{condition}] Docking {label}: {smiles} "
                  f"(prior score {best_row['docking']}, rep {best_row['replicate']}) ...")
            score, aux = dm.scoring_function(smiles)
            if aux is None:
                print(f"  DOCKING FAILED for {condition}/{label}, skipping.")
                continue

            pose_mol = aux['ligand']
            pose_mol.SetProp('_Name', f"{slug(label)}_{subdir}")
            pose_mol.SetProp('condition', condition)
            pose_mol.SetProp('proposer_adversary', label)
            pose_mol.SetProp('replicate', best_row['replicate'])
            pose_mol.SetProp('docking_score_redocked', str(score))
            pose_mol.SetProp('docking_score_analysis_csv', best_row['docking'])
            pose_mol.SetProp('QED', best_row['qed'])
            pose_mol.SetProp('canonical_smiles', smiles)

            out_path = os.path.join(out_dir, f"{slug(label)}_top_binder.sdf")
            w = Chem.SDWriter(out_path)
            w.write(pose_mol)
            w.close()
            print(f"  wrote {out_path} (re-docked score {score}, prior {best_row['docking']})")


if __name__ == '__main__':
    main()