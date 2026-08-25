#!/usr/bin/env python3
"""make_baseline_top_binder_sdfs.py - Real-dock the single top binder from each
of the 5 proposers in the non-agentic zero-shot / few-shot baseline
(run_zero_few_shot.py) and write each one's actual Vina docking pose to its
own SDF, paralleling make_top_binder_sdfs.py's output for the agentic 5x4
study.

Unlike the agentic study (one proposer per best_per_replicate CSV),
run_zero_few_shot.py writes all 5 models into a single best_per_replicate_
<shot>.csv with a set_label column, so this reads one CSV and groups on that
column instead of reading 5 separate per-proposer CSVs.

Writes into results/batches/sdf/<zero_shot|few_shot>/ (private repo -- data
only, no code), one file per proposer: <slug>_top_binder.sdf.

Re-runs Vina once per compound (5 total per shot type) -- do NOT run this
while any other docking job is active; Vina uses all available CPUs.

Usage:
  python3 code/make_baseline_top_binder_sdfs.py --shot zero
  python3 code/make_baseline_top_binder_sdfs.py --shot few
"""
import argparse
import csv
import os

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root
BATCHES = os.path.join(_ROOT, "results", "batches")

import docking_module as dm
from rdkit import Chem

dm.scoring_args[1] = 'HMGCR'  # default target is DRD2; this study docks to HMGCR

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
    "gpt-oss-20b": "gpt-oss-20b",
    "gpt-oss-120b": "gpt-oss-120b",
}


def slug(label):
    return DISPLAY[label].split(" ")[0].replace(".", "")


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shot", choices=["zero", "few", "frag"], required=True)
    args = p.parse_args()

    shot_dir = f"{args.shot}_shot"
    best_csv = os.path.join(BATCHES, shot_dir, "analysis", f"best_per_replicate_{shot_dir}.csv")
    out_dir = os.path.join(BATCHES, "sdf", shot_dir)
    os.makedirs(out_dir, exist_ok=True)

    rows = read_csv(best_csv)

    for label in LABELS:
        label_rows = [r for r in rows if r["set_label"] == label and r["docking"]]
        if not label_rows:
            print(f"  (no scored compounds for {label}, skipping)")
            continue
        best_row = min(label_rows, key=lambda r: float(r["docking"]))
        smiles = best_row["canonical_smiles"]
        display_label = DISPLAY[label]

        print(f"Docking {display_label} ({shot_dir}): {smiles} "
              f"(prior score {best_row['docking']}, rep {best_row['replicate']}) ...")
        score, aux = dm.scoring_function(smiles)
        if aux is None:
            print(f"  DOCKING FAILED for {display_label}, skipping.")
            continue

        pose_mol = aux['ligand']
        pose_mol.SetProp('_Name', slug(label))
        pose_mol.SetProp('proposer', display_label)
        pose_mol.SetProp('shot_type', shot_dir)
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
