#!/usr/bin/env python3
"""make_first_response_top_sdfs.py - Re-dock the single top binder from each of
the 5 proposers in the first-response baseline and write each one's actual Vina
docking pose to its own SDF, paralleling make_baseline_top_binder_sdfs.py but
reading the first-response baseline's CSV.

Docking-arm only by design: an SDF pose is a docking artifact; the HL-gap arm
has no pose concept and gets no SDFs.

Writes into results/batches/sdf/first_response_5x4/ (private repo -- data only,
no code), one file per proposer: <Slug>_top_binder.sdf.

Re-runs Vina once per compound (5 total) -- do NOT run this while any other
docking job is active; Vina uses all available CPUs.

Usage:
  python3 code/make_first_response_top_sdfs.py
"""
import csv
import os

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root
BATCHES = os.path.join(_ROOT, "results", "batches")

import docking_module as dm
from rdkit import Chem

dm.scoring_args[1] = 'HMGCR'  # default target is DRD2; this study docks to HMGCR

LABELS = ["openai", "anthropic", "gemini", "kimi", "deepseek"]
DISPLAY = {
    "openai": "OpenAI gpt-5.2",
    "anthropic": "Anthropic haiku-4.5",
    "gemini": "Gemini 3-flash",
    "kimi": "kimi k2.6",
    "deepseek": "deepseek v4-pro",
}


def slug(label):
    return DISPLAY[label].split(" ")[0].replace(".", "")


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def main():
    root = os.path.join(BATCHES, "first_response_5x4")
    out_dir = os.path.join(BATCHES, "sdf", "first_response_5x4")
    os.makedirs(out_dir, exist_ok=True)

    rows = [r for r in read_csv(os.path.join(root, "analysis", "best_per_replicate_first_response_5x4.csv"))
            if r["docking"] not in ("", None)]

    for label in LABELS:
        label_rows = [r for r in rows if r["set_label"] == label]
        if not label_rows:
            print(f"  (no scored compounds for {label}, skipping)")
            continue
        best_row = min(label_rows, key=lambda r: float(r["docking"]))
        smiles = best_row["canonical_smiles"]
        display_label = DISPLAY[label]

        print(f"Docking {display_label} (first_response_5x4): {smiles} "
              f"(prior score {best_row['docking']}, rep {best_row['replicate']}) ...")
        score, aux = dm.scoring_function(smiles)
        if aux is None:
            print(f"  DOCKING FAILED for {display_label}, skipping.")
            continue

        pose_mol = aux['ligand']
        pose_mol.SetProp('_Name', slug(label))
        pose_mol.SetProp('proposer', display_label)
        pose_mol.SetProp('source', 'first_response_5x4')
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