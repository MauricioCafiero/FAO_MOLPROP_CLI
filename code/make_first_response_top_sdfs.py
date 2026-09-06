#!/usr/bin/env python3
"""make_first_response_top_sdfs.py - 3D SDFs for the first-response baseline's
top compound per proposer, paralleling make_baseline_top_binder_sdfs.py
(docking arm) with an HL-arm counterpart (no HL SDF script existed before).

  --study dock  re-docks each proposer's single best first-reply compound with
                the real Vina engine and writes the actual docking pose
                (docking_module.scoring_function), exactly like
                make_baseline_top_binder_sdfs.py. Re-runs Vina once per
                compound (5 total) -- do NOT run while another docking job is
                active; Vina uses all available CPUs.
  --study hl    the HL arm has no pose concept and hl_gap_module.scoring_function
                is single-point (no geometry optimization), so the SDF is an
                RDKit ETKDG-embedded, MMFF-optimized 3D structure of the
                smallest-gap first-reply compound -- a depiction, not a
                computed geometry.

Writes one SDF per proposer:
  results/batches/sdf/first_response_5x4/<Slug>_top_binder.sdf   (--study dock)
  results/batches/hl_batches/sdf/first_response_5x4/<Slug>_top_gap.sdf (--study hl)

Usage:
  python code/make_first_response_top_sdfs.py --study dock
  python code/make_first_response_top_sdfs.py --study hl
"""
import argparse
import csv
import os

_HERE = os.path.dirname(os.path.abspath(__file__))  # code/
_ROOT = os.path.dirname(_HERE)  # repo root
BATCHES = os.path.join(_ROOT, "results", "batches")

from rdkit import Chem
from rdkit.Chem import AllChem

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
    p = argparse.ArgumentParser()
    p.add_argument("--study", choices=["dock", "hl"], required=True)
    args = p.parse_args()

    if args.study == "dock":
        root, out_dir = os.path.join(BATCHES, "first_response_5x4"), os.path.join(BATCHES, "sdf", "first_response_5x4")
        metric, stem = "docking", "top_binder"
    else:
        root, out_dir = (os.path.join(BATCHES, "hl_batches", "first_response_5x4"),
                         os.path.join(BATCHES, "hl_batches", "sdf", "first_response_5x4"))
        metric, stem = "gap", "top_gap"
    os.makedirs(out_dir, exist_ok=True)

    rows = [r for r in read_csv(os.path.join(root, "analysis", "best_per_replicate_first_response_5x4.csv"))
            if r[metric] not in ("", None)]

    if args.study == "dock":
        import docking_module as dm
        dm.scoring_args[1] = 'HMGCR'  # default target is DRD2; this study docks to HMGCR

    for label in LABELS:
        label_rows = [r for r in rows if r["set_label"] == label]
        if not label_rows:
            print(f"  (no scored compounds for {label}, skipping)")
            continue
        best_row = min(label_rows, key=lambda r: float(r[metric]))
        smiles = best_row["canonical_smiles"]
        display_label = DISPLAY[label]
        print(f"{display_label} ({args.study}): {smiles} "
              f"(prior {metric} {best_row[metric]}, rep {best_row['replicate']}) ...")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  SMILES PARSE FAILED for {display_label}, skipping.")
            continue

        if args.study == "dock":
            score, aux = dm.scoring_function(smiles)
            if aux is None:
                print(f"  DOCKING FAILED for {display_label}, skipping.")
                continue
            out_mol = aux['ligand']
            out_mol.SetProp('docking_score_redocked', str(score))
            out_mol.SetProp('docking_score_analysis_csv', best_row[metric])
            out_mol.SetProp('QED', best_row['qed'])
        else:
            out_mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 0xC0FFEE
            if AllChem.EmbedMolecule(out_mol, params) != 0:
                print(f"  EMBEDDING FAILED for {display_label}, skipping.")
                continue
            AllChem.MMFFOptimizeMolecule(out_mol)
            out_mol = Chem.RemoveHs(out_mol)
            out_mol.SetProp('gap_analysis_csv', best_row[metric])
            out_mol.SetProp('SAS', best_row['sas'])

        out_mol.SetProp('_Name', slug(label))
        out_mol.SetProp('proposer', display_label)
        out_mol.SetProp('source', 'first_response_5x4')
        out_mol.SetProp('replicate', best_row['replicate'])
        out_mol.SetProp('canonical_smiles', smiles)

        out_path = os.path.join(out_dir, f"{slug(label)}_{stem}.sdf")
        w = Chem.SDWriter(out_path)
        w.write(out_mol)
        w.close()
        extra = f", re-docked score {score}" if args.study == "dock" else ""
        print(f"  wrote {out_path}{extra}")


if __name__ == '__main__':
    main()