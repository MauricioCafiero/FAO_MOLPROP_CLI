#!/usr/bin/env python3
"""make_orca_validation_subset.py - Prepare the ORCA cross-validation subset for
the HL-gap agentic 5x4 self-critic study: the best molecule of every replicate
from each of the five proposer batches, deduplicated by InChIKey.

For every unique molecule this writes (into
results/batches/hl_batches/orca_validation/agentic_5x4_top/<InChIKey>/):
  geom.xyz  -- the SAME geometry the study's GFN2-xTB gaps were computed on
               (RDKit ETKDGv3 seed 0xf00d + MMFF optimize, via
               hl_gap_module._smiles_to_xyz_and_charge)
  <key>.inp -- ORCA 6 singlepoint: WB97M-V def2-tzvp RIJCOSX defgrid3 tightscf,
               4 procs, maxcore 4000 (per-machine limits: 4 cores / 4 GB)

plus a manifest CSV (agentic_5x4_top/manifest.csv) with one row per unique
molecule: proposer(s), replicate(s), min xTB gap, atom count, directory.
The manifest is sorted ascending by xTB gap so the run order puts the
best molecules first. Running is a separate step (run_orca_validation.sh).

No ORCA runs happen here -- geometry generation and input writing only.
"""
import csv
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))  # code_hl/
_ROOT = os.path.dirname(_HERE)  # FAO_MOLPROP_CLI root
BATCHES = os.path.join(_ROOT, 'results', 'batches', 'hl_batches')
OUT = os.path.join(BATCHES, 'orca_validation', 'agentic_5x4_top')

sys.path.insert(0, _HERE)
import hl_gap_module as h  # geometry + charge/spin convention

# (batch dir suffix, analysis-tag suffix, proposer display)
BATCHES_5X4 = [
    ('gpt-5.2', 'gpt-5.2', 'gpt-5.2'),
    ('claude-haiku-4-5', 'claude-haiku-4-5', 'haiku-4.5'),
    ('gemini-3-flash-preview', 'gemini', 'Gemini 3-flash'),
    ('kimi-k2.6', 'kimi-k2.6', 'kimi k2.6'),
    ('deepseek-v4-pro', 'deepseek-v4-pro', 'deepseek v4-pro'),
]

_INP_TEMPLATE = """! WB97M-V def2-tzvp RIJCOSX defgrid3 tightscf
%pal nprocs 4 end
%scf maxcore 4000 end
%output print[P_Basis] 2 print[P_Mulliken] 1 end
* xyzfile {charge} {mult} geom.xyz
"""

# rough valence-electron table for the multiplicity check (matches hl_gap_module)
_VALENCE = {'H': 1, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Si': 4, 'P': 5,
            'S': 6, 'Cl': 7, 'Br': 7, 'I': 7}


def main():
    os.makedirs(OUT, exist_ok=True)
    seen = {}  # inchikey -> row
    for sfx, tag, disp in BATCHES_5X4:
        path = os.path.join(BATCHES, f'hl_{sfx}_vs_{tag}_5x4', 'analysis',
                            f'best_per_replicate_hl_{sfx}_vs_{tag}_5x4.csv')
        with open(path) as f:
            for r in csv.DictReader(f):
                if not r['canonical_smiles'] or not r['gap']:
                    continue
                key = r['canonical_smiles']
                ik, charge, mult, symbols, xyz_text = mol_from_smiles(key)
                if ik is None:
                    print(f'  (skip, RDKit parse failed: {key})')
                    continue
                gap = float(r['gap'])
                if ik not in seen:
                    seen[ik] = dict(inchikey=ik, smiles=key, xtb_gap=gap,
                                    proposers=disp, replicates=str(r['replicate']),
                                    charge=charge, mult=mult, symbols=symbols,
                                    xyz=xyz_text)
                else:
                    s = seen[ik]
                    s['xtb_gap'] = min(s['xtb_gap'], gap)
                    s['proposers'] += ';' + disp
                    s['replicates'] += ';' + str(r['replicate'])
    rows = sorted(seen.values(), key=lambda r: r['xtb_gap'])

    manifest = []
    for r in rows:
        slug = r['inchikey'].split('-')[0]
        d = os.path.join(OUT, slug)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'geom.xyz'), 'w') as f:
            f.write(r['xyz'])
        with open(os.path.join(d, f'{slug}.inp'), 'w') as f:
            f.write(_INP_TEMPLATE.format(charge=r['charge'], mult=r['mult']))
        manifest.append(dict(inchikey=r['inchikey'], dir=slug,
                             canonical_smiles=r['smiles'],
                             proposers=r['proposers'], replicates=r['replicates'],
                             xtb_gap=f"{r['xtb_gap']:.3f}",
                             charge=r['charge'], mult=r['mult'],
                             n_atoms=len(r['symbols'])))
    with open(os.path.join(OUT, 'manifest.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(manifest[0].keys()))
        w.writeheader()
        w.writerows(manifest)
    print(f'{len(manifest)} unique molecules prepared in {OUT} (sorted by xTB gap ascending)')


def mol_from_smiles(smiles):
    """Return (inchikey, charge, mult, symbols, xyz_text) using the pipeline
    geometry; (None, 0, 1, [], '') if anything fails."""
    from rdkit import Chem
    from rdkit.Chem import inchi
    try:
        symbols, coords_bohr, charge = h._smiles_to_xyz_and_charge(smiles)
        mol = Chem.MolFromSmiles(smiles)
        ik = inchi.MolToInchiKey(mol)
        n_electrons = sum(_VALENCE.get(s, 0) for s in symbols) - charge
        mult = 1 if n_electrons % 2 == 0 else 2
        a = h._A_TO_BOHR
        lines = [str(len(symbols)),
                 f'{smiles}  (pipeline geometry: ETKDGv3 seed 0xf00d + MMFF)']
        for s, c in zip(symbols, coords_bohr):
            lines.append(f'{s} {c[0]/a:.6f} {c[1]/a:.6f} {c[2]/a:.6f}')
        return ik, charge, mult, symbols, '\n'.join(lines) + '\n'
    except Exception as err:
        print(f'  geometry/embed failure for {smiles}: {err}')
        return None, 0, 1, [], ''


if __name__ == '__main__':
    main()