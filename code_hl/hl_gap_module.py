from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import RDConfig
import sys, os
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
sys.path.append(os.path.join(RDConfig.RDContribDir, 'NP_Score'))
import sascorer, npscorer

import numpy as np
from tblite.interface import Calculator

scoring_args = [os.cpu_count(), 'GFN2-xTB']

# tblite's raw interface works in atomic units: positions in Bohr, energies in
# Hartree (see tblite.interface.Calculator docstring and tblite.ase's
# positions/Bohr conversion). RDKit conformers are in Angstrom.
_A_TO_BOHR = 1.889726125
_HARTREE_TO_EV = 27.211386


def _smiles_to_xyz_and_charge(smiles: str):
  '''RDKit 3D embed + MMFF optimise; returns (atom symbols, Bohr coords, formal charge).'''
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError(f'invalid SMILES: {smiles}')
  molH = Chem.AddHs(mol)
  params = AllChem.ETKDGv3()
  params.randomSeed = 0xf00d  # reproducible embeds across replicates
  if AllChem.EmbedMolecule(molH, params) != 0:
    raise ValueError(f'3D embed failed: {smiles}')
  AllChem.MMFFOptimizeMolecule(molH)
  conf = molH.GetConformer()
  symbols, coords = [], []
  for atom in molH.GetAtoms():
    symbols.append(atom.GetSymbol())
    pos = conf.GetAtomPosition(atom.GetIdx())
    coords.append([pos.x * _A_TO_BOHR, pos.y * _A_TO_BOHR, pos.z * _A_TO_BOHR])
  # Formal charge from the parsed molecule, NOT from scanning the SMILES string:
  # the original HL_gap_module counted '-' branch bonds as negative charges,
  # which misfires on druglike SMILES (e.g. a '-' branch bond not followed by 'c').
  charge = Chem.GetFormalCharge(molH)
  return symbols, np.array(coords), charge


def scoring_function(smiles: str):
  '''
    Receives a SMILES string and returns the HOMO-LUMO gap calculated with the
    GFN2-xTB tight-binding method (tblite). Same (score, aux) contract as
    docking_module.scoring_function: if the calculation fails, returns
    (100.0, None) -- the sentinel the proposer loop treats as a bad score.
  '''
  try:
    symbols, coords, charge = _smiles_to_xyz_and_charge(smiles)
    numbers = np.array([_PERIODIC_TABLE[s] for s in symbols])

    # Closed shell unless the electron count is odd (radical): then one unpaired
    # electron. tblite's spin knob is `uhf` (number of unpaired electrons).
    n_electrons = int(sum(_VALENCE_ELECTRONS.get(s, 0) for s in symbols)) - charge
    uhf = 0 if n_electrons % 2 == 0 else 1

    calc = Calculator(scoring_args[1], numbers, coords, charge=charge, uhf=uhf)
    res = calc.singlepoint()
    occ = res.get("orbital-occupations")
    orb = res.get("orbital-energies")  # Hartree
    homo = orb[occ > 0].max()
    lumo = orb[occ == 0].min()
    gap = (lumo - homo) * _HARTREE_TO_EV  # eV
    print(f'The HOMO-LUMO gap for {smiles} is: {gap:.3f} eV')
  except Exception as err:
    print(f'Could not calculate gap for {smiles}: {err}')
    return 100.0, None

  return float(gap), None


def calculate_HL_gap(smiles: str) -> str:
  '''
    Calculates the HOMO-LUMO gap of a single molecule with GFN2-xTB and returns it
    with the HOMO/LUMO orbital energies. Use on a molecule already deemed to have
    a small gap (score below ~4 eV) to characterise the frontier orbitals before
    proposing analogues around it.
  '''
  gap, _ = scoring_function(smiles)
  if gap == 100.0:
    return f'Gap calculation failed for {smiles} (invalid SMILES, 3D embed failure, or open-shell species that will not converge).'
  symbols, coords, charge = _smiles_to_xyz_and_charge(smiles)
  n_electrons = int(sum(_VALENCE_ELECTRONS.get(s, 0) for s in symbols)) - charge
  out = f'HOMO-LUMO gap for {smiles}: {gap:.3f} eV (GFN2-xTB; charge {charge}, {n_electrons} valence electrons, {len(symbols)} atoms).\n'
  out += ('The lower the gap, the better. A gap above ~5 eV is a poor candidate; '
          'below ~3 eV is a strong result. Compare with the other molecules in the list '
          'to learn which structural features shrink the gap.')
  return out


def calculate_SAS_and_NP(smiles_list: list):
  '''
  Calculate SAS and NP scores for a list of SMILES strings. SAS score is a measure
  of synthetic accessibility, and a value of 1 indicates that the molecule is easy to synthesize,
  while a value of 10 indicates that it is difficult to synthesize.
  NP score is a measure of natural product-likeness, and a higher score indicates that the
  molecule is more similar to natural products; the score runs from -5 to 5, with higher scores
  indicating greater similarity to natural products.

    Args:
        smiles_list (list): A list of SMILES strings representing the molecules to be scored.

    Returns:
        out_string (str): A string containing the SMILES, SAS score, and NP score for each molecule in the list.
  '''
  fscore = npscorer.readNPModel()

  out_string = '| SMILES | SAS Score | NP Score |\n'
  out_string += '|---------|-----------|----------|\n'
  for smiles in smiles_list:
      mol = Chem.MolFromSmiles(smiles)
      if mol is not None:
          sas_score = sascorer.calculateScore(mol)
          np_score = npscorer.scoreMol(mol, fscore)
          out_string += f'| {smiles} | {sas_score:.2f} | {np_score:.2f} |\n'
      else:
          out_string += f'| {smiles} | {"Invalid SMILES"} | {"Invalid SMILES"} |\n'
  return out_string


task_specific_prompt = '''# You are a materials science assistant. In the first user
message you will see a list of molecule SMILES strings and their corresponding HOMO-LUMO gaps.
Your task is to use the information in the list to learn trends about what makes a molecule
have a small or large HOMO-LUMO gap, and then use those trends to suggest new molecules
that should have the smallest possible HOMO-LUMO gap.
'''

task_specific_tools = '''
calculate_HL_gap(smiles: str) -> str: Calculates the HOMO-LUMO gap (GFN2-xTB) for a single
molecule and returns the gap in eV. Use on promising molecules to get their exact gap, and to
check whether a structural modification shrank or widened the gap.

calculate_SAS_and_NP(smiles_list: list) -> str: Calculates the SAS and NP scores for a list of
SMILES strings. SAS score is a measure of synthetic accessibility, and a value of 1 indicates
that the molecule is easy to synthesize, while a value of 10 indicates that it is difficult to
synthesize. NP score is a measure of natural product-likeness, and a higher score indicates
that the molecule is more similar to natural products; the score runs from -5 to 5, with higher
scores indicating greater similarity to natural products. Should be called for promising
molecules with small gaps and good synthetic accessibility (low SAS).
'''

# Element symbol -> nuclear charge. tblite takes atomic numbers; the interface module
# also exposes symbols_to_numbers(), but a tiny local map keeps scoring_function free
# of an import cycle and covers everything organic druglike/materials candidates need.
_PERIODIC_TABLE = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
    'K': 19, 'Ca': 20, 'Zn': 30, 'Se': 34, 'Br': 35, 'I': 53,
}

# Valence electron counts for the GFN2-xTB basis (light main-group elements).
_VALENCE_ELECTRONS = {
    'H': 1, 'He': 2, 'Li': 1, 'Be': 2, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8,
    'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Ar': 8,
    'K': 1, 'Ca': 2, 'Zn': 2, 'Se': 6, 'Br': 7, 'I': 7,
}

auxilliary_functions = [calculate_HL_gap, calculate_SAS_and_NP]