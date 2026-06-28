from dockstring import load_target
from rdkit import Chem
from rdkit.Chem import RDConfig
import sys, os
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
sys.path.append(os.path.join(RDConfig.RDContribDir, 'NP_Score'))
import sascorer, npscorer
# ODDT 0.7 is unmaintained and calls np.in1d, which was removed in NumPy 2.x
# (deprecated in 2.0, dropped thereafter). np.isin is a drop-in replacement with
# identical signature and semantics, so alias it back before oddt is imported.
# This must run before `import oddt` so oddt sees the attribute on the numpy module.
import numpy as np
if not hasattr(np, "in1d"):
    np.in1d = np.isin

import oddt
from oddt.interactions import (close_contacts,
                               hbonds,
                               distance,
                               halogenbonds,
                               halogenbond_acceptor_halogen,
                               pi_stacking,
                               salt_bridges,
                               pi_cation,
                               hydrophobic_contacts)

import oddt.toolkits.ob as ob_toolkit  # noqa: F401  (kept for back-compat imports)

scoring_args = [os.cpu_count(),'DRD2']

def scoring_function(smiles: str):
  '''
    docks a molecule to the target and returns the docking score. If the docking fails, returns 0.0.
  '''
  target = load_target(scoring_args[1])
  try:
    score, aux = target.dock(smiles, scoring_args[0])
  except:
    score = 0.0
    aux = None
  return score, aux

task_specific_prompt = '''# You are a drug design assistant. In the first user message you will
see a list of molecule SMILES strings and docking scores.
The lower the docking score (the more negative), the more affinity the
molecule has for the protein in question. Your task is to use the information 
in the list to learn trends about what makes a molecule a good binder, and then 
use those trends to suggest new molecules that should have better docking scores 
(more negative) than the ones in the list.'''

task_specific_tools = '''
calculate_SAS_and_NP(smiles_list: list[str]) -> str: Calculates the SAS and NP scores for a list of SMILES strings. 
SAS score is a measure of synthetic accessibility, and a value of 1 indicates that the molecule is easy to synthesize, 
while a value of 10 indicates that it is difficult to synthesize. NP score is a measure of natural product-likeness, 
and a higher score indicates that the molecule is more similar to natural products; the score runs from -5 to 5, with higher scores 
indicating greater similarity to natural products. Should be called for promising molecules with good docking scores and
good Lipinski Properties.

dock_and_get_interacting_residues(smiles: str) -> str: Returns the docking score and types of interactions between a docked molecule and
residues in the target protein. To be used to evaluate whether or not a promising molecule is docked to the expected binding site.
Not to be used until a molecule unless a molecule has been deemed to have a low docking score and  good Lipinsky properties. 
For the system studied here, HMGCR, the interactions for the known binder Rosuvastatin are:
## Contacts between the ligand and proteins residues -------------------------
HBONDS interactions:
Atom O-25 forms a hbonds with GLU119
Atom O-25 forms a hbonds with GLU119
Atom O-28 forms a hbonds with ALA311
Atom O-28 forms a hbonds with ASP671
Atom F-14 forms a hbonds with ARG571
Atom F-14 forms a hbonds with SER642
Atom O-25 forms a hbonds with GLU119
Atom O-25 forms a hbonds with ASN315
Atom O-25 forms a hbonds with LYS672
Atom O-28 forms a hbonds with LYS673
Atom O-31 forms a hbonds with LYS295
HYDROPHOBIC_CONTACTS interactions:
Atom C-20 forms a hydrophobic_contacts with GLU119
Atom C-21 forms a hydrophobic_contacts with LEU122
Atom C-26 forms a hydrophobic_contacts with HIS312
Atom C-0 forms a hydrophobic_contacts with LEU413
Atom C-9 forms a hydrophobic_contacts with LEU413
Atom C-9 forms a hydrophobic_contacts with LEU417
Atom C-10 forms a hydrophobic_contacts with LEU417
Atom C-10 forms a hydrophobic_contacts with VAL664
PI_STACKING interactions:

PI_CATION interactions:

HALOGENBONDS interactions:
Atom F-14 forms a halogenbonds with SER642
SALT_BRIDGES interactions:
Atom O-32 forms a salt_bridges with LYS295
Atom O-32 forms a salt_bridges with ARG571
Atom O-32 forms a salt_bridges with LYS673
List of all interacting residues:
 ALA311, ARG571, ASN315, ASP671, GLU119, HIS312, LEU122, LEU413, LEU417, LYS295, LYS672, LYS673, SER642, VAL664 
'''

def calculate_SAS_and_NP(smiles_list: list[str]):
  '''
  Calculate SAS and NP scores for a list of SMILES strings. SAS score is a measure 
  of synthetic accessibility, and a value of 1 indicates that the molecule is easy to synthesize, 
  while a value of 10 indicates that it is difficult to synthesize. 
  NP score is a measure of natural product-likeness, and a higher score indicates that the 
  molecule is more similar to natural products; the score runs from -5 to 5, with higher scores 
  indicating greater similarity to natural products.

    Args:
        smiles_list (list[str]): A list of SMILES strings representing the molecules to be scored.

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

def dock_and_get_interacting_residues(smiles: str) -> str:
  '''
    Docks a molecule to the target and returns the score and interacting residues. If the docking fails, returns an empty list.

      Args:
          smiles (str): the SMILES string of the molecule to dock and get interacting residues.
      Returns:
          contacts_results (str): a string containing the types of interactions between the docked molecule and 
          residues in the target protein. If the docking fails, returns a string indicating that the docking
  '''
 
  score, aux = scoring_function(smiles)
  #print(f'Docking score: {score}')

  if aux is None:
    return "Docking failed. No interacting residues found."
  
  pose_mol = aux['ligand']
  pose_mol.SetProp('_Name',str(score))
  sdf_filename = "test_mol.sdf"
  w = Chem.SDWriter(sdf_filename)
  w.write(pose_mol)
  w.close()

  protein_file = 'HMGCR_dude_receptor_2.pdb'
  ligand_file = sdf_filename

  pro = next(oddt.toolkit.readfile('pdb',protein_file))
  lig = next(oddt.toolkit.readfile('sdf',ligand_file))
  pro.protein = True

  contacts_results = f'Docking score for the ligand: {score}\n'
  contacts_results += find_contacts(pro, lig)

  return contacts_results

def find_contacts(pro, lig) -> str:
  '''
    Finds the interactions between a docked molecule and residues in the target protein. 

      Args:
          pro: the protein object from oddt
          lig: the ligand object from oddt
      Returns:
          output_string (str): a string containing the types of interactions between the docked molecule and 
          residues in the target protein.
  '''

  int_types = ['hbonds', 'hydrophobic_contacts', 'pi_stacking', 'pi_cation', 'halogenbonds', 'salt_bridges']
  int_functions = [hbonds, hydrophobic_contacts, pi_stacking, pi_cation, halogenbonds, salt_bridges]

  output_string = '## Contacts between the ligand and proteins residues -------------------------'
  interacting_residues = []

  for int_type, int_function in zip(int_types, int_functions): 
    output_string += f'\n\n{int_type.upper()} interactions:\n'
    if int_type == 'pi_stacking':
      pro_atoms, lig_atoms, strict_parallel, strict_perpendicular = int_function(pro, lig)
    elif int_type == 'salt_bridges' or int_type == 'hydrophobic_contacts':
      pro_atoms, lig_atoms = int_function(pro, lig)
    elif int_type == 'halogenbonds' or int_type == 'pi_cation':
      pro_atoms, lig_atoms, strict = int_function(pro, lig, tolerance=30)
    else:
      pro_atoms, lig_atoms, strict = int_function(pro, lig)
  
    if int_type not in ['pi_stacking']:
      assert len(pro_atoms) == len(lig_atoms)
      atoms = []
      residues = []
      for atom, res in zip(lig_atoms, pro_atoms):
        atom_name = str(atom['atomtype']).split('.')[0]
        atom_number = str(atom['id'])
        res_name = str(res['resname'])
        res_number = str(res['resnum'])
        res_id = res_name+res_number
        output_string += f'Atom {atom_name}-{atom_number} forms a {int_type} with {res_id}\n'
        atoms.append(atom_number)
        residues.append(res_id)

      atoms = list(set(atoms))
      residues = list(set(residues))
      interacting_residues.append(residues)

    if int_type == 'pi_stacking':
      residues = []
      for res in pro_atoms:
        res_name = str(res['resname'])
        res_number = str(res['resnum'])
        res_id = res_name+res_number
        output_string += f'{res_id} forms pi-stacking interactions with the Ligand\n'
        residues.append(res_id)
        para = strict_parallel.sum()
        perp = strict_perpendicular.sum()

      residues = list(set(residues))
      interacting_residues.append(residues)

      try:
        output_string += f'The total number of parallel interactions are {para}\n'
        output_string += f'The total number of perpendicular interactions are {perp}\n'
      except:
        pass

  interacting_residues = [item for sublist in interacting_residues for item in sublist]
  interacting_residues = list(set(interacting_residues))
  interacting_residues.sort()

  output_string += f'List of all interacting residues:\n {', '.join(interacting_residues)}'

  return output_string


auxilliary_functions = [dock_and_get_interacting_residues, calculate_SAS_and_NP]
