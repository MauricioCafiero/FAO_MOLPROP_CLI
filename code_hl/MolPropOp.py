from rdkit import Chem
from rdkit.Chem import Draw, QED
import os, re, random, sys
from typing import Optional
import pubchempy as pcp

from hl_gap_module import scoring_function, scoring_args

''' define functions'''
# define a scoring function that can be substituted into the cycle functions.
# This allows us to easily switch between gap calculations and other scoring functions

base_rings = ['c1ccccc1', #benzene
               'n1ccccc1', #pyridine
               'o1cccc1',  #furan
               's1cccc1',  #thiophene
               '[nH]1cccc1', #pyrrole
               'n1c[nH]cc1', #imidazole
               'c1ccc2ccccc2c1', #naphthalene
               'c1ccc2cc3ccccc3cc2c1', #anthracene
               'O=c1cc(-c2ccccc2)oc2ccccc12'] #flavone

free_carbon_search_patterns = ['c[0-9]c', r'1c\[n', 'cc', r'c\[nH\]']
free_carbon_insert_points = [2, 2, 1, 1]

clean_ring_locations = [
  [1], #benzene
  [2, 3, 4], #pyridine
  [2, 3], #furan
  [2, 3], #thiophene
  [5, 6], #pyrrole
  [2, 7, 9], #imidazole
  [1, 3], #naphthalene
  [1, 3, 6],  #anthracene
  [4, 10, 11, 12, 13, 21, 22, 23] #flavone
]

e_withdraw = [
    'I',              #iodo
    'Br',             #bromo 
    'Cl',             #chloro
    'F',              #fluoro
    'C(Cl)(Cl)(Cl)',  #trichloromethyl
    'C(F)(F)(F)',     #trifluoromethyl
    'C(=O)[O-]',      #carboxylate
    'C(=O)',          #carbonyl
    'C#N',            #nitrile
    '[N+](=O)[O-]',   #nitro
    '[NH3+]']         #ammonium

e_donate = [
    "[O-]",          # Phenoxide (strongest resonance donor)
    "N(C)C",         # Dimethylamino (-NMe2)
    "NC",            # Methylamino (-NHMe)
    "N",             # Amino (-NH2)
    "O",             # Hydroxy (-OH)
    "OC",            # Methoxy (-OMe)
    "NC(=O)C",       # Acetamido (-NHCOCH3)
    "SC",            # Methylthio (-SMe)
    "OC(=O)C",       # Acetoxy (-OCOCH3)
    "C(C)(C)C",      # tert-Butyl (-C(CH3)3)
    "C(C)C",         # Isopropyl (-CHMe2)
    "CC",            # Ethyl (-Et)
    "C",             # Methyl (-Me)
    "c5ccccc5",      # Phenyl (-Ph)
    "C=C"           # Vinyl (-CH=CH2)
    #"[Si](C)(C)C"    # Trimethylsilyl (-SiMe3)
]

linkers = [
    "C",             # Methylene (-CH2-)
    "CC",            # Ethylene (-CH2CH2-)
    "CCC",           # Propylene (-CH2CH2CH2-)
    "C=C",           # Vinylene (-CH=CH-)
    "C#C",           # Acetylene (-C≡C-)
    "CC=C",          # Allylene (-CH2CH=CH-)
    "C=CC",          # Propenylene (-CH=CHCH2-)
    "O",             # Oxygen (-O-)
    "S",             # Sulfur (-S-)
    "N",             # Nitrogen (-NH-)
    "C(=O)",          # Carbonyl (-C(=O)-)
    "C(=O)O",        # Ester (-C(=O)O-)
    "C(=O)N",        # Amide (-C(=O)N-)
    "C(=O)C"       # Ketone (-C(=O)C-)
]

withdraw_with_linkers = [f'{linker}({e})' for linker in linkers for e in e_withdraw] 
donate_with_linkers = [f'{linker}({e})' for linker in linkers for e in e_donate]

def make_random_list(num_items: int) -> tuple[list, list]:
    '''
    selects num_items from the lists e_withdraw, e_donate, withdraw_with_linkers, and 
    donate_with_linkers and returns them as a single list along with the remaining items.
    Args:
      num_items: the number of items to select from the lists

    Returns:
      tuple: (selected_items: list, remaining_items: list)
        - selected_items: a list of num_items items randomly selected from the lists
        - remaining_items: a list of items that were not selected
    '''
    # Combine all substituent lists
    all_substituents = e_withdraw + e_donate + withdraw_with_linkers + donate_with_linkers
    
    # Randomly select num_items (without replacement if num_items <= total available)
    if num_items >= len(all_substituents):
        return all_substituents, []
    else:
        selected = random.sample(all_substituents, num_items)
        remaining = [item for item in all_substituents if item not in selected]
        return selected, remaining

def add_ring(ring_smiles: str, locations: list, comment: str = ''):
    '''
    Add a new ring to both base_rings and clean_ring_locations lists.
    
    Args:
      ring_smiles: SMILES string for the ring structure
      locations: list of integers representing clean substitution locations
      comment: optional comment describing the ring (e.g., 'benzene', 'pyridine')
    
    Returns:
      tuple: (success: bool, message: str)
    '''
    # Validate the SMILES string
    mol = Chem.MolFromSmiles(ring_smiles)
    if mol is None:
        return False, f"Invalid SMILES string: {ring_smiles}"
    
    # Check if ring already exists
    if ring_smiles in base_rings:
        return False, f"Ring {ring_smiles} already exists in base_rings"
    
    # Add to both lists
    base_rings.append(ring_smiles)
    clean_ring_locations.append(locations)
    
    message = f"Successfully added ring: {ring_smiles}"
    if comment:
        message += f" ({comment})"
    message += f" with locations {locations}"
    
    return True, message

def remove_ring(ring_smiles: Optional[str] = None, index: Optional[int] = None):
    '''
    Remove a ring from both base_rings and clean_ring_locations lists.
    Can remove by SMILES string or by index.
    
    Args:
      ring_smiles: SMILES string of the ring to remove
      index: index of the ring to remove (0-based)
    
    Returns:
      tuple: (success: bool, message: str)
    '''
    if ring_smiles is None and index is None:
        return False, "Must provide either ring_smiles or index"
    
    if ring_smiles is not None and index is not None:
        return False, "Provide only one of ring_smiles or index, not both"
    
    try:
        idx: int
        if ring_smiles is not None:
            # Find by SMILES string
            if ring_smiles not in base_rings:
                return False, f"Ring {ring_smiles} not found in base_rings"
            idx = base_rings.index(ring_smiles)
        else:
            # Use provided index
            assert index is not None  # Type checker satisfaction
            if index < 0 or index >= len(base_rings):
                return False, f"Index {index} out of range (0-{len(base_rings)-1})"
            idx = index
        
        # Remove from both lists
        removed_ring = base_rings.pop(idx)
        removed_locs = clean_ring_locations.pop(idx)
        
        return True, f"Successfully removed ring: {removed_ring} with locations {removed_locs}"
    
    except Exception as e:
        return False, f"Error removing ring: {str(e)}"

def list_rings():
    '''
    Display all rings in base_rings with their corresponding clean_ring_locations.
    
    Returns:
      str: formatted string showing all rings and their locations
    '''
    output = "Current rings in base_rings:\n"
    output += "=" * 80 + "\n"
    for i, (ring, locs) in enumerate(zip(base_rings, clean_ring_locations)):
        output += f"{i}: {ring:35s} -> {locs}\n"
    output += "=" * 80
    return output

def get_name(var):
    name = [k for k, v in globals().items() if v is var][0]
    return name

def convert_smiles(smiles: str) -> str:
  '''
    converts a SMILES string to a canonical form. This is useful for ensuring that
    the substituents are added to the correct positions on the ring.
    Args:
      smiles: a SMILES string for a molecule

    Returns:
      new_smiles: a canonical SMILES string for the same molecule
  '''
  new_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
  return new_smiles

def sub_cycle(substituents: list = e_withdraw, scoring_args: list = scoring_args):
  '''
  add substituents to free carbons in the molecule. Only adds to symmetrically unique 
  positions on the ring, as defined in clean_ring_locations.
  
  Args:
    substituents : a list of SMILES strings for substituents to add

  Returns:
    total_list : a list of tuples containing the new SMILES strings and their corresponding scores
  '''
  print('=============================================================================')
  print(f"Starting sub cycle on base rings for method {scoring_args[1]}.")

  best_score = 0.0
  best_smiles = ''
  total_list = []
  for ring, locs in zip(base_rings, clean_ring_locations):
    for loc in locs:

      for e in substituents:

        new_smiles = f'{ring[:loc+1]}({e}){ring[loc+1:]}'
        new_mol = Chem.MolFromSmiles(new_smiles)
        if new_mol != None:
          try:
              score, aux = scoring_function(new_smiles)
              total_list.append((new_smiles, score))
          except:
              print(f"Error scoring {new_smiles}")
              score = 0.0
              total_list.append((new_smiles, 'scoring error, skipped'))
          print(f"{new_smiles}: {score}")
          if score < best_score:
            best_score = score
            best_smiles = new_smiles    
            print(f"=========== New best score: {best_score} for {new_smiles} ======")
        else:
          print(new_smiles, 'bad')
          total_list.append((new_smiles, 'invalid SMILES, skipped'))

  print(f"Best score: {best_score} for {best_smiles}")
  print('=============================================================================')

  return total_list

def grow_cycle(best_smiles: str = 'c1ccccc1', best_score: float = 0.0, substituents: list[str] = e_withdraw):
    '''
    add substituents to free carbons in the molecule. 

    Args:
      best_smiles : the current best molecule, as a SMILES string
      best_score : the current best HOMO-LUMO gap, as a float
      substituents : a list of SMILES strings for substituents to add

    Returns:
      total_list : a list of tuples containing the new SMILES strings and their corresponding scores
    '''
    print('=============================================================================')
    print(f"Starting grow cycle with best score {best_score} for {best_smiles}.")
  
    total_list = []
    for pattern, insert_point in zip(free_carbon_search_patterns, free_carbon_insert_points):
        current_smiles = best_smiles
        for match in re.finditer(f'(?={pattern})', current_smiles):
            for e in substituents:
                new_smiles = f'{current_smiles[:match.start() + insert_point]}({e}){current_smiles[match.start() + insert_point:]}'
                mol = Chem.MolFromSmiles(new_smiles)
                if mol != None:
                    #print(new_smiles)
                    try:
                        score, aux = scoring_function(new_smiles)
                        total_list.append((new_smiles, score))
                    except:
                        print(f"Error scoring {new_smiles}")
                        score = 0.0
                        total_list.append((new_smiles, 'scoring error, skipped'))
                    print(f"{new_smiles}: {score}")
                    if score < best_score:
                        best_score = score
                        best_smiles = new_smiles
                        print(f"=========== New best score: {best_score} for {best_smiles} ========")
                else:
                    print(new_smiles, 'bad')
                    total_list.append((new_smiles, 'invalid SMILES, skipped'))
    print(f"Best score: {best_score} for {best_smiles}")
    print('=============================================================================')
    return total_list

def replace_groups(orig_smiles: str = 'c1ccccc1', best_score: float = 0.0, substituents_to_replace: list[str] = e_withdraw, 
                   new_substituents: list[str] = e_donate):
    '''
    replace existing substituents in the molecule with new ones. 

    Args:
      orig_smiles: the current best molecule, as a SMILES string
      best_score: the current best HOMO-LUMO gap, as a float
      substituents_to_replace: a list of SMILES strings for substituents to replace
      new_substituents: a list of SMILES strings for substituents to add

    Returns:
      total_list: a list of tuples containing the new SMILES strings and their corresponding scores
    '''
    print('=============================================================================')
    print(f"Starting replace cycle with best score {best_score} for {orig_smiles}.")
    best_smiles = orig_smiles
    #look in best_smiles for substituents in the substituents_to_replace list and replace them with substituents in the new_substituents list
    total_list = []
    for old in substituents_to_replace:
        if old in orig_smiles:
            print(f"Found {old} in {orig_smiles}")
            for new in new_substituents:
                new_smiles = orig_smiles.replace(old, new)
                mol = Chem.MolFromSmiles(new_smiles)
                if mol != None:
                    try:
                        score, aux = scoring_function(new_smiles)
                        total_list.append((new_smiles, score))
                    except:
                        print(f"Error scoring {new_smiles}")
                        score = 0.0
                        total_list.append((new_smiles, 'scoring error, skipped'))
                    #print(f"{new_smiles}: {score}")
                    if score < best_score:
                        best_score = score
                        best_smiles = new_smiles
                        print(f"=========== New best score: {best_score} for {best_smiles} ========")
                else:
                    print(new_smiles, 'bad')
                    total_list.append((new_smiles, 'invalid SMILES, skipped'))

    print(f"Best score: {best_score} for {best_smiles}")
    print('=============================================================================')
    return total_list

def related(smiles_list: list[str]):
  '''
    Queries Pubchem for similar molecules based on the smiles string or name

      Args:
        smiles_list: the input smiles strings

      Returns:
        related_string: a string containing the related molecules and their properties, 
        formatted for display in the UI
  '''
  print("related tool")
  print('===================================================')

  related_string = ''
  for smiles in smiles_list:
    try:
        res = pcp.get_compounds(smiles, "smiles", searchtype="similarity",listkey_count=25)
        related_string += f'## The following molecules are similar to {smiles}: \n'
        print('got related molecules with smiles')

        sub_smiles = []

        i = 0
        for compound in res:
            if i == 0:
                print(compound.iupac_name)
                i+=1
            sub_smiles.append(compound.smiles)
            related_string += f'- Name: {compound.iupac_name}\n'
            related_string += f'- SMILES: {compound.smiles}\n'
            related_string += f'- Molecular Weight: {compound.molecular_weight}\n'
            related_string += f'- LogP: {compound.xlogp}\n'
            related_string += '===================\n'
    except:
        related_string += f'{smiles}: Fail\n'

  return related_string

def lipinski(smiles_list: list[str] = ['c1ccccc1']):
  '''
    A tool to calculate QED and other lipinski properties of a molecule.

      Args:
        smiles_list: the input smiles strings

      Returns:
        total_lipinski_string: a string of the QED and other lipinski properties of the molecules,
                      including Molecular Weight, LogP, HBA, HBD, Polar Surface Area,
                      Rotatable Bonds, Aromatic Rings and Undesireable Moieties.
  '''
  print("lipinski tool")
  print('===================================================')

  total_lipinski_string = ''

  for smiles in smiles_list:
    for ion in ['.[Na+]', '.[K+]', '.[Cl-]', '.[Br-]', '[Na+].', '[K+].', '[Cl-].', '[Br-].']:
        smiles = smiles.replace(ion, '')
    lipinski_list = []
    try:
        mol = Chem.MolFromSmiles(smiles)
        qed = Chem.QED.default(mol)

        p = Chem.QED.properties(mol)
        mw = p[0]
        logP = p[1]
        hba = p[2]
        hbd = p[3]
        psa = p[4]
        rb = p[5]
        ar = p[6]
        um = p[7]

        lipinski_list.append(qed)
        lipinski_list.append(mw)
        lipinski_list.append(logP)
        lipinski_list.append(hba)
        lipinski_list.append(hbd)
        lipinski_list.append(psa)
        lipinski_list.append(rb)
        lipinski_list.append(ar)
        lipinski_list.append(um)

        total_lipinski_string += f"Properties of SMILES: {smiles}: QED: {qed:.3f}\n"
        total_lipinski_string += f"Molecular Weight: {mw:.3f}, LogP: {logP:.3f}\n"
        total_lipinski_string += f"Hydrogen bond acceptors: {hba}, Hydrogen bond donors: {hbd}\n" 
        total_lipinski_string += f"Polar Surface Area: {psa:.3f}, Rotatable Bonds: {rb}\n" 
        total_lipinski_string += f"Aromatic Rings: {ar}, Undesireable moieties: {um}\n"
        total_lipinski_string += "===================================================\n"    
    except:
        total_lipinski_string += f"SMILES: {smiles}, Could not get properties\n"
  
  return total_lipinski_string

''' batch versions of the above functions. These generate all SMILES first, then score them in batch with a single function call to scoring_function. This is more efficient than scoring each SMILES individually, especially when using gap calculations as the scoring function.'''

def sub_cycle_batch(substituents: list = e_withdraw, scoring_args: list = scoring_args):
    '''
    Batch version of sub_cycle. Add substituents to free carbons in the molecule. 
    Only adds to symmetrically unique positions on the ring, as defined in clean_ring_locations.
    Generates all SMILES first, then scores them in batch.
    
    Args:
      substituents: a list of SMILES strings for substituents to add
      scoring_args: arguments to pass to the scoring function

    Returns:
      total_list: a list of tuples containing the new SMILES strings and their corresponding scores
    '''
    print('=============================================================================')
    print(f"Starting batch sub cycle on base rings for method {scoring_args[1]}.")

    # First pass: generate all valid SMILES
    intermediate_list = []
    for ring, locs in zip(base_rings, clean_ring_locations):
        for loc in locs:
            for e in substituents:
                new_smiles = f'{ring[:loc+1]}({e}){ring[loc+1:]}'
                new_mol = Chem.MolFromSmiles(new_smiles)
                if new_mol != None:
                    intermediate_list.append(new_smiles)
                else:
                    print(new_smiles, 'bad')

    print(f"Generated {len(intermediate_list)} valid SMILES strings. Scoring in batch...")

    # Second pass: score all SMILES in batch with single function call
    try:
        scores = scoring_function(intermediate_list)
    except:
        print(f"Error in batch scoring, returning zero scores")
        scores = [0.0] * len(intermediate_list)

    # Combine SMILES with scores
    total_list = list(zip(intermediate_list, scores))

    # Find best score
    best_score = 0.0
    best_smiles = ''
    for smiles, score in total_list:
        print(f"{smiles}: {score}")
        if score < best_score:
            best_score = score
            best_smiles = smiles
            print(f"=========== New best score: {best_score} for {best_smiles} ======")

    print(f"Best score: {best_score} for {best_smiles}")
    print('=============================================================================')

    return total_list

def grow_cycle_batch(best_smiles: str = 'c1ccccc1', best_score: float = 0.0, substituents: list = e_withdraw):
    '''
    Batch version of grow_cycle. Add substituents to free carbons in the molecule.
    Generates all SMILES first, then scores them in batch.

    Args:
      best_smiles: the current best molecule, as a SMILES string
      best_score: the current best HOMO-LUMO gap, as a float
      substituents: a list of SMILES strings for substituents to add

    Returns:      total_list: a list of tuples containing the new SMILES strings and their corresponding scores
    '''
    print('=============================================================================')
    print(f"Starting batch grow cycle with best score {best_score} for {best_smiles}.")
    
    # First pass: generate all valid SMILES
    intermediate_list = []
    for pattern, insert_point in zip(free_carbon_search_patterns, free_carbon_insert_points):
        current_smiles = best_smiles
        for match in re.finditer(f'(?={pattern})', current_smiles):
            for e in substituents:
                new_smiles = f'{current_smiles[:match.start() + insert_point]}({e}){current_smiles[match.start() + insert_point:]}'
                mol = Chem.MolFromSmiles(new_smiles)
                if mol != None:
                    intermediate_list.append(new_smiles)
                else:
                    print(new_smiles, 'bad')

    print(f"Generated {len(intermediate_list)} valid SMILES strings. Scoring in batch...")

    # Second pass: score all SMILES in batch with single function call
    try:
        scores = scoring_function(intermediate_list)
    except:
        print(f"Error in batch scoring, returning zero scores")
        scores = [0.0] * len(intermediate_list)

    # Combine SMILES with scores
    total_list = list(zip(intermediate_list, scores))

    # Find best score
    for smiles, score in total_list:
        print(f"{smiles}: {score}")
        if score < best_score:
            best_score = score
            best_smiles = smiles
            print(f"=========== New best score: {best_score} for {best_smiles} ========")

    print(f"Best score: {best_score} for {best_smiles}")
    print('=============================================================================')
    
    return total_list

def replace_groups_batch(orig_smiles: str = 'c1ccccc1', best_score: float = 0.0, 
                         substituents_to_replace: list = e_withdraw, 
                         new_substituents: list = e_donate):
    '''
    Batch version of replace_groups. Replace existing substituents in the molecule with new ones.
    Generates all SMILES first, then scores them in batch.

    Args:
      orig_smiles: the current best molecule, as a SMILES string
      best_score: the current best HOMO-LUMO gap, as a float
      substituents_to_replace: a list of SMILES strings for substituents to replace
      new_substituents: a list of SMILES strings for substituents to add

    Returns:      total_list: a list of tuples containing the new SMILES strings and their corresponding scores
    '''
    print('=============================================================================')
    print(f"Starting batch replace cycle with best score {best_score} for {orig_smiles}.")
    best_smiles = orig_smiles
    
    # First pass: generate all valid SMILES
    intermediate_list = []
    for old in substituents_to_replace:
        if old in orig_smiles:
            print(f"Found {old} in {orig_smiles}")
            for new in new_substituents:
                new_smiles = orig_smiles.replace(old, new)
                mol = Chem.MolFromSmiles(new_smiles)
                if mol != None:
                    intermediate_list.append(new_smiles)
                else:
                    print(new_smiles, 'bad')

    print(f"Generated {len(intermediate_list)} valid SMILES strings. Scoring in batch...")

    # Second pass: score all SMILES in batch with single function call
    try:
        scores = scoring_function(intermediate_list)
    except:
        print(f"Error in batch scoring, returning zero scores")
        scores = [0.0] * len(intermediate_list)

    # Combine SMILES with scores
    total_list = list(zip(intermediate_list, scores))

    # Find best score
    for smiles, score in total_list:
        if score < best_score:
            best_score = score
            best_smiles = smiles
            print(f"=========== New best score: {best_score} for {best_smiles} ========")

    print(f"Best score: {best_score} for {best_smiles}")
    print('=============================================================================')
    
    return total_list
