import sys
import os

# Add code folder to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'code'))

from docking_module import dock_and_get_interacting_residues

print('imported')

results = dock_and_get_interacting_residues('c1ccc(O)cc1')
print(results)
