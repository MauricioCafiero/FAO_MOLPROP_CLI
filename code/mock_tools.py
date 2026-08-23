"""mock_tools.py - synthetic stand-in for docking_module.scoring_function, used by
molopt.py/molopt_oa.py's --mock-tools flag to smoke-test the LLM tool-calling loop
without paying for real Vina docking.

scoring_function is the single choke point all three docking-dependent tools call
through: MolPropOp.grow_cycle / MolPropOp.replace_groups (imported their own
`scoring_function` binding at import time) and docking_module.dock_and_get_interacting_residues
(uses its own module-level binding). Patching both modules' `scoring_function` name
covers all three without touching their bodies. The other tools (make_random_list,
related, lipinski, calculate_SAS_and_NP) don't dock, so they're left alone.

aux=None makes dock_and_get_interacting_residues take its existing "Docking failed.
No interacting residues found." branch -- a real docked pose can't be synthesized
without actually docking, so this is the honest stand-in rather than fabricated
contacts.
"""

import random


def mock_scoring_function(smiles: str):
    """Instant synthetic (score, aux) in place of a real dockstring/Vina call."""
    return round(random.uniform(-9.5, -5.5), 2), None


def install(mock: bool) -> None:
    """If mock, monkeypatch scoring_function in both modules that call it."""
    if not mock:
        return
    import MolPropOp
    import docking_module
    MolPropOp.scoring_function = mock_scoring_function
    docking_module.scoring_function = mock_scoring_function
