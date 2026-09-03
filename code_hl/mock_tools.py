"""mock_tools.py - synthetic stand-in for hl_gap_module.scoring_function, used by
molopt.py/molopt_oa.py's --mock-tools flag to smoke-test the LLM tool-calling loop
without paying for real GFN2-xTB gap calculations.

scoring_function is the single choke point all three gap-dependent tools call
through: MolPropOp.grow_cycle / MolPropOp.replace_groups (imported their own
`scoring_function` binding at import time) and hl_gap_module.calculate_HL_gap
(uses its own module-level binding). Patching both modules' `scoring_function` name
covers all three without touching their bodies. The other tools (make_random_list,
related, calculate_SAS_and_NP) don't score, so they're left alone.

aux=None makes calculate_HL_gap take its existing "calculation failed" branch --
a real frontier-orbital analysis can't be synthesized, so this is the honest
stand-in rather than fabricated orbital energies.
"""

import random


def mock_scoring_function(smiles: str):
    """Instant synthetic (score, aux) in place of a real GFN2-xTB call.

    Gaps are minimized, so plausible synthetic values sit in the 1.5-6 eV range
    (100.0 is the failure sentinel, never mocked)."""
    return round(random.uniform(1.5, 6.0), 2), None


def install(mock: bool) -> None:
    """If mock, monkeypatch scoring_function in both modules that call it."""
    if not mock:
        return
    import MolPropOp
    import hl_gap_module
    MolPropOp.scoring_function = mock_scoring_function
    hl_gap_module.scoring_function = mock_scoring_function