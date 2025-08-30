from __future__ import annotations
from typing import Tuple
import math

def spectral_gap_from_tube_cost(tau0: float) -> Tuple[float, float]:
    """
    If every time-slice tube has cost >= tau0 > 0, the transfer operator T
    has spectrum contained in {1} ∪ [exp(-tau0), 1).
    The Hamiltonian gap m0 satisfies m0 >= tau0.
    Returns (lambda_max_below_1, m0_lower_bound).
    """
    if tau0 <= 0:
        raise ValueError("tau0 must be positive.")
    lam = math.exp(-tau0)
    return lam, tau0