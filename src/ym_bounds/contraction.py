from __future__ import annotations
from typing import List, Tuple
import math

def contraction_sequence(eta0: float, A: float, steps: int = 20) -> List[float]:
    """
    Generate η_k with quadratic contraction: η_{k+1} <= A * η_k^2.
    We saturate the inequality (worst case) to compute rigorous upper bounds.
    """
    if eta0 <= 0 or A <= 0:
        raise ValueError("eta0 and A must be positive.")
    if steps <= 0:
        raise ValueError("steps must be positive.")
    seq = [eta0]
    for _ in range(steps):
        seq.append(A * seq[-1] * seq[-1])
    return seq

def check_summability_and_product(etas: List[float], C: float) -> Tuple[float, float, bool, bool]:
    """
    Return (sum_eta, prod_term, sum_finite?, prod_positive?)
    where prod_term = Π_k (1 - C * η_k) (clipped at 0 for numerics but flag indicates theoretical).
    """
    if C <= 0:
        raise ValueError("C must be positive.")
    s = 0.0
    prod = 1.0
    prod_pos = True
    for e in etas:
        if e < 0 or not math.isfinite(e):
            raise ValueError("All η_k must be finite and nonnegative.")
        s += e
        term = 1.0 - C * e
        if term <= 0.0:
            prod_pos = False
            term = 0.0
        prod *= term
    return s, prod, math.isfinite(s), prod_pos

# -------- Analytic, infinite-steps certificates (no dependence on 'steps') --------

def infinite_sum_upper_bound(eta0: float, A: float) -> float:
    """
    Analytic bound for Σ_{k≥0} η_k assuming η_{k+1} <= A η_k^2 and z0 := A η0 < 1.
    Let z_k := A η_k, then z_{k+1} <= z_k^2 ⇒ z_k <= z0^(2^k).
    Hence Σ η_k <= (1/A) Σ z0^(2^k) <= (1/A) * z0/(1 - z0) for 0<z0<1.
    """
    if eta0 <= 0 or A <= 0:
        raise ValueError("eta0 and A must be positive.")
    z0 = A * eta0
    if not (0.0 < z0 < 1.0):
        raise ValueError("Requires A*eta0 < 1 for contraction; got A*eta0 = %.6g" % z0)
    return (1.0 / A) * (z0 / (1.0 - z0))

def infinite_square_sum_upper_bound(eta0: float, A: float) -> float:
    """
    Trivial but useful bound: Σ η_k^2 <= η0 * Σ η_k since η_k is nonincreasing.
    Combined with infinite_sum_upper_bound gives a fully explicit constant.
    """
    S1 = infinite_sum_upper_bound(eta0, A)
    return eta0 * S1