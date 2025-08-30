from __future__ import annotations

from typing import Iterable, Tuple
import math
import mpmath as mp
from .contraction import infinite_sum_upper_bound, infinite_square_sum_upper_bound

def _validate_etas(etas: Iterable[float]) -> None:
    for e in etas:
        if not isinstance(e, (int, float)):
            raise TypeError("η_k must be real numbers.")
        if not math.isfinite(e):
            raise ValueError("η_k must be finite.")
        if e < 0:
            raise ValueError("η_k must be nonnegative.")

def admissible_collar_constant_upper(etas: Iterable[float]) -> float:
    """
    Strict upper bound on C for positivity of every factor (1 - C η_k) > 0:
    C < 1 / max_k η_k (if max η_k > 0). Returns +∞ if all η_k=0.
    """
    _validate_etas(etas)
    max_eta = max((float(e) for e in etas), default=0.0)
    if max_eta == 0.0:
        return math.inf
    return 1.0 / max_eta

def collar_product_logsum(
    etas: Iterable[float],
    C: float,
    *,
    mp_dps: int = 80,
) -> Tuple[float, bool]:
    """
    High-precision log-sum of Π_k (1 - C η_k). Returns (log_product, positive_flag).
    If any factor is ≤ 0, returns (-inf, False).
    """
    _validate_etas(etas)
    if not isinstance(C, (int, float)) or not math.isfinite(C) or C <= 0:
        raise ValueError("C must be a positive finite number.")
    mp.mp.dps = mp_dps
    log_sum = mp.mpf("0.0")
    for e in etas:
        term = 1.0 - C * float(e)
        if term <= 0.0:
            return float("-inf"), False
        log_sum += mp.log(term)
    return float(log_sum), True

def collar_product_bound(
    etas: Iterable[float],
    C: float,
    *,
    mp_dps: int = 80,
) -> float:
    """
    Rigorous lower bound for Π_k (1 - C η_k). Zero if any factor is ≤ 0.
    """
    log_prod, positive = collar_product_logsum(etas, C, mp_dps=mp_dps)
    if not positive:
        return 0.0
    val = mp.e**(mp.mpf(log_prod))
    return float(val if val > 0.0 else 0.0)

# -------- Analytic (infinite-steps) collar certificate --------

def collar_product_lower_bound_analytic(eta0: float, A: float, C: float) -> float:
    """
    Analytic, step-free lower bound using:
      log(1 - x) >= -x - x^2/(1 - x), for 0 <= x < 1.
    With x_k = C η_k, and bounds:
      Σ η_k <= (1/A) * z0/(1 - z0),  z0 = A η0  (requires z0<1),
      Σ η_k^2 <= η0 * Σ η_k.

    Returns exp( -C Ση_k - (C^2/(1 - C η0)) Ση_k^2 ), provided C η0 < 1 and z0 < 1.
    """
    if eta0 <= 0 or A <= 0 or C <= 0:
        raise ValueError("eta0, A, C must be positive.")
    z0 = A * eta0
    if not (0.0 < z0 < 1.0):
        raise ValueError("Requires A*eta0 < 1; got A*eta0 = %.6g" % z0)
    if not (C * eta0 < 1.0):
        # If this fails, the first factor (1 - C η0) is nonpositive; product lower bound is 0.
        return 0.0
    S1 = infinite_sum_upper_bound(eta0, A)
    S2 = infinite_square_sum_upper_bound(eta0, A)  # <= eta0 * S1
    # log lower bound:
    log_lb = -C * S1 - (C * C) * S2 / (1.0 - C * eta0)
    return float(math.exp(log_lb))