from __future__ import annotations
from typing import Tuple
import mpmath as mp

def dim_su3(p: int, q: int) -> int:
    """
    SU(3) irrep (p,q) dimension: d = (1/2)*(p+1)*(q+1)*(p+q+2)
    """
    if p < 0 or q < 0:
        raise ValueError("Dynkin labels (p,q) must be nonnegative.")
    return ((p + 1) * (q + 1) * (p + q + 2)) // 2

def casimir_su3(p: int, q: int) -> float:
    """
    Quadratic Casimir with C2(fundamental)=4/3:
    C2(p,q) = (p^2 + q^2 + p q)/3 + (p + q)
    """
    return (p*p + q*q + p*q)/3.0 + (p + q)

def _shell_dim_envelope(n: int) -> float:
    """
    Safe polynomial envelope for total dimension on shell p+q=n:
      sum_{p+q=n} d_{p,q} <= K * (n+2)^4, with K=1/2 conservative.
    """
    return 0.5 * (n + 2)**4

def _c2_shell_min(n: int) -> float:
    """
    Exact minimum of C2(p,q) over p,q >= 0 with p+q=n:
      min C2 = n^2/4 + n  (attained near p=q=n/2).
    """
    return (n*n)/4.0 + n

def tail_bound_su3_sharp(beta: float, p0: int, q0: int, *, mp_dps: int = 80) -> Tuple[mp.mpf, mp.mpf]:
    """
    Upper bound on the SU(3) character/Casimir tail:
      S_tail = sum_{p>=p0 or q>=q0} d_{p,q} * exp(-beta * C2(p,q) / 6)

    Strategy:
      • Let N = max(p0, q0). Split:
          (i) exact partial over triangle p+q < N with (p>=p0 or q>=q0);
         (ii) shell tail for n >= N via sharp shell minimum C2_min(n) = n^2/4 + n
             and envelope K (n+2)^4. Integrate continuous relaxation n->x from N-1 to ∞.

    mp_dps controls mpmath precision.
    """
    if beta <= 0:
        raise ValueError("beta must be positive.")
    if p0 < 0 or q0 < 0:
        raise ValueError("p0,q0 must be nonnegative.")

    mp.mp.dps = mp_dps
    N = max(p0, q0)
    partial = mp.mpf("0.0")

    # exact partial on p+q < N with (p>=p0 or q>=q0)
    for p in range(0, N):
        for q in range(0, N - p):
            if p >= p0 or q >= q0:
                d = dim_su3(p, q)
                c = casimir_su3(p, q)
                partial += mp.e**(-beta * c / 6.0)

    # tail integral bound
    alpha = mp.mpf(beta) / 24.0  # from n^2/4, divided by 6
    gamma = mp.mpf(beta) / 6.0   # from +n, divided by 6
    K = mp.mpf("0.5")

    def envelope(x):
        return K * (x + 2.0)**4 * mp.e**(-alpha * x * x - gamma * x)

    tail = mp.quad(envelope, [N - 1, mp.inf])
    return partial, tail

# Back-compat alias (older imports)
tail_bound_arxiv_style = tail_bound_su3_sharp