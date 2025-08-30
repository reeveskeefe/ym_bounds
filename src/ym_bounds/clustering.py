from __future__ import annotations
from typing import Callable
import math

def exp_clustering_bound(distance: float, m_star: float, prefactor: float = 1.0) -> float:
    """
    |<F(x)F(0)>| <= prefactor * exp(-m_star * |x|)
    """
    if distance < 0 or m_star <= 0 or prefactor <= 0:
        raise ValueError("distance >= 0, m_star > 0, prefactor > 0 required.")
    return prefactor * math.exp(-m_star * distance)