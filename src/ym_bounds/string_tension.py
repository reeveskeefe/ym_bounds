from __future__ import annotations
from typing import Tuple
import math

def sigma_phys_from_lattice(sigma_lat: float, a: float) -> float:
    """
    Convert lattice string tension (in lattice units) to physical units:
    sigma_phys = sigma_lat / a^2.
    """
    if sigma_lat <= 0 or a <= 0:
        raise ValueError("sigma_lat and a must be positive.")
    return sigma_lat / (a * a)

def area_law_upper_bound(area: float, sigma_phys: float) -> float:
    """
    Upper bound on <W(C)> via area law: exp(-sigma_phys * Area(C)).
    """
    if area < 0 or sigma_phys <= 0:
        raise ValueError("area >= 0, sigma_phys > 0 required.")
    return math.exp(-sigma_phys * area)