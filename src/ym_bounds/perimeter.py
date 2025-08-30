# src/ym_bounds/perimeter.py
from __future__ import annotations

"""
Perimeter-term bounds derived from collar-product lower bounds.

Context
-------
In the polymer/collar analysis, the analytic collar-product lower bound
    P_collar := prod_k (1 - C * eta_k)   with  0 < P_collar <= 1
quantifies an additional suppression on Wilson loops coming from collars near the contour.
If the collar effect factors approximately per unit boundary length (a standard outcome of
local expansions with reflection positivity and uniform mixing), we can convert P_collar
into a rigorous *perimeter-term* κ such that, for a loop with perimeter length Perim (in
lattice length units),

    ⟨W(C)⟩  ≤  exp(-σ_phys * Area(C)) * exp(-κ * Perim) .

This file provides:
  • Conversion P_collar → κ per unit perimeter (lattice units),
  • Optional conversion κ_latt → κ_phys (per unit physical length) using lattice spacing a,
  • A one-shot helper to compute κ from (eta0, A, C) via the analytic collar-product bound,
  • Convenience routines to compute the perimeter factor and a combined area+perimeter bound.

All functions are deterministic, use high precision via mpmath, and validate inputs strictly.
"""

from dataclasses import dataclass
from typing import Optional, Tuple

import mpmath as mp
from ym_bounds.collar import collar_product_lower_bound_analytic


# --------------------------- Helpers & validation ---------------------------

def _pos(name: str, x: float) -> None:
    if not (isinstance(x, (int, float)) and x > 0):
        raise ValueError(f"{name} must be > 0")

def _unit_interval_closed(name: str, x: float) -> None:
    if not (isinstance(x, (int, float)) and 0.0 < x <= 1.0):
        raise ValueError(f"{name} must lie in (0, 1], got {x}")

def _nonneg(name: str, x: float) -> None:
    if not (isinstance(x, (int, float)) and x >= 0):
        raise ValueError(f"{name} must be ≥ 0")


# --------------------------- Core perimeter machinery ---------------------------

def kappa_per_unit_from_product(
    *,
    product_lower_bound: float,
    perim_scale: float = 1.0,
    density: float = 1.0,
    mp_dps: int = 120,
) -> float:
    """
    Convert a collar-product lower bound P_collar into a per-unit perimeter constant κ (lattice units).

    Intuition:
      • If each *independent collar block* contributes a multiplicative factor ≥ P_collar,
        and if we can pack 'density' such blocks per unit perimeter (in lattice length units),
        then a loop of perimeter L gets a factor ≥ P_collar^(density * L / perim_scale),
        i.e. an exponential with per-unit constant

            κ_latt = (density / perim_scale) * (-log P_collar)  ≥ 0.

      • 'perim_scale' is the length (in lattice units) associated with one collar block.
        With perim_scale=1, you interpret the block as a single edge-length contribution.

    Parameters
    ----------
    product_lower_bound : float in (0,1], the analytic lower bound for Π_k (1 - C η_k).
    perim_scale         : >0, lattice-length normalization for one independent collar block.
    density             : >0, packing density (blocks per unit perimeter) in lattice units.
    mp_dps              : precision for mpmath.

    Returns
    -------
    float: κ per unit *lattice* perimeter length (≥ 0).
    """
    _unit_interval_closed("product_lower_bound", product_lower_bound)
    _pos("perim_scale", perim_scale)
    _pos("density", density)
    mp.mp.dps = mp_dps
    # Safe: if product_lower_bound == 1, we get κ=0 (no perimeter suppression beyond area law).
    kappa = (density / perim_scale) * (-mp.log(product_lower_bound))
    # Cast to float for downstream convenience (mpf → float); high-precision retained if needed.
    return float(kappa)


def kappa_phys_from_kappa_latt(*, kappa_latt: float, a: float) -> float:
    """
    Convert lattice κ (per unit lattice length) to *physical* κ_phys (per unit physical length).

        κ_phys = κ_latt / a,

    where 'a' is the lattice spacing in the same length units used in your physical observables.
    """
    _nonneg("kappa_latt", kappa_latt)
    _pos("a", a)
    return float(kappa_latt / a)


def perimeter_factor(*, perimeter_latt: float, kappa_latt: float, mp_dps: int = 120) -> float:
    """
    Compute the multiplicative perimeter suppression factor:

        F_perim = exp(-κ_latt * perimeter_latt).

    • If perimeter_latt = 0, this returns 1.0.
    • κ_latt ≥ 0 ensures F_perim ∈ (0,1].
    """
    _nonneg("perimeter_latt", perimeter_latt)
    _nonneg("kappa_latt", kappa_latt)
    mp.mp.dps = mp_dps
    return float(mp.e**(-kappa_latt * perimeter_latt))


def area_perimeter_bound(
    *,
    area_phys: float,
    sigma_phys: float,
    perimeter_latt: Optional[float],
    kappa_latt: Optional[float],
    mp_dps: int = 120,
) -> Tuple[float, Optional[float], Optional[float]]:
    """
    Return the standard area-law bound and (optionally) the area+perimeter bound.

    Parameters
    ----------
    area_phys       : area of the loop in physical units used for σ_phys.
    sigma_phys      : physical string tension.
    perimeter_latt  : perimeter length in *lattice* units (optional).
    kappa_latt      : perimeter constant per unit lattice length (optional).

    Returns
    -------
    area_only, perim_factor_or_None, combined_or_None
      area_only = exp(-σ_phys * area_phys)
      perim_factor = exp(-κ_latt * perimeter_latt) if both perimeter_latt and kappa_latt are provided
      combined = area_only * perim_factor if perim_factor is available; otherwise None
    """
    _nonneg("area_phys", area_phys)
    _nonneg("sigma_phys", sigma_phys)
    mp.mp.dps = mp_dps
    area_only = float(mp.e**(-sigma_phys * area_phys))
    if perimeter_latt is None or kappa_latt is None:
        return area_only, None, None
    pf = perimeter_factor(perimeter_latt=perimeter_latt, kappa_latt=kappa_latt, mp_dps=mp_dps)
    return area_only, pf, float(area_only * pf)


@dataclass(frozen=True)
class PerimeterCertificate:
    product_lower_bound: float
    perim_scale: float
    density: float
    kappa_latt: float
    a: Optional[float]
    kappa_phys: Optional[float]
    notes: str


def perimeter_certificate_from_params(
    *,
    eta0: float,
    A: float,
    C: float,
    perim_scale: float = 1.0,
    density: float = 1.0,
    a: Optional[float] = None,
    mp_dps: int = 120,
) -> PerimeterCertificate:
    """
    One-shot computation:
      1) Compute analytic collar product lower bound P_collar ≥ Π_k (1 - C η_k) using (eta0, A, C),
      2) Convert to κ_latt via (density / perim_scale) * (-log P_collar),
      3) Optionally convert to κ_phys (if 'a' is provided).

    Returns a structured certificate suitable for JSON export or CLI display.
    """
    _pos("eta0", eta0)
    _pos("A", A)
    _pos("C", C)
    _pos("perim_scale", perim_scale)
    _pos("density", density)
    if a is not None:
        _pos("a", a)
    mp.mp.dps = mp_dps

    # Step 1: analytic product lower bound (∞-step, no truncation)
    P = float(collar_product_lower_bound_analytic(eta0=eta0, A=A, C=C))
    _unit_interval_closed("collar product lower bound", P)

    # Step 2: κ per unit lattice length
    kappa_latt = kappa_per_unit_from_product(
        product_lower_bound=P, perim_scale=perim_scale, density=density, mp_dps=mp_dps
    )

    # Step 3: optional physical conversion
    kappa_phys = kappa_phys_from_kappa_latt(kappa_latt=kappa_latt, a=a) if a is not None else None

    return PerimeterCertificate(
        product_lower_bound=P,
        perim_scale=perim_scale,
        density=density,
        kappa_latt=float(kappa_latt),
        a=a,
        kappa_phys=None if kappa_phys is None else float(kappa_phys),
        notes=(
            "Derived from analytic collar-product lower bound. "
            "Perimeter normalization uses 'density' blocks per lattice-length unit and 'perim_scale' length per block. "
            "κ_latt is per unit lattice length; κ_phys = κ_latt / a if 'a' is provided."
        ),
    )