# src/ym_bounds/cli.py
from __future__ import annotations

import argparse
import csv
import json
import sys
import subprocess
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any, Optional

from ym_bounds.contraction import (
    contraction_sequence,
    check_summability_and_product,
    infinite_sum_upper_bound,
    infinite_square_sum_upper_bound,
)
from ym_bounds.collar import (
    collar_product_bound,
    collar_product_lower_bound_analytic,
)
from ym_bounds.su3 import tail_bound_su3_sharp
from ym_bounds.tube_cost import spectral_gap_from_tube_cost
from ym_bounds.string_tension import sigma_phys_from_lattice, area_law_upper_bound
from ym_bounds.clustering import exp_clustering_bound
from ym_bounds.rg import (
    rp_block_certificate,
    eta0_bare_upper_bound,
)
from ym_bounds.perimeter import (
    kappa_per_unit_from_product,
    kappa_phys_from_kappa_latt,
    area_perimeter_bound,
)

TOOL_NAME = "ym-bounds"
TOOL_VERSION = "0.5.0"  # perimeter integration


# ----------------------------- Data models -----------------------------

@dataclass(frozen=True)
class ReportInputs:
    beta: float
    p0: Optional[int]
    q0: Optional[int]
    sym_N: Optional[int]
    eta0: float
    A: float
    C: float
    steps: int
    tau0: float
    sigma_lat: float
    a: float
    area: float
    mstar: float
    pref: float
    mp_dps: int
    perim_scale: float
    perim_density: float
    perimeter: Optional[float]

@dataclass(frozen=True)
class ReportOutputs:
    # contraction / collar (numeric)
    sum_eta_numeric: float
    collar_prod_numeric: float
    sum_eta_finite: bool
    collar_prod_positive: bool
    # contraction / collar (analytic infinite-steps)
    sum_eta_inf_bound: float
    sum_eta_sq_inf_bound: float
    collar_prod_analytic_lb: float
    # su3 tail
    low_shell_sum: float
    high_shell_tail: float
    # tube-cost => spectral gap
    lam_below_1: float
    m0_lb: float
    # string tension & area law
    sigma_phys: float
    wl_bound_area_only: float
    # clustering
    corr_bound_at_1: float
    # perimeter term
    kappa_latt: float
    kappa_phys: float
    wl_perim_factor: Optional[float]
    wl_bound_area_perim: Optional[float]

@dataclass(frozen=True)
class ReportMeta:
    timestamp_utc: str
    tool: str
    tool_version: str
    python_version: str
    mp_dps: int
    git_commit: str


# ----------------------------- Utilities -----------------------------

def _git_commit() -> str:
    try:
        out = subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL)
        return out.decode("utf-8").strip()
    except Exception:
        return "unknown"

def _validate_report_inputs(inp: ReportInputs) -> None:
    if inp.beta <= 0:
        raise ValueError("--beta must be > 0")
    if inp.sym_N is None:
        if inp.p0 is None or inp.q0 is None:
            raise ValueError("Either specify --sym-N, or specify both --p0 and --q0.")
        if inp.p0 < 0 or inp.q0 < 0:
            raise ValueError("--p0 and --q0 must be >= 0")
    else:
        if inp.sym_N < 1:
            raise ValueError("--sym-N must be >= 1")
    if min(inp.eta0, inp.A, inp.C, inp.tau0, inp.sigma_lat, inp.a, inp.mstar, inp.pref, inp.perim_scale, inp.perim_density) <= 0:
        raise ValueError("eta0, A, C, tau0, sigma-lat, a, mstar, pref, perim-scale, perim-density must be > 0")
    if inp.steps <= 0:
        raise ValueError("--steps must be a positive integer")
    if inp.area < 0:
        raise ValueError("--area must be >= 0")
    if inp.mp_dps < 50:
        raise ValueError("--mp-dps should be at least 50 for reliability")
    if inp.perimeter is not None and inp.perimeter < 0:
        raise ValueError("--perimeter must be >= 0 if provided")


# ----------------------------- Core computations -----------------------------

def _compute_report(inp: ReportInputs) -> ReportOutputs:
    # 1) Contraction numeric (finite steps) + numeric collar
    etas = contraction_sequence(inp.eta0, inp.A, steps=inp.steps)
    sum_eta_num, prod_num, sum_ok, prod_ok = check_summability_and_product(etas, inp.C)

    # 2) Analytic, infinite-step certificates (contraction/collar)
    sum_eta_inf = infinite_sum_upper_bound(inp.eta0, inp.A)
    sum_eta_sq_inf = infinite_square_sum_upper_bound(inp.eta0, inp.A)
    collar_lb_analytic = collar_product_lower_bound_analytic(inp.eta0, inp.A, inp.C)

    # 3) SU(3) tail with optional symmetric-cut semantics
    if inp.sym_N is not None:
        p0 = q0 = inp.sym_N
    else:
        p0, q0 = inp.p0, inp.q0  # type: ignore[assignment]
    partial, tail = tail_bound_su3_sharp(inp.beta, p0, q0, mp_dps=inp.mp_dps)

    # 4) Tube-cost ⇒ spectral gap
    lam_below_1, m0_lb = spectral_gap_from_tube_cost(inp.tau0)

    # 5) String tension + area law (area-only)
    sigma_phys = sigma_phys_from_lattice(inp.sigma_lat, inp.a)
    wl_area_only = area_law_upper_bound(inp.area, sigma_phys)

    # 6) Clustering at |x|=1
    corr_bound_at_1 = exp_clustering_bound(distance=1.0, m_star=inp.mstar, prefactor=inp.pref)

    # 7) Perimeter constant from analytic collar product
    kappa_latt = kappa_per_unit_from_product(
        product_lower_bound=float(collar_lb_analytic),
        perim_scale=inp.perim_scale,
        density=inp.perim_density,
        mp_dps=inp.mp_dps,
    )
    kappa_phys = kappa_phys_from_kappa_latt(kappa_latt=kappa_latt, a=inp.a)

    # 8) Optional perimeter factor and combined area+perimeter bound
    wl_perim_factor = None
    wl_area_perim = None
    if inp.perimeter is not None:
        area_only, perim_factor, combined = area_perimeter_bound(
            area_phys=inp.area, sigma_phys=sigma_phys,
            perimeter_latt=inp.perimeter, kappa_latt=kappa_latt,
            mp_dps=inp.mp_dps,
        )
        # area_only should equal wl_area_only; keep prints consistent with prior section
        wl_perim_factor = perim_factor
        wl_area_perim = combined

    return ReportOutputs(
        sum_eta_numeric=float(sum_eta_num),
        collar_prod_numeric=float(prod_num),
        sum_eta_finite=bool(sum_ok),
        collar_prod_positive=bool(prod_ok),
        sum_eta_inf_bound=float(sum_eta_inf),
        sum_eta_sq_inf_bound=float(sum_eta_sq_inf),
        collar_prod_analytic_lb=float(collar_lb_analytic),
        low_shell_sum=float(partial),
        high_shell_tail=float(tail),
        lam_below_1=float(lam_below_1),
        m0_lb=float(m0_lb),
        sigma_phys=float(sigma_phys),
        wl_bound_area_only=float(wl_area_only),
        corr_bound_at_1=float(corr_bound_at_1),
        kappa_latt=float(kappa_latt),
        kappa_phys=float(kappa_phys),
        wl_perim_factor=None if wl_perim_factor is None else float(wl_perim_factor),
        wl_bound_area_perim=None if wl_area_perim is None else float(wl_area_perim),
    )


def _print_report_console(inp: ReportInputs, out: ReportOutputs, meta: ReportMeta) -> None:
    # Identify the cut semantics for display
    if inp.sym_N is not None:
        cut_desc = f"(symmetric cut N={inp.sym_N})"
    else:
        cut_desc = f"(p0,q0)=({inp.p0},{inp.q0})"

    print(f"\n=== {meta.tool} Report (v{meta.tool_version}) ===")
    print(f"[Meta] time_utc={meta.timestamp_utc} python={meta.python_version} mp_dps={meta.mp_dps} git={meta.git_commit}")

    print(f"\n[Inputs]")
    print(f"  β={inp.beta:.6g}  {cut_desc}")
    print(f"  η0={inp.eta0:.6g}  A={inp.A:.6g}  C={inp.C:.6g}  steps={inp.steps}")
    print(f"  τ0={inp.tau0:.6g}  σ_lat={inp.sigma_lat:.6g}  a={inp.a:.6g}  area={inp.area:.6g}")
    print(f"  clustering: m*={inp.mstar:.6g}, pref={inp.pref:.6g}")
    print(f"  perimeter options: perim-scale={inp.perim_scale:.6g}, perim-density={inp.perim_density:.6g}, perimeter={inp.perimeter}")

    print(f"\n[Contraction — numeric vs analytic]")
    print(f"  Ση_k (numeric, {inp.steps} steps) = {out.sum_eta_numeric:.6g}  (finite? {out.sum_eta_finite})")
    print(f"  Ση_k (analytic, ∞ steps) ≤ {out.sum_eta_inf_bound:.6g}")
    print(f"  Ση_k^2 (analytic, ∞ steps) ≤ {out.sum_eta_sq_inf_bound:.6g}")

    print(f"\n[Collar product]")
    print(f"  Π(1 - C η_k) (numeric, {inp.steps} steps) = {out.collar_prod_numeric:.6g}  (positive? {out.collar_prod_positive})")
    print(f"  Π(1 - C η_k) (analytic LB, ∞ steps) ≥ {out.collar_prod_analytic_lb:.6g}")

    print(f"\n[SU(3) Tail] (sharp shell-min bound, mp_dps={meta.mp_dps})")
    print(f"  Low-shell sum ≤ {out.low_shell_sum}")
    print(f"  High-shell tail ≤ {out.high_shell_tail}")

    print(f"\n[Tube-Cost ⇒ Gap]")
    print(f"  Spec(T) ⊂ {{1}} ∪ [{out.lam_below_1:.6g}, 1)  ⇒  m0 ≥ {out.m0_lb:.6g}")

    print(f"\n[String Tension]")
    print(f"  σ_phys = {out.sigma_phys:.6g}")
    print(f"  Area-law bound for Area={inp.area:.6g}:  ⟨W⟩ ≤ {out.wl_bound_area_only:.6g}")

    print(f"\n[Perimeter]")
    print(f"  κ_latt (per unit lattice length) = {out.kappa_latt:.6g}")
    print(f"  κ_phys (per unit physical length) = {out.kappa_phys:.6g}")
    if inp.perimeter is not None:
        print(f"  Perimeter factor for Perim={inp.perimeter:.6g}:  ≤ {out.wl_perim_factor:.6g}")
        if out.wl_bound_area_perim is not None:
            print(f"  Combined area+perimeter bound:  ⟨W⟩ ≤ {out.wl_bound_area_perim:.6g}")

    print(f"\n[Clustering]")
    print(f"  |⟨F(x)F(0)⟩| ≤ {inp.pref:.6g}·e^(-{inp.mstar:.6g}·|x|)  ⇒  at |x|=1: ≤ {out.corr_bound_at_1:.6g}")

    print("\nAll quantities computed with rigorous inequalities or saturated worst-cases;")
    print("analytic certificates do not depend on finite-step truncation.")


def _write_report_csv(path: Path, inp: ReportInputs, out: ReportOutputs, meta: ReportMeta) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row: Dict[str, Any] = {
        **asdict(inp),
        **asdict(out),
        "timestamp_utc": meta.timestamp_utc,
        "tool": meta.tool,
        "tool_version": meta.tool_version,
        "python_version": meta.python_version,
        "git_commit": meta.git_commit,
    }
    cols = [
        # meta
        "timestamp_utc","tool","tool_version","python_version","git_commit",
        # inputs
        "beta","p0","q0","sym_N","eta0","A","C","steps","tau0","sigma_lat","a","area","mstar","pref","mp_dps",
        "perim_scale","perim_density","perimeter",
        # outputs
        "sum_eta_numeric","sum_eta_inf_bound","sum_eta_sq_inf_bound",
        "collar_prod_numeric","collar_prod_analytic_lb","sum_eta_finite","collar_prod_positive",
        "low_shell_sum","high_shell_tail","lam_below_1","m0_lb",
        "sigma_phys","wl_bound_area_only",
        "kappa_latt","kappa_phys","wl_perim_factor","wl_bound_area_perim",
        "corr_bound_at_1",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerow({k: row.get(k, "") for k in cols})


def _write_report_tex(path: Path, inp: ReportInputs, out: ReportOutputs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cut_line = f"(symmetric N={inp.sym_N})" if inp.sym_N is not None else f"(p0={inp.p0}, q0={inp.q0})"
    tex = rf"""\begin{{tabular}}{{l r}}
\hline
Cut {cut_line} & -- \\
$\sum_k \eta_k$ (numeric) & {out.sum_eta_numeric:.6g} \\
$\sum_k \eta_k$ (analytic, $\infty$) & \le {out.sum_eta_inf_bound:.6g} \\
$\sum_k \eta_k^2$ (analytic, $\infty$) & \le {out.sum_eta_sq_inf_bound:.6g} \\
$\prod_k(1 - C\,\eta_k)$ (numeric) & {out.collar_prod_numeric:.6g} \\
$\prod_k(1 - C\,\eta_k)$ (analytic LB) & \ge {out.collar_prod_analytic_lb:.6g} \\
Low-shell sum & {out.low_shell_sum:.6g} \\
High-shell tail & {out.high_shell_tail:.6g} \\
$m_0$ (from $\tau_0$) & \ge {out.m0_lb:.6g} \\
$\sigma_{{\rm phys}}$ & {out.sigma_phys:.6g} \\
Area-law bound (Area={inp.area:.6g}) & {out.wl_bound_area_only:.6g} \\
$\kappa_{{\rm latt}}$ (per unit) & {out.kappa_latt:.6g} \\
$\kappa_{{\rm phys}}$ (per unit) & {out.kappa_phys:.6g} \\
\hline
\end{{tabular}}
"""
    path.write_text(tex)


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


# ----------------------------- Subcommands -----------------------------

def _cmd_report(ns: argparse.Namespace) -> int:
    # Load constants from JSON if provided
    if ns.constants:
        print(f"📥 Loading constants from: {ns.constants}")
        try:
            with open(ns.constants, 'r') as f:
                constants = json.load(f)
            
            # Override CLI arguments with loaded constants
            if ns.eta0 is None:
                ns.eta0 = constants.get("eta0_estimate", constants.get("eta0"))
            if ns.A is None:
                ns.A = constants.get("A")
            if ns.C is None:
                ns.C = constants.get("C")
            if ns.tau0 is None:
                ns.tau0 = constants.get("tau0")
                
            print(f"   ✅ Loaded: A={ns.A}, C={ns.C}, tau0={ns.tau0}, eta0={ns.eta0}")
            if "metadata" in constants:
                meta = constants["metadata"]
                print(f"   📋 Source: {meta.get('tool', 'unknown')} v{meta.get('version', 'unknown')}")
                
        except Exception as e:
            print(f"   ❌ Error loading constants: {e}")
            return 1
    
    # Validate required arguments
    missing = []
    if ns.eta0 is None: missing.append("--eta0")
    if ns.A is None: missing.append("--A") 
    if ns.C is None: missing.append("--C")
    if ns.tau0 is None: missing.append("--tau0")
    if missing:
        print(f"❌ Missing required arguments: {', '.join(missing)}")
        print("   Either provide them explicitly or use --constants <file.json>")
        return 1
    
    p0 = ns.p0 if ns.p0 is not None else None
    q0 = ns.q0 if ns.q0 is not None else None
    sym_N = ns.sym_N if ns.sym_N is not None else None

    inp = ReportInputs(
        beta=ns.beta, p0=p0, q0=q0, sym_N=sym_N,
        eta0=ns.eta0, A=ns.A, C=ns.C, steps=ns.steps,
        tau0=ns.tau0, sigma_lat=ns.sigma_lat, a=ns.a,
        area=ns.area, mstar=ns.mstar, pref=ns.pref, mp_dps=ns.mp_dps,
        perim_scale=ns.perim_scale, perim_density=ns.perim_density, perimeter=ns.perimeter,
    )
    _validate_report_inputs(inp)
    out = _compute_report(inp)
    meta = ReportMeta(
        timestamp_utc=datetime.now(tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        tool=TOOL_NAME,
        tool_version=TOOL_VERSION,
        python_version=f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        mp_dps=ns.mp_dps,
        git_commit=_git_commit(),
    )
    _print_report_console(inp, out, meta)
    if ns.csv:
        _write_report_csv(Path(ns.csv), inp, out, meta)
        print(f"[CSV] wrote {ns.csv}")
    if ns.tex:
        _write_report_tex(Path(ns.tex), inp, out)
        print(f"[TeX] wrote {ns.tex}")
    if ns.json:
        _write_json(Path(ns.json), {"meta": asdict(meta), "inputs": asdict(inp), "outputs": asdict(out)})
        print(f"[JSON] wrote {ns.json}")
    return 0


def _cmd_smallness_lemma(ns: argparse.Namespace) -> int:
    eta0_bound = eta0_bare_upper_bound(beta=ns.beta, blocks=ns.blocks, sym_N=ns.sym_N, mp_dps=ns.mp_dps)
    ok = (eta0_bound <= ns.threshold)
    meta = {
        "timestamp_utc": datetime.now(tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "tool": TOOL_NAME,
        "tool_version": TOOL_VERSION,
        "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "git_commit": _git_commit(),
        "mp_dps": ns.mp_dps,
    }
    print("\n=== Initial Smallness Lemma ===")
    print(f"[Meta] time_utc={meta['timestamp_utc']} python={meta['python_version']} mp_dps={meta['mp_dps']} git={meta['git_commit']}")
    print(f"[Inputs] β={ns.beta:.6g}  blocks={ns.blocks}  symmetric cut N={ns.sym_N}  threshold={ns.threshold:.6g}")
    print(f"[Result] η0 ≤ {float(eta0_bound):.12g}   ⇒  OK vs threshold? {ok}")

    payload = {
        "meta": meta,
        "inputs": {"beta": ns.beta, "blocks": ns.blocks, "sym_N": ns.sym_N, "threshold": ns.threshold, "mp_dps": ns.mp_dps},
        "outputs": {"eta0_upper_bound": float(eta0_bound), "ok": bool(ok)},
        "notes": "Bound derived from nontrivial SU(3) heat-kernel coefficients after RP-preserving blocks; symmetric shell cut used for tail bound.",
    }

    if ns.json:
        _write_json(Path(ns.json), payload)
        print(f"[JSON] wrote {ns.json}")
    if ns.csv:
        cols = ["beta","blocks","sym_N","threshold","eta0_upper_bound","ok","mp_dps","git_commit","timestamp_utc","tool_version"]
        row = {
            "beta": ns.beta, "blocks": ns.blocks, "sym_N": ns.sym_N, "threshold": ns.threshold,
            "eta0_upper_bound": float(eta0_bound), "ok": ok,
            "mp_dps": ns.mp_dps, "git_commit": meta["git_commit"], "timestamp_utc": meta["timestamp_utc"], "tool_version": TOOL_VERSION,
        }
        path = Path(ns.csv); path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerow(row)
        print(f"[CSV] wrote {ns.csv}")
    return 0


def _cmd_rp_cert(ns: argparse.Namespace) -> int:
    cert = rp_block_certificate(beta=ns.beta, blocks=ns.blocks, sym_N=ns.sym_N, mp_dps=ns.mp_dps)
    print("\n=== RP-Preserving RG Certificate ===")
    print(f"β_in={cert['beta_input']:.6g}  blocks={cert['blocks']}  β_blocked={cert['beta_blocked']:.6g}  sym_N={cert['sym_N']}  mp_dps={cert['mp_dps']}")
    print(f"coeff_nonneg={cert['coeff_nonneg']}  class_function={cert['class_function']}  gauge_invariant={cert['gauge_invariant']}  rp_preserved={cert['rp_preserved']}")
    print(f"locality_radius={cert['locality_radius']}  nontrivial_weight={cert['nontrivial_weight']:.12g}")
    print(f"notes: {cert['notes']}")
    if ns.json:
        payload = {
            "meta": {
                "timestamp_utc": datetime.now(tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "tool": TOOL_NAME, "tool_version": TOOL_VERSION,
                "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
                "git_commit": _git_commit(),
            },
            "certificate": cert,
        }
        _write_json(Path(ns.json), payload)
        print(f"[JSON] wrote {ns.json}")
    return 0


# ----------------------------- Parser wiring -----------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=TOOL_NAME,
        description="Compute explicit bounds and export machine-verifiable constants for the SU(3) mass-gap pipeline.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # report
    rp = sub.add_parser("report", help="Compute a full bounds report.")
    rp.add_argument("--beta", type=float, required=True)
    # Either symmetric cut OR explicit p0,q0
    rp.add_argument("--sym-N", type=int, help="Use symmetric cut p0=q0=N.")
    rp.add_argument("--p0", type=int, help="Lower bound p0 (ignored if --sym-N is set).")
    rp.add_argument("--q0", type=int, help="Lower bound q0 (ignored if --sym-N is set).")
    rp.add_argument("--eta0", type=float, help="Initial KP norm (required unless --constants provided)")
    rp.add_argument("--A", type=float, help="Contraction constant A (required unless --constants provided)")
    rp.add_argument("--C", type=float, help="Collar constant C (required unless --constants provided)")
    rp.add_argument("--steps", type=int, default=20)
    rp.add_argument("--tau0", type=float, help="Tube cost parameter (required unless --constants provided)")
    rp.add_argument("--sigma-lat", type=float, required=True)
    rp.add_argument("--a", type=float, required=True)
    rp.add_argument("--area", type=float, default=1.0)
    rp.add_argument("--mstar", type=float, default=0.3)
    rp.add_argument("--pref", type=float, default=1.0)
    # NEW: constants import
    rp.add_argument("--constants", type=Path, help="Import constants from JSON file (from RG step analysis)")
    # NEW: perimeter options
    rp.add_argument("--perim-scale", type=float, default=1.0, help="Lattice-length normalization per independent collar block.")
    rp.add_argument("--perim-density", type=float, default=1.0, help="Independent collar blocks per unit lattice perimeter.")
    rp.add_argument("--perimeter", type=float, help="Loop perimeter in lattice units (optional; enables combined area+perimeter bound).")
    rp.add_argument("--mp-dps", type=int, default=80, help="mpmath precision (decimal digits).")
    rp.add_argument("--csv", type=Path)
    rp.add_argument("--tex", type=Path)
    rp.add_argument("--json", type=Path)

    # smallness-lemma
    sl = sub.add_parser("smallness-lemma", help="Derive an analytic upper bound on η0 from the Wilson/heat-kernel action after RP-preserving blocks.")
    sl.add_argument("--beta", type=float, required=True, help="Base β in exp(-β C2/6).")
    sl.add_argument("--blocks", type=int, default=1, help="Number of RP-preserving convolution blocks.")
    sl.add_argument("--sym-N", type=int, default=8, help="Symmetric shell cut N (p0=q0=N) for tail control.")
    sl.add_argument("--threshold", type=float, default=0.05, help="Pass/fail threshold for η0 (default 0.05).")
    sl.add_argument("--mp-dps", type=int, default=120, help="mpmath precision (decimal digits).")
    sl.add_argument("--json", type=Path)
    sl.add_argument("--csv", type=Path)

    # rp-cert
    rc = sub.add_parser("rp-cert", help="Emit a machine-checkable RP-preserving RG certificate (convolution of SU(3) heat kernels).")
    rc.add_argument("--beta", type=float, required=True)
    rc.add_argument("--blocks", type=int, default=1)
    rc.add_argument("--sym-N", type=int, default=8)
    rc.add_argument("--mp-dps", type=int, default=120)
    rc.add_argument("--json", type=Path)

    return p


def main() -> None:
    parser = _build_parser()
    try:
        ns = parser.parse_args()
        if ns.cmd == "report":
            sys.exit(_cmd_report(ns))
        elif ns.cmd == "smallness-lemma":
            sys.exit(_cmd_smallness_lemma(ns))
        elif ns.cmd == "rp-cert":
            sys.exit(_cmd_rp_cert(ns))
        else:
            parser.error("Unknown command")
    except KeyboardInterrupt:
        print("\nAborted by user.", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()