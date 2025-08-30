# src/ym_research/BoundDiscovererToolCLI.py
from __future__ import annotations

"""
Research CLI for exploring and stress-testing the SU(3) mass-gap bounds.

Commands:
  - sweep-contraction : grid sweep over (eta0, A) and check summability/product positivity
  - optimize-tail     : choose (p0,q0) or symmetric N so SU(3) tail ≤ tolerance
  - target-gap        : compute (λ_below_1, m0) from a chosen τ0
  - string-table      : tabulate σ_phys and area-law bounds for a list/range of spacings a
  - a-sweep           : demonstrate uniformity in a, and print η0 bound & collar LB independent of a
  - recipe            : show curated examples

All outputs are deterministic and exportable (CSV/JSON/TeX where appropriate).
"""

import argparse
import csv
import json
import sys
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any, Tuple, List, Optional

# --- Imports from your ym_bounds package ---
from ym_bounds.contraction import contraction_sequence, check_summability_and_product
from ym_bounds.su3 import tail_bound_su3_sharp
from ym_bounds.tube_cost import spectral_gap_from_tube_cost
from ym_bounds.string_tension import sigma_phys_from_lattice, area_law_upper_bound
from ym_bounds.clustering import exp_clustering_bound
from ym_bounds.rg import eta0_bare_upper_bound
from ym_bounds.collar import collar_product_lower_bound_analytic

TOOL = "ym-research"
VERSION = "0.3.0"  # includes a-sweep, optimize-tail precision control, and robust exports


# ------------------------- Utilities & meta -------------------------

@dataclass(frozen=True)
class Meta:
    timestamp_utc: str
    tool: str
    version: str
    python: str
    mp_dps: int

def _meta(mp_dps: int) -> Meta:
    return Meta(
        timestamp_utc=datetime.now(tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        tool=TOOL,
        version=VERSION,
        python=f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        mp_dps=mp_dps,
    )

def _write_csv(path: Path, rows: List[Dict[str, Any]], header: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in header})

def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))

def _pos(name: str, x: float) -> None:
    if not (isinstance(x, (int, float)) and x > 0):
        raise ValueError(f"{name} must be > 0")

def _nonneg(name: str, x: float) -> None:
    if not (isinstance(x, (int, float)) and x >= 0):
        raise ValueError(f"{name} must be >= 0")

def _parse_a_list(ns: argparse.Namespace) -> List[float]:
    vals: List[float] = []
    if ns.a_list:
        parts = [t.strip() for t in ns.a_list.split(",") if t.strip()]
        if not parts:
            raise ValueError("--a-list provided but empty")
        for t in parts:
            try:
                v = float(t)
            except Exception:
                raise ValueError(f"--a-list contains non-numeric value: {t}")
            _pos("a", v)
            vals.append(v)
    else:
        # range mode
        if ns.a_min is None or ns.a_max is None or ns.a_num is None:
            raise ValueError("Provide either --a-list OR all of --a-min, --a-max, --a-num")
        _pos("--a-min", ns.a_min)
        _pos("--a-max", ns.a_max)
        if not (isinstance(ns.a_num, int) and ns.a_num >= 2):
            raise ValueError("--a-num must be an integer >= 2")
        if ns.a_max <= ns.a_min:
            raise ValueError("--a-max must be > --a-min")
        step = (ns.a_max - ns.a_min) / (ns.a_num - 1)
        for i in range(ns.a_num):
            vals.append(ns.a_min + i * step)
    return vals


# ------------------------- Commands -------------------------

def cmd_sweep_contraction(ns: argparse.Namespace) -> int:
    _pos("--C", ns.C)
    if ns.steps <= 0:
        raise ValueError("--steps must be positive")
    if ns.grid_eta <= 0 or ns.grid_A <= 0:
        raise ValueError("--grid-eta and --grid-A must be positive integers")

    eta_vals = [ns.eta0_min + i*(ns.eta0_max - ns.eta0_min)/(ns.grid_eta - 1)
                for i in range(ns.grid_eta)]
    A_vals = [ns.A_min + j*(ns.A_max - ns.A_min)/(ns.grid_A - 1)
              for j in range(ns.grid_A)]
    rows: List[Dict[str, Any]] = []
    ok_count = 0
    for e in eta_vals:
        _pos("eta0", e)
        for A in A_vals:
            _pos("A", A)
            etas = contraction_sequence(e, A, steps=ns.steps)
            s, prod, s_ok, p_ok = check_summability_and_product(etas, ns.C)
            rows.append({
                "eta0": e, "A": A, "steps": ns.steps, "C": ns.C,
                "sum_eta": float(s),
                "product": float(prod),
                "sum_finite": bool(s_ok),
                "product_positive": bool(p_ok),
            })
            if s_ok and p_ok:
                ok_count += 1

    meta = _meta(mp_dps=80)
    print(f"\n[{TOOL}] sweep-contraction: {ok_count}/{len(rows)} grid points OK")
    print(f"  grid_eta={ns.grid_eta} grid_A={ns.grid_A} steps={ns.steps} C={ns.C}")
    if ns.csv:
        header = ["eta0", "A", "steps", "C", "sum_eta", "product", "sum_finite", "product_positive"]
        _write_csv(Path(ns.csv), rows, header)
        print(f"[CSV] wrote {ns.csv}")
    if ns.json:
        _write_json(Path(ns.json), {"meta": asdict(meta), "rows": rows})
        print(f"[JSON] wrote {ns.json}")
    return 0


def _tail_ok(beta: float, p0: int, q0: int, tol: float, mp_dps: int) -> Tuple[float, float, bool]:
    partial, tail = tail_bound_su3_sharp(beta, p0, q0, mp_dps=mp_dps)
    return float(partial), float(tail), (tail <= tol)

def cmd_optimize_tail(ns: argparse.Namespace) -> int:
    _pos("--beta", ns.beta)
    _pos("--tol", ns.tol)
    if ns.strategy not in {"symmetric", "strip"}:
        raise ValueError("--strategy must be 'symmetric' or 'strip'")
    if ns.N_max < ns.N_min:
        raise ValueError("N-max must be >= N-min")
    mp_dps = ns.mp_dps if ns.mp_dps else 80

    N = max(0, ns.N_min)
    best: Dict[str, Any] | None = None
    while N <= ns.N_max:
        p0, q0 = (N, N) if ns.strategy == "symmetric" else (N, 0)
        partial, tail, ok = _tail_ok(ns.beta, p0, q0, ns.tol, mp_dps)
        if best is None or (ok and N < best["N"]):
            best = {"N": N, "p0": p0, "q0": q0, "partial": partial, "tail": tail, "ok": ok}
        if ok:
            break
        N += 1

    meta = _meta(mp_dps=mp_dps)
    print(f"\n[{TOOL}] optimize-tail (β={ns.beta}, tol={ns.tol}, strategy={ns.strategy}, mp_dps={mp_dps})")
    if best:
        print(f"  result: N={best['N']}  (p0={best['p0']}, q0={best['q0']})  tail≤tol? {best['ok']}")
        print(f"  partial≤ {best['partial']}   tail≤ {best['tail']}")
    else:
        print("  no candidate satisfied tolerance in given range")
    if ns.json:
        _write_json(Path(ns.json), {"meta": asdict(meta), "result": best})
        print(f"[JSON] wrote {ns.json}")
    return 0


def cmd_target_gap(ns: argparse.Namespace) -> int:
    _pos("--m0", ns.m0)
    lam, m0_lb = spectral_gap_from_tube_cost(ns.m0)
    print(f"\n[{TOOL}] target-gap: m0_target={ns.m0}")
    print(f"  choose τ0={ns.m0} ⇒ Spec(T) ⊂ {{1}} ∪ [{lam:.6g}, 1)")
    print(f"  certified m0 ≥ {m0_lb:.6g}")
    return 0


def cmd_string_table(ns: argparse.Namespace) -> int:
    _pos("--sigma-lat", ns.sigma_lat)
    a_vals: List[float] = []
    if ns.a_list:
        parts = [t.strip() for t in ns.a_list.split(",") if t.strip()]
        for t in parts:
            v = float(t)
            _pos("a", v)
            a_vals.append(v)
    else:
        _pos("--a", ns.a)
        a_vals = [ns.a]

    rows: List[Dict[str, Any]] = []
    for a in a_vals:
        sigma_phys = sigma_phys_from_lattice(ns.sigma_lat, a)
        wl = area_law_upper_bound(ns.area, sigma_phys)
        rows.append({"sigma_lat": ns.sigma_lat, "a": a, "sigma_phys": sigma_phys, "area": ns.area, "wl_bound": wl})

    meta = _meta(mp_dps=80)
    print(f"\n[{TOOL}] string-tension table  (σ_lat={ns.sigma_lat}, area={ns.area})")
    for r in rows:
        print(f"  a={r['a']:.6g}  σ_phys={r['sigma_phys']:.6g}  ⟨W⟩≤{r['wl_bound']:.6g}")

    if ns.csv:
        _write_csv(Path(ns.csv), rows, ["sigma_lat", "a", "sigma_phys", "area", "wl_bound"])
        print(f"[CSV] wrote {ns.csv}")
    if ns.json:
        _write_json(Path(ns.json), {"meta": asdict(meta), "rows": rows})
        print(f"[JSON] wrote {ns.json}")
    return 0


def cmd_recipe(_: argparse.Namespace) -> int:
    print(f"""
[{TOOL}] recipes — exploration commands

1) Map safe contraction region:
   ym-research sweep-contraction --eta0-min 0.02 --eta0-max 0.08 --grid-eta 7 \\
     --A-min 2.0 --A-max 5.0 --grid-A 7 --steps 20 --C 0.2 --csv sweep.csv

2) Choose (p0,q0) so SU(3) tail ≤ 1e-6 at β=6 (symmetric cut, higher precision):
   ym-research optimize-tail --beta 6.0 --tol 1e-6 --strategy symmetric --N-min 4 --N-max 30 --mp-dps 120

3) Target spectral gap:
   ym-research target-gap --m0 0.5

4) String-tension table across 'a':
   ym-research string-table --sigma-lat 0.045 --a-list 0.12,0.10,0.08 --area 1 --csv sigma.csv

5) Uniform-in-a sweep with Initial Smallness and collar LB:
   ym-research a-sweep --beta 6.0 --blocks 2 --sym-N 8 --A 3.0 --C 0.2 --tau0 0.4 \\
     --sigma-lat 0.045 --a-list 0.12,0.10,0.08 --area 1 --mp-dps 120 --csv a_sweep.csv
""".rstrip())
    return 0


def cmd_a_sweep(ns: argparse.Namespace) -> int:
    """
    Sweep lattice spacings 'a' and report, for each:
      - sigma_phys = sigma_lat / a^2
      - Wilson area-law upper bound for the given area
      - Initial Smallness Lemma bound: eta0 ≤ S_nontriv(beta_b) with beta_b = blocks*beta, symmetric cut N
      - Analytic collar product lower bound using eta0_bound (step-free)
      - Tube-cost ⇒ spectral gap: (lambda_below_1, m0_lb)
    """
    _pos("--beta", ns.beta)
    if ns.blocks <= 0:
        raise ValueError("--blocks must be a positive integer")
    if ns.sym_N < 1:
        raise ValueError("--sym-N must be >= 1")
    _pos("--A", ns.A)
    _pos("--C", ns.C)
    _pos("--tau0", ns.tau0)
    _pos("--sigma-lat", ns.sigma_lat)
    _nonneg("--area", ns.area)
    if ns.mp_dps < 60:
        raise ValueError("--mp-dps should be at least 60 for reliable tail/eta0 bounds")

    a_vals = _parse_a_list(ns)

    # Constants independent of 'a'
    meta = _meta(mp_dps=ns.mp_dps)
    eta0_bound = eta0_bare_upper_bound(beta=ns.beta, blocks=ns.blocks, sym_N=ns.sym_N, mp_dps=ns.mp_dps)
    collar_lb = collar_product_lower_bound_analytic(eta0=float(eta0_bound), A=ns.A, C=ns.C)
    lam, m0_lb = spectral_gap_from_tube_cost(ns.tau0)

    rows: List[Dict[str, Any]] = []
    for a in a_vals:
        _pos("a", a)
        sigma_phys = sigma_phys_from_lattice(ns.sigma_lat, a)
        wl = area_law_upper_bound(ns.area, sigma_phys)
        rows.append({
            "a": a,
            "sigma_lat": ns.sigma_lat,
            "sigma_phys": float(sigma_phys),
            "area": ns.area,
            "wl_bound": float(wl),
            "eta0_upper_bound": float(eta0_bound),
            "collar_prod_analytic_lb": float(collar_lb),
            "lam_below_1": float(lam),
            "m0_lb": float(m0_lb),
            "beta": ns.beta,
            "blocks": ns.blocks,
            "sym_N": ns.sym_N,
            "A": ns.A,
            "C": ns.C,
            "tau0": ns.tau0,
            "mp_dps": ns.mp_dps,
        })

    print(f"\n[{TOOL}] a-sweep — uniformity in a (β={ns.beta}, blocks={ns.blocks}, sym_N={ns.sym_N}, A={ns.A}, C={ns.C}, τ0={ns.tau0})")
    print(f"  (η0_upper_bound and collar LB are 'a'-independent; m0_lb from τ0 is 'a'-independent)")
    for r in rows:
        print(f"  a={r['a']:.6g}  σ_phys={r['sigma_phys']:.6g}  ⟨W⟩≤{r['wl_bound']:.6g}  "
              f"η0≤{r['eta0_upper_bound']:.6g}  Π(1−Cη)≥{r['collar_prod_analytic_lb']:.6g}  "
              f"Spec(T)⊂{{1}}∪[{r['lam_below_1']:.6g},1) ⇒ m0≥{r['m0_lb']:.6g}")

    if ns.csv:
        header = ["a","sigma_lat","sigma_phys","area","wl_bound","eta0_upper_bound",
                  "collar_prod_analytic_lb","lam_below_1","m0_lb",
                  "beta","blocks","sym_N","A","C","tau0","mp_dps"]
        _write_csv(Path(ns.csv), rows, header)
        print(f"[CSV] wrote {ns.csv}")
    if ns.json:
        payload = {"meta": asdict(meta), "rows": rows,
                   "notes": "Uniform-in-a sweep: η0 bound and collar LB are independent of a; σ_phys scales as 1/a^2; gap from τ0 is constant."}
        _write_json(Path(ns.json), payload)
        print(f"[JSON] wrote {ns.json}")
    if ns.tex:
        path = Path(ns.tex); path.parent.mkdir(parents=True, exist_ok=True)
        lines = ["\\begin{tabular}{r r r r r}",
                 "\\hline",
                 " $a$ & $\\sigma_{\\rm phys}$ & $\\langle W\\rangle$ bound & $\\eta_0$ upper & $\\prod(1-C\\eta)$ LB \\\\",
                 "\\hline"]
        for r in rows:
            lines.append(f" {r['a']:.6g} & {r['sigma_phys']:.6g} & {r['wl_bound']:.6g} & {r['eta0_upper_bound']:.6g} & {r['collar_prod_analytic_lb']:.6g} \\\\")
        lines.append("\\hline")
        lines.append("\\end{tabular}")
        path.write_text("\n".join(lines))
        print(f"[TeX] wrote {ns.tex}")
    return 0


# ------------------------- Parser wiring -------------------------

def _build() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=TOOL,
        description="Research CLI to explore and optimize explicit bounds for the SU(3) mass-gap pipeline.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("sweep-contraction", help="Grid-sweep (eta0,A) and check summability/product positivity.")
    sp.add_argument("--eta0-min", type=float, required=True)
    sp.add_argument("--eta0-max", type=float, required=True)
    sp.add_argument("--grid-eta", type=int, required=True, help="grid points along eta0")
    sp.add_argument("--A-min", type=float, required=True)
    sp.add_argument("--A-max", type=float, required=True)
    sp.add_argument("--grid-A", type=int, required=True, help="grid points along A")
    sp.add_argument("--steps", type=int, default=20)
    sp.add_argument("--C", type=float, required=True)
    sp.add_argument("--csv", type=Path)
    sp.add_argument("--json", type=Path)

    op = sub.add_parser("optimize-tail", help="Find minimal (p0,q0) so the SU(3) tail ≤ tol.")
    op.add_argument("--beta", type=float, required=True)
    op.add_argument("--tol", type=float, required=True)
    op.add_argument("--strategy", choices=["symmetric", "strip"], default="strip")
    op.add_argument("--N-min", type=int, default=0)
    op.add_argument("--N-max", type=int, default=40)
    op.add_argument("--mp-dps", type=int, default=80)
    op.add_argument("--json", type=Path)

    tg = sub.add_parser("target-gap", help="Compute τ0 to meet a target spectral gap m0.")
    tg.add_argument("--m0", type=float, required=True)

    st = sub.add_parser("string-table", help="Table of σ_phys and area-law bounds for a list of 'a'.")
    st.add_argument("--sigma-lat", type=float, required=True)
    st.add_argument("--a", type=float)
    st.add_argument("--a-list", type=str, help="Comma-separated list of lattice spacings.")
    st.add_argument("--area", type=float, default=1.0)
    st.add_argument("--csv", type=Path)
    st.add_argument("--json", type=Path)

    rc = sub.add_parser("recipe", help="Show curated example commands and what they are for.")

    # New: a-sweep
    asw = sub.add_parser("a-sweep", help="Uniform-in-a sweep: σ_phys, area-law bound, Initial Smallness (η0) & collar LB, and τ0⇒gap across spacings.")
    asw.add_argument("--beta", type=float, required=True, help="Base β in exp(-β C2/6).")
    asw.add_argument("--blocks", type=int, default=1, help="Number of RP-preserving heat-kernel convolutions.")
    asw.add_argument("--sym-N", type=int, default=8, help="Symmetric SU(3) shell cut N for tail control in η0 bound.")
    asw.add_argument("--A", type=float, required=True, help="Quadratic contraction constant A (for collar bound).")
    asw.add_argument("--C", type=float, required=True, help="Collar loss constant C (for collar bound).")
    asw.add_argument("--tau0", type=float, required=True, help="Tube-cost lower bound for the spectral gap (m0≥tau0).")
    # spacing selection
    asw.add_argument("--a-list", type=str, help="Comma-separated list of lattice spacings.")
    asw.add_argument("--a-min", type=float, help="Start of spacing range (use with --a-max and --a-num).")
    asw.add_argument("--a-max", type=float, help="End of spacing range (use with --a-min and --a-num).")
    asw.add_argument("--a-num", type=int, help="Number of points in spacing range (>=2).")
    asw.add_argument("--sigma-lat", type=float, required=True, help="Lattice string tension (assumed fixed as a varies).")
    asw.add_argument("--area", type=float, default=1.0, help="Wilson loop area for the area-law bound.")
    asw.add_argument("--mp-dps", type=int, default=120, help="mpmath precision (decimal digits) for η0/tails.")
    asw.add_argument("--csv", type=Path)
    asw.add_argument("--json", type=Path)
    asw.add_argument("--tex", type=Path)

    return p


def main() -> None:
    parser = _build()
    try:
        ns = parser.parse_args()
        if ns.cmd == "sweep-contraction":
            sys.exit(cmd_sweep_contraction(ns))
        elif ns.cmd == "optimize-tail":
            sys.exit(cmd_optimize_tail(ns))
        elif ns.cmd == "target-gap":
            sys.exit(cmd_target_gap(ns))
        elif ns.cmd == "string-table":
            sys.exit(cmd_string_table(ns))
        elif ns.cmd == "recipe":
            sys.exit(cmd_recipe(ns))
        elif ns.cmd == "a-sweep":
            sys.exit(cmd_a_sweep(ns))
        else:
            parser.error("Unknown command")
    except KeyboardInterrupt:
        print("\nAborted by user.", file=sys.stderr); sys.exit(130)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr); sys.exit(2)


if __name__ == "__main__":
    main()