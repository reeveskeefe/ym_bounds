#!/usr/bin/env python3
"""
main.py - Unified entry point for ym_bounds pipeline

This script provides the main interface described in the proof pipeline:
    python main.py --beta 6 --constants input_constants.json

It's a convenience wrapper around the full CLI in src/ym_bounds/cli.py
"""

import sys
import argparse
from pathlib import Path

# Add src to path to import ym_bounds
sys.path.insert(0, str(Path(__file__).parent / "src"))

from ym_bounds.cli import main as cli_main


def main():
    """Main entry point for unified proof pipeline"""
    parser = argparse.ArgumentParser(
        description="Yang-Mills mass gap bounds computation - unified pipeline interface",
        epilog="This is a convenience wrapper. For full options, use: python -m ym_bounds.cli"
    )
    
    # Core parameters
    parser.add_argument("--beta", type=float, required=True, help="Coupling constant β")
    parser.add_argument("--constants", type=Path, help="JSON file with constants from RG step analysis")
    
    # Optional overrides
    parser.add_argument("--eta0", type=float, help="Override initial KP norm")
    parser.add_argument("--A", type=float, help="Override contraction constant A")
    parser.add_argument("--C", type=float, help="Override collar constant C") 
    parser.add_argument("--tau0", type=float, help="Override tube cost τ₀")
    
    # Physical parameters (required)
    parser.add_argument("--sigma-lat", type=float, default=0.045, help="Lattice string tension")
    parser.add_argument("--a", type=float, default=0.08, help="Lattice spacing")
    
    # Output options
    parser.add_argument("--csv", type=Path, help="Export results to CSV")
    parser.add_argument("--tex", type=Path, help="Export results to TeX")
    parser.add_argument("--json", type=Path, help="Export results to JSON")
    
    # Advanced options
    parser.add_argument("--steps", type=int, default=20, help="RG steps for computation")
    parser.add_argument("--sym-N", type=int, default=8, help="Symmetric shell cut")
    parser.add_argument("--mp-dps", type=int, default=120, help="Multiprecision decimal places")
    
    args = parser.parse_args()
    
    # Build arguments for the full CLI
    cli_args = [
        "report",
        "--beta", str(args.beta),
        "--sigma-lat", str(args.sigma_lat),
        "--a", str(args.a),
        "--steps", str(args.steps),
        "--sym-N", str(args.sym_N),
        "--mp-dps", str(args.mp_dps),
    ]
    
    # Add constants file if provided
    if args.constants:
        cli_args.extend(["--constants", str(args.constants)])
    
    # Add any explicit overrides
    if args.eta0 is not None:
        cli_args.extend(["--eta0", str(args.eta0)])
    if args.A is not None:
        cli_args.extend(["--A", str(args.A)])
    if args.C is not None:
        cli_args.extend(["--C", str(args.C)])
    if args.tau0 is not None:
        cli_args.extend(["--tau0", str(args.tau0)])
    
    # Add output options
    if args.csv:
        cli_args.extend(["--csv", str(args.csv)])
    if args.tex:
        cli_args.extend(["--tex", str(args.tex)])
    if args.json:
        cli_args.extend(["--json", str(args.json)])
    
    print("🔗 Yang-Mills Mass Gap Proof Pipeline")
    print(f"   Running: ym_bounds {' '.join(cli_args)}")
    print()
    
    # Replace sys.argv and call the main CLI
    original_argv = sys.argv[:]
    try:
        sys.argv = ["ym_bounds"] + cli_args
        return cli_main()
    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    sys.exit(main())
