# `ym-research` — Research CLI for Exploring & Optimizing Mass-Gap Bounds

`ym-research` is a researcher-facing command-line tool that builds on the rigorous primitives in
[`ym-bounds`](../README.md). It helps you **sweep**, **optimize**, and **export** explicit,
machine-verifiable constants for the SU(3) mass-gap constructive pipeline.

- Explore safe contraction regions (η₀, A)
- Choose SU(3) tail cutoffs (p₀,q₀) for a target tolerance at fixed β
- Back-solve tube-cost for a desired spectral gap
- Generate string-tension tables across lattice spacings
- Export CSV/JSON for reproducibility

---

## Installation (editable dev)

From your repo root (with a Python 3.10+ venv recommended):

```bash
cd ym_bounds
pip install -e .
pip install mpmath
```
This installs both CLIs:
	•	ym-bounds — prints a single, self-contained constants report
	•	ym-research — exploration & optimization workflows (this tool)

⸻

Quick Start

List recipes (what to run and why):

ym-research recipe

Run a full constants snapshot (via ym-bounds) and export artifacts:
```
ym-research sweep-contraction \
  --eta0-min 0.04 --eta0-max 0.06 --grid-eta 9 \
  --A-min 2.5  --A-max 3.5  --grid-A 11 \
  --steps 20 --C 0.2 \
  --csv src/ym_research/bounds.csv \
  --json src/ym_research/bounds.json
```

⸻

Perfect — since your test run of

ym-research sweep-contraction ...

succeeded, we now know exactly how the CLI is supposed to be used. Let’s re-check the rest of the commands you pasted and adjust them so they’re consistent with the working syntax.

⸻

Command Set

1) sweep-contraction


Here’s the dense version they suggest in the “workflows” section:
```
ym-research sweep-contraction \
  --eta0-min 0.02 --eta0-max 0.08 --grid-eta 13 \
  --A-min 2.0  --A-max 5.0  --grid-A 13 \
  --steps 20 --C 0.2 \
  --csv src/ym_research/sweep_dense.csv
```

⸻

2) optimize-tail

Your tool expects flags: --beta, --tol, --strategy, --N-min, --N-max, and optionally --json.

Strip strategy example:

ym-research optimize-tail \
  --beta 6.0 --tol 1e-2 \
  --strategy strip --N-min 4 --N-max 20 \
  --json src/ym_research/tail.json

Symmetric strategy example:

ym-research optimize-tail \
  --beta 6.0 --tol 5e-3 \
  --strategy symmetric --N-min 4 --N-max 30


⸻

3) target-gap

This one is simple and  as-is:
```
ym-research target-gap --m0 0.5

```
⸻

4) string-table

This one was failing earlier because of --input (not valid). The right flags are: --sigma-lat, --a or --a-list, --area, and output flags (--csv, --json, --tex).

Single spacing:
```
ym-research string-table \
  --sigma-lat 0.045 --a 0.08 --area 1 \
  --csv src/ym_research/sigma_single.csv
```
Multiple spacings:
```
ym-research string-table \
  --sigma-lat 0.045 --a-list 0.12,0.10,0.08 --area 1 \
  --csv src/ym_research/sigma_multi.csv
```

⸻

5) recipe

This is just:
```
ym-research recipe
```

⸻


Reproducibility
	•	Outputs are deterministic for given inputs (no randomness).
	•	ym-bounds report also records a metadata header (UTC time, Python version).
	•	Recommended: commit the CSV/JSON artifacts alongside your manuscript for auditability.

⸻

Theory Notes (what’s under the hood)
	•	Contraction & Collar: uses explicit quadratic contraction ( \eta_{k+1} \le A\eta_k^2 ) and
stable computation of ( \prod_k (1 - C \eta_k) ).
	•	SU(3) Tail: rigorous shell-wise bound with the exact shell minimum
( C_2^{\min}(n) = n^2/4 + n ); polynomial envelope for shell multiplicities; integral over-bound.
	•	Tube-Cost ⇒ Gap: certifies Spec(T) ⊂ {1} ∪ [e^{−τ₀},1) and hence m₀ ≥ τ₀.
	•	String Tension: ( \sigma_{\rm phys} = \sigma_{\rm lat} / a^2 ); area-law upper bound ( e^{-\sigma_{\rm phys} \cdot \text{Area}} ).
	•	Clustering: exports a bound of the form ( C e^{-m^\ast |x|} ) (you choose (m^\ast)).

⸻

Example Session

# Explore recipes
ym-research recipe

# Sweep contraction region
ym-research sweep-contraction \
  --eta0-min 0.02 --eta0-max 0.08 --grid-eta 7 \
  --A-min 2.5 --A-max 4.5 --grid-A 9 \
  --steps 20 --C 0.2 --csv sweep.csv --json sweep.json

# Optimize SU(3) tail
ym-research optimize-tail --beta 6.0 --tol 1e-2 --strategy strip --N-min 4 --N-max 20 --json tail.json

# Target a 0.5 gap
ym-research target-gap --m0 0.5

# Build a string table
ym-research string-table --sigma-lat 0.045 --a-list 0.12,0.10,0.08 --area 1 --csv sigma.csv

# Freeze a paper-ready constants snapshot
ym-bounds report \
  --beta 6.0 --p0 6 --q0 0 \
  --eta0 0.05 --A 3.0 --C 0.2 --steps 20 \
  --tau0 0.4 --sigma-lat 0.045 --a 0.08 --area 1 \
  --csv bounds.csv --tex bounds.tex --json bounds.json


⸻

Troubleshooting
	•	“command not found” — ensure your venv is active and package installed: python -m pip install -e .
	•	ImportError for console script — check pyproject.toml has:

[project.scripts]
ym-research = "ym_research.cli:main"
ym-bounds = "ym_bounds.cli:main"


	•	Mac zsh ignores comments — don’t paste shell comments (# ...) directly; remove them or run commands line by line.

⸻

License

MIT (see repository LICENSE). Contributions welcome—please include tests for new commands or bound refinements.

⸻


