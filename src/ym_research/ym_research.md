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
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .

This installs both CLIs:
	•	ym-bounds — prints a single, self-contained constants report
	•	ym-research — exploration & optimization workflows (this tool)

⸻

Quick Start

List recipes (what to run and why):

ym-research recipe

Run a full constants snapshot (via ym-bounds) and export artifacts:

ym-bounds report \
  --beta 6.0 --p0 6 --q0 0 \
  --eta0 0.05 --A 3.0 --C 0.2 --steps 20 \
  --tau0 0.4 --sigma-lat 0.045 --a 0.08 --area 1.0 \
  --csv bounds.csv --tex bounds.tex --json bounds.json


⸻

Commands

1) sweep-contraction

Map safe regions of initial smallness η0 and contraction constant A where:
	•	the correction tail Σ η_k is finite; and
	•	the collar product ∏(1 − C η_k) stays positive.

ym-research sweep-contraction \
  --eta0-min 0.02 --eta0-max 0.08 --grid-eta 7 \
  --A-min 2.0 --A-max 5.0 --grid-A 7 \
  --steps 20 --C 0.2 \
  --csv sweep.csv --json sweep.json

Output columns (CSV):
	•	eta0, A, steps, C, sum_eta, product, sum_finite, product_positive

Interpretation:
	•	sum_finite==True + product_positive==True → region satisfies the contraction and collar constraints.

⸻

2) optimize-tail

Pick (p₀,q₀) so the SU(3) rep tail is below a given tolerance at fixed β.

Two strategies:
	•	strip: keeps q0=0, increases p0
	•	symmetric: increases p0=q0 together

Examples

# Fundamental strip (q0=0), target tail ≤ 1e-2 at β=6
ym-research optimize-tail --beta 6.0 --tol 1e-2 \
  --strategy strip --N-min 4 --N-max 20 --json tail.json

# Symmetric (p0=q0), tighter tolerance
ym-research optimize-tail --beta 6.0 --tol 5e-3 \
  --strategy symmetric --N-min 4 --N-max 30

Interpretation:
	•	The tool returns the smallest N that meets tail ≤ tol, along with the partial and tail values.
	•	Uses a sharp shell minimum: ( \min_{p+q=n} C_2 = \frac{n^2}{4} + n ), ensuring a rigorous bound.

⸻

3) target-gap

Back out τ₀ from a desired spectral gap m₀ using the certified relation m₀ ≥ τ₀.

ym-research target-gap --m0 0.5

Interpretation:
	•	Prints Spec(T) ⊂ {1} ∪ [e^{−τ₀}, 1) with τ₀ = m₀_target.
	•	Use this to align spectral claims with tube-cost lower bounds in your text.

⸻

4) string-table

Convert lattice string tension to physical units across one or multiple spacings a,
and compute the area-law upper bound for a given loop area.

# Single spacing
ym-research string-table --sigma-lat 0.045 --a 0.08 --area 1 --csv sigma_single.csv

# Multiple spacings
ym-research string-table --sigma-lat 0.045 --a-list 0.12,0.10,0.08 --area 1 \
  --csv sigma_multi.csv

Output columns (CSV):
	•	sigma_lat, a, sigma_phys, area, wl_bound

Interpretation:
	•	sigma_phys = sigma_lat / a^2
	•	wl_bound = exp(−sigma_phys * area) (area-law upper bound for ⟨W(C)⟩)

⸻

5) recipe

Show curated example commands and what they’re for:

ym-research recipe


⸻

Practical Workflows

A. Certify a constants snapshot (paper appendix)
	1.	Choose parameters (β, p₀,q₀, η₀, A, C, τ₀, σ_lat, a).
	2.	Run:

ym-bounds report \
  --beta 6.0 --p0 6 --q0 0 \
  --eta0 0.05 --A 3.0 --C 0.2 \
  --tau0 0.4 --sigma-lat 0.045 --a 0.08 --area 1 \
  --csv bounds.csv --tex bounds.tex --json bounds.json


	3.	Commit bounds.csv and bounds.json for reproducibility; paste bounds.tex into the appendix.

B. Find a safe RG region (summability + collar positivity)

ym-research sweep-contraction \
  --eta0-min 0.02 --eta0-max 0.08 --grid-eta 13 \
  --A-min 2.0 --A-max 5.0 --grid-A 13 \
  --steps 20 --C 0.2 \
  --csv sweep_dense.csv

Sort by product descending to highlight robust zones.

C. Tighten SU(3) representation tails

ym-research optimize-tail --beta 6.0 --tol 1e-3 --strategy strip --N-min 4 --N-max 40

Use the reported N as (p0=N,q0=0) in your ym-bounds run.

⸻

Input Validation & Failure Modes
	•	All commands validate positivity/finite-ness of inputs and fail fast with clear errors.
	•	If CSV/JSON paths point to non-existent directories, they will be created.
	•	For string-table, use either --a or --a-list. If both are omitted, the command fails.

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


