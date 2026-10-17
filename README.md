# ContourRepair

Reference implementation for **Contour Quadrature with Missing Evaluations: Sharp Annular Recovery and the Cost of One Erasure**, by Shivam Gupta.

Read the [paper](paper/main.pdf). The method recovers a circular contour average from an incomplete equally spaced grid under a specified Laurent coefficient-energy assumption. It uses FFTs and a solve whose dimension equals the number of missing evaluations. It returns weights, the sharp unit-ball error, a normal-equation diagnostic, and exact-arithmetic one-retry gains.

The general optimal recovery and kernel quadrature principles are classical. The paper develops an explicit one-erasure formula, its convergence exponent and lower bound, and a complementary computation for incomplete contour information. See [the literature record](LITERATURE.md) for scope and prior work.

## Reproduce

Python 3.9 or newer is required. The recorded run used Python 3.9.6 on an Apple M4 Pro; exact package versions are in `requirements.txt` and `results/environment.json`.

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 run_experiments.py
```

This regenerates the CSV/JSON evidence under `results/` and four publication figures under `paper/figures/`. Timing values depend on the host; the committed timing file records the manuscript run. No downloads, remote services, paid APIs, or human data are used by the experiments.

To compile the manuscript with [Tectonic](https://tectonic-typesetting.github.io/):

```sh
cd paper
tectonic main.tex
```

## Example

```python
import numpy as np
from src.contour_recovery import recover, observed_error_bound

n, q = 64, 0.8
missing = [0, 1, 2, 3]
rule = recover(n, q, missing)
z = np.exp(2j * np.pi * np.arange(n) / n)
available = np.array([j for j in range(n) if j not in missing])
values = 1 / (z[available] - 0.4) + 0.2
estimate = rule.weights[available] @ values
print(estimate, rule.error, rule.gains)
# A bound on the actual function's Laurent norm must be supplied separately.
M = np.sqrt(0.2**2 + 1 / (q - 0.4**2))
print(observed_error_bound(rule.weights, q, M))
```

Use only available indices when applying weights: multiplying a zero weight by a stored NaN still produces NaN. The code never needs values at erased nodes.

## Evidence and limits

- All 256 masks at N=8 are compared against a separate dense kernel implementation. Every one-retry gain is checked by restoring its node.
- Exact rational arithmetic checks the four-node failure of diminishing returns.
- A residual extremizer attains the error bound and vanishes at the available nodes.
- Forty rational functions and three mask patterns yield 600 method records; a 48-dimensional manufactured nonnormal matrix gives a projector-action test.
- At fixed q, one erasure halves the asymptotic minimax convergence exponent. This is a function-class worst-case result, not a claim about every integrand.
- Some Fourier baselines are more accurate on the smoother nominal tests; these outcomes are retained in the paper.
- The floating-point output evaluates rigorous error formulas but is not an outward-rounded numerical enclosure. The analytic scale and amplitude bound must be justified. The implementation rejects a full-grid spectral condition number above 1e12.
- The noise contribution is a sufficient additive bound. The implementation does not optimize weights for noisy observations or certify a complete eigensolver.

`COMMERCIAL_CONTEXT.md` distinguishes an integration hypothesis from demonstrated market demand. No acceptance, patentability, or financial claim is made.

## Layout

| Path | Purpose |
|---|---|
| `src/contour_recovery.py` | FFT recovery, independent dense comparator, bounds, retry selection |
| `tests/test_recovery.py` | Independent numerical and analytic verification |
| `run_experiments.py` | All experiments and plot generation, with fixed seeds |
| `results/` | Raw data, environment and exact counterexample |
| `paper/` | TeX source, PDF and four figures |

Code: MIT license. Manuscript and figures: copyright Shivam Gupta, 2026; see `LICENSE-MANUSCRIPT`.
