# 📊 Selected CI Extrapolation Tool

[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey)](https://creativecommons.org/licenses/by-sa/4.0/)
[![Last Update](https://img.shields.io/github/last-commit/pfloos/SCI_extrapolator?label=last%20update)](https://github.com/pfloos/SCI_extrapolator/commits/main)

A small Python tool to extrapolate selected configuration-interaction (SCI/CIPSI) energies to the **full-CI limit within a fixed one-electron basis**, using Quantum Package rPT2 data. It also estimates excitation energies and tracks electronic states between iterations.

## Features

- Linear extrapolation of variational energy against rPT2, emphasizing points near rPT2 = 0.
- Backward state tracking from the largest CI expansion using a fingerprint cost and the Hungarian algorithm.
- Configurable compact output and a detailed diagnostic table.
- A **heuristic extrapolation-stability indicator**, incorporating the sensitivity of the fitted intercept to different fitting windows.
- Spin labels inferred from the ⟨S²⟩ values of the largest expansion.
- Conservative recovery of JSON files truncated at the end, without overwriting the input.

**Important:** State tracking is heuristic rather than overlap-based. Near-degenerate same-spin states can be difficult to distinguish. The reported `u` values are not calibrated statistical standard deviations or confidence intervals.

## Installation

Requires **Python 3.10+**, NumPy and SciPy:

```bash
python3 -m pip install numpy scipy
```

No package installation is needed; the script can be run directly.

## Usage

```bash
python3 SCI_extrapolator.py example/00003.json
python3 SCI_extrapolator.py --compact example/00003.json
python3 SCI_extrapolator.py --compact --diagnostics example/00003.json
python3 SCI_extrapolator.py --no-colors example/00003.json
```

The standard output includes a run summary, an extrapolated energy table and fit diagnostics. Compact mode prints the summary and a shorter energy table; `--diagnostics` also shows diagnostics in compact mode. Exit status is nonzero for invalid data or failed numerical processing.

## Input format

The JSON file should contain SCI iterations, commonly under the `fci` key. Quantum Package 2 typically generates these files under the calculation's `ezfio/json/` directory. See the [Quantum Package 2 repository](https://github.com/QuantumPackage/qp2).

Each iteration requires:

- `n_det`: positive integer number of determinants.
- `states`: nonempty list of state objects.

Each state requires:

- `energy`: variational energy in Hartree.
- `rpt2`: renormalized second-order perturbative correction in Hartree.
- `s2`: expectation value of S².
- `variance`: wavefunction variance.

The optional `ex_energy` list is used in state matching. Other Quantum Package fields, including `pt2`, `pt2_err`, `rpt2_err` and `variance_err`, are allowed but not currently used by the extrapolation fit.

Minimal example:

```json
{
  "fci": [
    {
      "n_det": 1000,
      "states": [
        {"energy": -100.10, "rpt2": -0.05, "s2": 0.0, "variance": 0.1, "ex_energy": []}
      ]
    }
  ]
}
```

At least **three usable points spanning sufficiently distinct rPT2 values** per state are required for extrapolation (unless the correction is exactly zero). The JSON may contain other metadata. A valid input is never modified. If a file is missing only closing containers at its end, a corrected copy is written as `<input>.fixed`; malformed internal content and unfinished strings are rejected rather than silently altered.

## Method

### State matching

Starting from the largest CI expansion, roots are matched backward through successive iterations using the Hungarian algorithm. The cost includes differences in variational energy, S², variance, and the shared part of the `ex_energy` fingerprint:

\[
C_{ij}=|E_i-E_j|+0.1|\langle S^2\rangle_i-\langle S^2\rangle_j|+0.1|v_i-v_j|+0.5\lVert\mathbf f_i-\mathbf f_j\rVert.
\]

Only the common length of the two fingerprint vectors is compared. State IDs are preserved even when the number of computed roots changes; a root that cannot be matched at an earlier iteration is not assigned to a different state's ID. **The largest wavefunction is always included** in the trajectory. The assignment is a matching heuristic and cannot guarantee correct identity across arbitrary avoided crossings.

### Extrapolation

For each state, fit

\[
E_{\mathrm{var}}(x)=a x+b,\qquad x=\mathrm{rPT2},\qquad E_{\mathrm{FCI}}\approx b.
\]

We minimize the weighted residual sum

\[
\sum_k\left[\frac{E_k-(a x_k+b)}{|x_k|}\right]^2.
\]

This corresponds to `numpy.polyfit(..., w=1/abs(x))`, **not** `w=1/x²`. The weighting is a heuristic for emphasizing the most converged expansions, not a measured inverse measurement variance. Near-duplicate or rank-deficient fit windows are rejected.

- **Primary fit:** five points with the smallest absolute rPT2 corrections, or three/four if fewer are available.
- **Diagnostics:** independent linear fits using the 3-, 4-, 5- and 6-point windows, when valid.
- **Window sensitivity:** largest absolute deviation of a diagnostic-window intercept from the primary intercept.
- **Residual-based intercept estimate:** calculated from the primary fit's weighted residual covariance, conditional on the heuristic weighting scheme.

The reported stability indicator is

\[
u(E)=\max\!\left(u_{\mathrm{fit}},\ \max_{k\in\{3,4,5,6\}}|b_k-b_{\mathrm{primary}}|\right).
\]

This is a **practical sensitivity indicator**, not a statistically calibrated one-sigma uncertainty. In particular, the different windows are correlated and are not independent observations. If an input rPT2 is *exactly* zero, its variational energy is used directly (and `u(E)=0`); this assumes that the zero correction is physically meaningful, rather than rounded or otherwise invalid data.

### Excitation energies

\[
\Delta E_i=(E_i-E_0)\times27.211386245988\ \mathrm{eV/Ha}.
\]

For the displayed heuristic indicator, the default propagation is

\[
u(\Delta E_i)=\sqrt{u(E_i)^2+u(E_0)^2}\times27.211386245988.
\]

This assumes independent energy errors and does not estimate correlations between state extrapolations. The ground-state excitation energy and its indicator are **exactly zero by definition**.

## Included HF example

The input file is [`example/00003.json`](example/00003.json), containing 9 SCI iterations and 12 electronic states. The table below is representative of the **corrected** primary five-point fits, including the last SCI iteration:

| State | Spin | E_tot (Ha) | u(E) (Ha) | ΔE (eV) | u(ΔE) (eV) |
|---:|:---|---:|---:|---:|---:|
| 0 | Singlet | -100.34935819 | 2.875e-04 | 0.000 | 0.000 |
| 1 | Triplet | -99.97919020 | 4.064e-04 | 10.073 | 0.014 |
| 2 | Triplet | -99.97915906 | 4.701e-04 | 10.074 | 0.015 |
| 3 | Singlet | -99.96500885 | 3.672e-04 | 10.459 | 0.013 |
| 4 | Singlet | -99.96507457 | 4.830e-04 | 10.457 | 0.015 |
| 5 | Triplet | -99.85388342 | 3.388e-04 | 13.483 | 0.012 |
| 6 | Triplet | -99.84827297 | 5.454e-04 | 13.635 | 0.017 |
| 7 | Triplet | -99.84835780 | 5.542e-04 | 13.633 | 0.017 |
| 8 | Singlet | -99.84342696 | 4.798e-04 | 13.767 | 0.015 |
| 9 | Singlet | -99.84340936 | 6.760e-04 | 13.768 | 0.020 |
| 10 | Triplet | -99.83248652 | 4.171e-04 | 14.065 | 0.014 |
| 11 | Triplet | -99.82047402 | 2.311e-04 | 14.392 | 0.010 |

The repository's original [`example/HF_aug-cc-pvtz.out`](example/HF_aug-cc-pvtz.out) was produced by the **previous** implementation, which omitted the final iteration and selected three-point fits by minimum residual; those older values should not be used as a regression reference. Regenerate the log with:

```bash
python3 SCI_extrapolator.py example/00003.json > example/HF_aug-cc-pvtz.out
```

## Tests

Run the included unit tests (no pytest dependency required):

```bash
python3 -m unittest discover -s tests -v
```

The tests cover state tracking with missing roots, inclusion of the final iteration, weighting, fit-window behavior, rank deficiency, zero rPT2, input validation, JSON recovery, and the HF benchmark trajectories.

## Author and license

[**Pierre-François Loos**](https://pfloos.github.io/WEB_LOOS) — [CC BY-SA 4.0](LICENSE).
