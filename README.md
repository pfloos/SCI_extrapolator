# Selected CI Extrapolation Tool

This script performs robust extrapolation of **selected CI (sCI) energies** using the renormalized second-order perturbative correction (rPT2). It is designed to process JSON outputs from quantum chemistry calculations and compute reliable excitation energies with associated uncertainty estimates.

---

## Overview

For each electronic state, the script:
1. Extracts **variational energies** and **rPT2 corrections** from a sequence of SCI iterations.
2. Tracks states consistently across iterations using a fingerprint-based matching procedure.
3. Extrapolates the energy to the **zero rPT2 limit** using weighted linear regression:
   $ E_\text{var}(rPT2) \rightarrow E(rPT2 = 0) $
4. Estimates uncertainties from fit residuals.
5. Computes **excitation energies relative to the ground state**, including propagated error bars.
6. Assigns approximate **spin character** based on ⟨S²⟩.

---

## Input Format

The script expects a JSON file containing multiple SCI iterations. Each iteration includes:
- `n_det`: number of determinants
- `states`: list of electronic states
- For each state:
  - `energy`: variational energy (Hartree)
  - `rpt2`: renormalized PT2 correction
  - `s2`: ⟨S²⟩ value
  - `variance`: variance of the wave function
  - `ex_energy`: auxiliary state-dependent energies (used for matching)
Example structure:
```json
{
  "fci": [
    {
      "n_det": 12345,
      "states": [
        {
          "energy": -100.123,
          "rpt2": -0.0012,
          "s2": 0.0,
          "variance": 0.01,
          "ex_energy": [...]
        }
      ]
    }
  ]
}
```

⸻

## Installation

No installation required beyond standard Python scientific libraries:

```bash
pip install numpy scipy
```

⸻

## Usage

Run the script directly on a JSON file:

```bash
python3 extrap.py path/to/file.json
```

Example:

```bash
python3 extrap.py HF/aug-cc-pvdz/HF/json/00003.json
```

⸻

## Output

The script prints two structured sections.

1. Run Summary

Provides metadata about the calculation:

* file name
* number of iterations
* number of tracked states
* determinant range
* extrapolation model details

⸻

2. Final Results Table

A formatted table containing:

Column	Description
#	State index
Spin	Approximate spin multiplicity from ⟨S²⟩
⟨S²⟩	Spin expectation value
E_tot (Ha)	Extrapolated total energy at rPT2 → 0
σ(E) (Ha)	Uncertainty from linear fit
ΔE (eV)	Excitation energy relative to ground state
σ(ΔE) (eV)	Propagated uncertainty

⸻

## Method Details

State Tracking

States are matched across iterations using a cost function combining:

* energy difference
* spin contamination (⟨S²⟩)
* variance
* excitation fingerprints (ex_energy)

Assignment is solved using the Hungarian algorithm.

⸻

## Extrapolation Model

For each state:

[
E(rPT2) = a \cdot rPT2 + b
]

The extrapolated energy is:

[
E_0 = b \quad \text{at} \quad rPT2 = 0
]

Weights used in the fit:

[
w = \frac{1}{rPT2^2}
]

Fits are performed using 3–6 most reliable data points.

⸻

## Error Estimate

Uncertainty is estimated from the residuals of the weighted linear regression and propagated for excitation energies as:

[
\sigma(\Delta E) = \sqrt{\sigma_i^2 + \sigma_0^2}
]

⸻

## Features

* Robust handling of state crossings
* Automatic detection of number of states
* Stable extrapolation with multiple fit windows
* Spin labeling based on ⟨S²⟩
* Publication-quality formatted output table

⸻

## Limitations

* Assumes linear dependence on rPT2 in the small-rPT2 regime
* Accuracy depends on quality of last SCI iterations
* Very noisy states may require manual inspection

⸻

## Author

Developed for selected CI post-processing workflows in quantum chemistry.

⸻

## License

Internal / academic use (adapt as needed).

---
