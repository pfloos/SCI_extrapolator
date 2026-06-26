# Selected CI Extrapolation Tool

This script performs robust extrapolation of **selected CI (sCI) energies** using the renormalized second-order perturbative correction (rPT2). It is designed to process JSON outputs from quantum chemistry calculations and compute reliable excitation energies with associated uncertainty estimates.

---

## Overview

For each electronic state, the script:

1. Extracts **variational energies** and **rPT2 corrections** from a sequence of SCI iterations.
2. Tracks states consistently across iterations using a fingerprint-based matching procedure.
3. Extrapolates the energy to the **zero rPT2 limit** using weighted linear regression:
   \[
   E_\text{var}(rPT2) \rightarrow E(rPT2 = 0)
   \]
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
