# SCI_extrapolator correction bundle — 5 November 2026

This bundle replaces the script and README, and adds unit tests. It does **not** contain the original Quantum Package example JSON or license. The GitHub connector did not permit creating a branch/commit (403), so the upstream repository was not changed.

## Applying to an existing checkout

1. Back up your existing `SCI_extrapolator.py` and `README.md` (or use Git to preserve the previous commit).
2. Copy the bundle's `SCI_extrapolator.py` and `README.md` to the repository root.
3. Copy the entire bundle's `tests/` folder to the repository root.
4. Install dependencies if necessary: `python3 -m pip install numpy scipy` (Python 3.10+).
5. Run `python3 -m unittest discover -s tests -v`.
6. Run `python3 SCI_extrapolator.py example/00003.json` on the original example input (it remains in your GitHub checkout).
7. When satisfied, regenerate the saved example log using
   `python3 SCI_extrapolator.py example/00003.json > example/HF_aug-cc-pvtz.out`.

## Corrections

- Adds the largest SCI expansion to every state's trajectory.
- Preserves original root IDs across rectangular Hungarian assignments.
- Uses an effective weighted least-squares objective with residuals divided by |rPT2|, rather than |rPT2|².
- Uses a fixed five-point primary window (or fewer if necessary) and compares 3–6-point intercepts.
- Reports residual-based intercept uncertainty and fit-window sensitivity separately. The final `u(E)` indicator is the maximum of the two; it is **not** a statistical confidence interval.
- Handles rPT2=0 and singular/rank-deficient windows without unhandled numerical crashes.
- Makes the ground-state excitation energy and its indicator exactly zero.
- Validates input and repairs only unclosed JSON container suffixes, without modifying the original.
- Fixes Python-version requirements, command-line description, script filename and license link in README.
- Adds 15 tests, including a regression check based on trajectories extracted from the HF example and a synthetic end-to-end JSON test.

## Scope / limitations

The test suite passed with NumPy 2.3.5 and SciPy 1.17.0. The full original 47 KB `example/00003.json` was inspected via GitHub but was not available inside the local execution environment. Therefore, the HF regression is based on the matched `(rPT2, energy)` trajectories extracted from that file, not a direct end-to-end run of the unmodified example file. You should run step 6 before merging.

State matching is still heuristic (not based on wavefunction overlaps). Excitation-energy stability indicators ignore possible correlations between states. The statistical interpretation of rPT2-based weights and rigorous SCI extrapolation errors remains a methodological question.
