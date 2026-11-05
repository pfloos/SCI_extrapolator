# SCI extrapolator — Version 2 paired-state excitation extrapolations

This package builds on the corrected SCI extrapolator already discussed. It contains the complete replacement script, both optional Matplotlib plotting modules, the paired-state excitation models, the updated README, and regression tests. The original `example/` inputs and the `LICENSE` are intentionally **not** bundled: keep those from your Git checkout.

## Installation / update

1. In your local clone of `pfloos/SCI_extrapolator`, keep your Git working tree clean (`git status`).
2. Copy the **contents** of this folder into the repository root, allowing `SCI_extrapolator.py` and `README.md` to be replaced. Preserve the repository's `example/` and `LICENSE` folders/files.
3. On Python 3.10+, install `numpy`, `scipy`, and `matplotlib` (the last package is optional unless you plot):

   ```bash
   python3 -m pip install numpy scipy matplotlib
   ```

4. Run the tests and then make a PDF from the included HF example:

   ```bash
   python3 -m unittest discover -s tests -v
   python3 SCI_extrapolator.py example/00003.json --plot
   python3 SCI_extrapolator.py example/00003.json --excitation
   python3 SCI_extrapolator.py example/00003.json --plot-excitation
   ```

   This produces `example/00003_SCI_extrapolations.pdf`, containing a 12-state overview and 12 detailed pages.

5. Optional separate PNG exports and custom paths:

   ```bash
   python3 SCI_extrapolator.py example/00003.json --plot --plot-png
   python3 SCI_extrapolator.py example/00003.json --plot-output figures/HF.pdf --plot-dir figures/HF_png
   ```

6. Inspect your changes and commit them once satisfied:

   ```bash
   git status
   git add SCI_extrapolator.py sci_plotting.py sci_excitation.py sci_excitation_plotting.py README.md INSTALL.md tests/
   git commit -m "Add paired-state SCI excitation-energy extrapolations and plots"
   git push origin main
   ```

## Version 2 implementation

- `--excitation` prints three methods: independent total-energy fits (unchanged), direct paired variational-gap fits vs `q`, and paired rPT2-corrected gaps vs the same `q`, where `q` is the mean absolute correction.
- `--plot-excitation` produces a multipage comparison PDF. `--excitation-plot-png` enables PNGs.
- Both new fits use `q = (|rPT2_i| + |rPT2_0|)/2` for the horizontal coordinate, select windows by smallest `q`, and use residual weights proportional to `1/q`; a differential rPT2 of zero **does not** imply FCI convergence. The differential correction only enters the corrected gap's vertical coordinate.
- `sci_excitation.py` handles numerical analysis without Matplotlib; `sci_excitation_plotting.py` is imported only for plotting.
- Pairing uses exact iteration IDs, so missing historical states are not accidentally joined with another iteration.

## Version 1 implementation (unchanged)

- `--plot` writes a multipage PDF with one overview page plus a page for each state.
- `--plot-png` writes the overview and every individual state page to PNGs.
- Points are colored by determinant count (logarithmic scale), with black outlines on points used for the primary fit. The fit-window lines (3, 4, 5, 6 points) and extrapolated intercept are also shown.
- Every detailed page contains a plot of extrapolated-intercept shifts for the 3–6-point windows, in meV.
- Graphs zoom to the six points nearest the rPT2=0 limit, and say how many earlier points lie beyond the zoom. Points are not discarded from state tracking.
- `sci_plotting.py` loads Matplotlib only when plotting is requested; `numpy` and `scipy` remain sufficient for computations without plots.
- Plotting does **not** change the numerical extrapolations, weights, window selection, or uncertainty indicators.

As inspected on 5 November 2026, GitHub main contains the earlier bug corrections. This package builds on the Version 1 visualization add-on. These Version 2 files have **not** been pushed to GitHub.

## Validation

- The original 19 tests are retained; new tests additionally validate matched-iteration pairing, zero-correction cancellation, exact synthetic models, numerical invariance of the default result, and Version 2 PDF/PNG creation.
- The PDF preview provided alongside the package was generated from the stored **HF trajectories fixture** extracted from the original example, not by directly re-reading `example/00003.json` in this local environment. Run the commands above against your original checkout to reproduce it end to end.
