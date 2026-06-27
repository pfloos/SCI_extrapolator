# 📊 Selected CI Extrapolation Tool

[![License](https://img.shields.io/badge/License-CC%20BY%20SA%204.0-lightgrey)](https://creativecommons.org/licenses/by/4.0/)
[![Last Update](https://img.shields.io/github/last-commit/pfloos/SCI_extrapolator?label=last%20update)](https://github.com/pfloos/SCI_extrapolator/commits/main)

[![GitHub Repo stars](https://img.shields.io/github/stars/pfloos/SCI_extrapolator?style=social)](https://github.com/pfloos/SCI_extrapolator/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/pfloos/SCI_extrapolator?style=social)](https://github.com/pfloos/SCI_extrapolator/network/members)
[![GitHub watchers](https://img.shields.io/github/watchers/pfloos/SCI_extrapolator?style=social)](https://github.com/pfloos/SCI_extrapolator/watchers)

---

## 📚 Table of Contents

- [✨ Key Features](#-key-features)
- [🎯 What This Tool Does](#-what-this-tool-does)
- [📥 Input Format](#-input-format)
- [⚙️ Installation](#️-installation)
- [🚀 Usage](#-usage)
- [📈 Output](#-output)
- [🔬 Method Details](#-method-details)
- [⚡ Features & Limitations](#-features--limitations)
- [🧠 Author](#-author)
- [📜 License](#-license)

---

## ✨ Key Features

- **🔮 Robust Extrapolation:**  
  Performs energy extrapolation to the **zero rPT2 limit** using weighted linear regression with automatic model selection.

- **🎯 Smart State Tracking:**  
  Matches electronic states consistently across iterations using a **fingerprint-based cost function** and the **Hungarian algorithm**.

- **📊 Publication-Ready Output:**  
  Generates formatted results tables with **spin labeling**, **uncertainties**, and **excitation energies** with propagated error bars.

- **🛡️ Stable & Reliable:**  
  Handles **state crossings**, automatically detects the number of states, and performs extrapolation with multiple fit windows for robustness.

- **⚡ Fast & Simple:**  
  No external dependencies beyond standard Python scientific libraries—just NumPy and SciPy.

---

## 🎯 What This Tool Does

For each electronic state in your SCI calculation, this script:

1. 🔍 **Extracts data** from a sequence of SCI iterations (variational energies and rPT2 corrections)
2. 🔗 **Tracks states** consistently across iterations using a sophisticated matching algorithm
3. 📈 **Extrapolates** the energy to the **zero rPT2 limit** using weighted linear regression:
   
   $$E_\text{var}(\text{rPT2}) \rightarrow E(\text{rPT2} = 0)$$

4. 📏 **Estimates uncertainties** from fit residuals
5. 🌟 **Computes excitation energies** relative to the ground state with **propagated error bars**
6. 🎲 **Assigns spin character** based on ⟨S²⟩ expectation value

---

## 📥 Input Format

The script expects a **JSON file** containing multiple SCI iterations. Each iteration must include:
- `n_det`: number of determinants
- `states`: list of electronic states

For each state:
- `energy`: variational energy (Hartree)
- `rpt2`: renormalized PT2 correction
- `s2`: ⟨S²⟩ value
- `variance`: variance of the wave function
- `ex_energy`: auxiliary state-dependent energies (used for matching)

### 📋 Example JSON Structure

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

---

## ⚙️ Installation

**No installation required!** Just ensure you have Python 3.6+ with standard scientific libraries:

```bash
pip install numpy scipy
```

---

## 🚀 Usage

### Basic Command

```bash
python3 SCI_extrapolator.py path/to/file.json
```

### Example

```bash
python3 SCI_extrapolator.py HF/aug-cc-pvdz/HF/json/00003.json
```

### Included example output

This repository also includes a sample run output file (human-readable log) at:

```
example/HF_aug-cc-pvtz.out
```

This file contains the script's run-time messages and a complete "Selected CI Extrapolation Summary" including the final extrapolated excitation energies table. You can view it on GitHub here:

https://github.com/pfloos/SCI_extrapolator/blob/main/example/HF_aug-cc-pvtz.out

Full results table from that run:

| State | Spin | ⟨S²⟩ | E_tot (Ha) | σ(E) (Ha) | ΔE (eV) | σ(ΔE) (eV) |
|-------:|:-----:|:-----:|-----------:|----------:|--------:|-----------:|
| 0 | Singlet | 0.0000 | -100.34944779 | 1.880e-04 | 0.000 | 0.007 |
| 1 | Triplet | 2.0000 | -99.97904441 | 1.386e-04 | 10.079 | 0.006 |
| 2 | Triplet | 2.0000 | -99.97920699 | 2.619e-04 | 10.075 | 0.009 |
| 3 | Singlet | 0.0000 | -99.96492204 | 9.306e-05 | 10.463 | 0.006 |
| 4 | Singlet | 0.0000 | -99.96510648 | 2.629e-04 | 10.458 | 0.009 |
| 5 | Triplet | 2.0000 | -99.85364321 | 2.075e-04 | 13.492 | 0.008 |
| 6 | Triplet | 2.0000 | -99.84794625 | 2.973e-04 | 13.647 | 0.010 |
| 7 | Triplet | 2.0000 | -99.84840278 | 3.136e-04 | 13.634 | 0.010 |
| 8 | Singlet | 0.0000 | -99.84317336 | 1.005e-04 | 13.776 | 0.006 |
| 9 | Singlet | 0.0000 | -99.84316692 | 5.240e-04 | 13.777 | 0.015 |
| 10 | Triplet | 2.0000 | -99.83249307 | 1.395e-04 | 14.067 | 0.006 |
| 11 | Triplet | 2.0000 | -99.82045549 | 1.668e-04 | 14.395 | 0.007 |

You can use this example to check formatting, understand the output layout, or as a test input when adapting the parser for different workflows.

### How to reproduce this example

To reproduce the included example output, run the script on the JSON file used for the run (named `00003.json` in the example). If you have the JSON in a folder such as `HF/.../json/00003.json`, run:

```bash
python3 SCI_extrapolator.py path/to/00003.json > example/HF_aug-cc-pvtz.out
```

Notes:
- The script prints human-readable logging messages; redirecting stdout (as shown) saves the full log to `example/HF_aug-cc-pvtz.out`.
- The script may attempt to auto-correct common JSON issues and will print a small diagnostic header (e.g. "Attempting to fix JSON file: 00003.json").
- If you prefer to inspect results interactively, omit the redirection and the summary table will be printed to the terminal.

---

## 📈 Output

The script produces **two structured sections**:

### 1️⃣ Run Summary

Metadata about your calculation:
- Input file name
- Number of iterations analyzed
- Number of tracked states
- Range of determinants
- Extrapolation model details

### 2️⃣ Final Results Table

A formatted table with all key results:

| Column | Description |
|--------|-------------|
| **#** | State index |
| **Spin** | Approximate spin multiplicity from ⟨S²⟩ |
| **⟨S²⟩** | Spin expectation value |
| **E_tot (Ha)** | Extrapolated total energy at rPT2 → 0 |
| **σ(E) (Ha)** | Uncertainty from linear fit |
| **ΔE (eV)** | Excitation energy relative to ground state |
| **σ(ΔE) (eV)** | Propagated uncertainty |

---

## 🔬 Method Details

### 🔗 State Tracking Algorithm

States are matched across iterations using a **weighted cost function** combining:

- Energy difference
- Spin contamination (⟨S²⟩)
- Wave function variance
- Excitation fingerprints (ex_energy)

The optimal assignment is solved using the **Hungarian algorithm** for maximum efficiency and robustness.

---

### 📐 Extrapolation Model

For each state, we fit a **linear model** in the small-rPT2 regime:

$$E(\text{rPT2}) = a \cdot \text{rPT2} + b$$

The extrapolated energy at the zero-rPT2 limit is:

$$E_0 = b \quad \text{at} \quad \text{rPT2} = 0$$

**Weighting scheme:**

$$w = \frac{1}{\text{rPT2}^2}$$

Fits automatically select **3–6 most reliable data points** from multiple fit windows to maximize stability.

---

### 📊 Error Estimates

Uncertainty is computed from **residuals of the weighted linear regression** and propagated for excitation energies as:

$$\sigma(\Delta E) = \sqrt{\sigma_i^2 + \sigma_0^2}$$

where $\sigma_i$ is the uncertainty of state $i$ and $\sigma_0$ is the uncertainty of the ground state.

---

## ⚡ Features & Limitations

### ✅ Strengths

- ✔️ Robust handling of **state crossings**
- ✔️ Automatic detection of **number of states**
- ✔️ Stable extrapolation with **multiple fit windows**
- ✔️ Spin-based state **labeling and identification**
- ✔️ Publication-quality **formatted output**

### ⚠️ Limitations

- Assumes **linear dependence** on rPT2 in the small-rPT2 regime
- Accuracy depends on quality of **last SCI iterations**
- Very noisy states may require **manual inspection**
- Best performance with **6+ iterations** in your dataset

---

## 🧠 Author

**[Pierre-François Loos](https://pfloos.github.io/WEB_LOOS)**

---

## 📜 License

See the repository for license information

---
