#!/usr/bin/env python3
"""
SCI Extrapolator
Extrapolate Selected CI energies to the complete basis set limit using rPT2.
"""

import sys
import json
import argparse
import numpy as np
from scipy.optimize import linear_sum_assignment

# Constants
HA_TO_EV = 27.211386245988
SPIN_REFERENCES = {
    "Singlet": 0.0,
    "Doublet": 0.75,
    "Triplet": 2.0,
    "Quartet": 3.75,
    "Quintet": 6.0
}

# ============================================================
# Data Extraction
# ============================================================

def extract_iterations(obj: dict | list) -> list[dict]:
    """
    Recursively extract all dictionaries containing 'states' and 'n_det'.
    
    Args:
        obj: Dictionary or list to search through
        
    Returns:
        List of iteration dictionaries
    """
    out = []
    if isinstance(obj, dict):
        if "states" in obj and "n_det" in obj:
            out.append(obj)
        for v in obj.values():
            out.extend(extract_iterations(v))
    elif isinstance(obj, list):
        for v in obj:
            out.extend(extract_iterations(v))
    return out


def load_data(filename: str) -> list[dict]:
    """
    Load and parse JSON data, returning sorted CI iterations.
    
    Args:
        filename: Path to JSON file
        
    Returns:
        List of iterations sorted by number of determinants
    """
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error loading file: {e}")
        sys.exit(1)

    iters = extract_iterations(data)
    if not iters:
        print("No valid CI iterations found in the file.")
        sys.exit(1)
    
    # Sort by number of determinants (smallest to largest)
    iters.sort(key=lambda x: x.get("n_det", 0))
    return iters


# ============================================================
# Fingerprint for State Matching
# ============================================================

def fingerprint(state: dict) -> dict:
    """
    Create a fingerprint dictionary for a single state.
    
    Args:
        state: State dictionary from JSON
        
    Returns:
        Dictionary containing energy, S2, variance, and excitation energies
    """
    return {
        "energy": state.get("energy", 0.0),
        "s2": state.get("s2", 0.0),
        "variance": state.get("variance", 0.0),
        "ex": np.array(state.get("ex_energy", []), dtype=float),
    }


def distance(a: dict, b: dict) -> float:
    """
    Calculate distance metric between two state fingerprints.
    
    Args:
        a: First fingerprint
        b: Second fingerprint
        
    Returns:
        Distance metric (lower = more similar)
    """
    energy_diff = abs(a["energy"] - b["energy"])
    s2_diff = abs(a["s2"] - b["s2"])
    variance_diff = abs(a["variance"] - b["variance"])

    # Compare excitation energies up to the minimum length
    ex_diff = 0.0
    if len(a["ex"]) and len(b["ex"]):
        min_len = min(len(a["ex"]), len(b["ex"]))
        ex_diff = np.linalg.norm(a["ex"][:min_len] - b["ex"][:min_len])

    return energy_diff + 0.1 * s2_diff + 0.1 * variance_diff + 0.5 * ex_diff


# ============================================================
# Spin Label
# ============================================================

def spin_label(s2: float) -> str:
    """
    Return the spin label closest to the given S2 value.
    
    Args:
        s2: S2 value
        
    Returns:
        Spin label (Singlet, Doublet, etc.)
    """
    return min(SPIN_REFERENCES, key=lambda x: abs(s2 - SPIN_REFERENCES[x]))


# ============================================================
# State Tracking
# ============================================================

def track_states(iters: list[dict]) -> tuple[dict, dict]:
    """
    Track states from largest to smallest CI expansion using Hungarian matching.
    
    Args:
        iters: List of iteration dictionaries (sorted by n_det)
        
    Returns:
        trajectories: Dictionary mapping state index to list of (rPT2, energy) points
        final_s2: Dictionary mapping state index to S2 value at smallest iteration
    """
    # Start with the largest wave function (most accurate)
    ref = [fingerprint(s) for s in iters[-1]["states"]]
    n_states = len(ref)
    
    trajectories = {i: [] for i in range(n_states)}
    final_s2 = {}
    
    # Start from the largest and work backwards to the smallest
    prev = ref
    
    # Iterate from second-largest to smallest
    for idx in range(len(iters) - 2, -1, -1):
        curr = [fingerprint(s) for s in iters[idx]["states"]]
        
        # Cost matrix: compare previous (larger) to current (smaller)
        C = np.zeros((len(prev), len(curr)))
        for i in range(len(prev)):
            for j in range(len(curr)):
                C[i, j] = distance(prev[i], curr[j])
        
        # Hungarian matching to find best state correspondences
        row_indices, col_indices = linear_sum_assignment(C)
        
        # Store the trajectory for this state
        for i, j in zip(row_indices, col_indices):
            rpt2_val = iters[idx]["states"][j].get("rpt2", 0.0)
            energy_val = iters[idx]["states"][j].get("energy", np.nan)
            trajectories[i].append((rpt2_val, energy_val))
            
            # If this is the smallest iteration, store its S2 value
            if idx == 0:
                final_s2[i] = curr[j]["s2"]
        
        # Update reference to the current (smaller) wave function
        prev = [curr[j] for _, j in zip(row_indices, col_indices)]
    
    # Each trajectory now has points from largest to smallest rPT2
    # Reverse to get smallest to largest for extrapolation
    for i in range(n_states):
        trajectories[i].reverse()
    
    return trajectories, final_s2


# ============================================================
# Extrapolation
# ============================================================

def extrapolate(points: list[tuple]) -> tuple | None:
    """
    Extrapolate energy to zero rPT2 using a weighted linear fit.
    
    Args:
        points: List of (rPT2, energy) tuples
        
    Returns:
        (E0, error) tuple or None if insufficient points
    """
    # Sort by absolute rPT2 value (smallest first)
    points_sorted = sorted(points, key=lambda x: abs(x[0]))
    if len(points_sorted) < 3:
        return None
    
    best_fit = None
    # Try fits with 3 to 6 points (or up to the available points)
    max_points = min(6, len(points_sorted))
    
    for k in range(3, max_points + 1):
        sub_points = np.array(points_sorted[:k])
        x = sub_points[:, 0].astype(float)
        y = sub_points[:, 1].astype(float)
        
        # Weights: 1 / rPT2² (gives more weight to smaller rPT2 values)
        weights = 1.0 / (x ** 2)
        # Weighted linear fit
        p = np.polyfit(x, y, 1, w=weights)
        y_fit = np.polyval(p, x)
        # Error estimate: standard deviation of residuals
        fit_error = np.std(y - y_fit)
        e0 = p[1]  # Intercept (E at rPT2 = 0)
        
        if best_fit is None or fit_error < best_fit[1]:
            best_fit = (e0, fit_error)
    
    return best_fit


# ============================================================
# Pretty Print Functions
# ============================================================

def print_summary(filename: str, iters: list[dict], n_states: int):
    """
    Print a summary of the extrapolation setup.
    
    Args:
        filename: Input file name
        iters: List of iterations
        n_states: Number of states tracked
    """
    ndets = [it["n_det"] for it in iters]
    
    print("\n" + "═" * 78)
    print(" Selected CI Extrapolation Summary ".center(78))
    print("═" * 78)
    
    print(f" File              : {filename}")
    print(f" Iterations        : {len(iters)}")
    print(f" States tracked    : {n_states}")
    print(f" Determinants      : {min(ndets):>10d} → {max(ndets):>10d}")
    print(f" Extrapolation     : linear fit of E_var vs rPT2")
    print(f" Fit points        : 3 to 6 best (minimum residual)")
    print(f" Weights           : 1 / rPT2²")
    print(f" State tracking    : Hungarian matching (largest → smallest)")
    print(f" Matching metric   : E + S² + variance + ex_energy fingerprint")
    
    print("═" * 78)


def print_table(energies: list[float], errors: list[float], s2vals: dict):
    """
    Print a formatted table of extrapolated energies and excitations.
    
    Args:
        energies: List of extrapolated total energies
        errors: List of extrapolation errors
        s2vals: Dictionary mapping state index to S2 value
    """
    E0 = energies[0]
    dE0 = errors[0]
    
    width = 98
    
    print("\n" + "═" * width)
    print(" Selected CI Extrapolated Excitation Energies ".center(width))
    print("═" * width)
    
    header = (
        f"│ {'#':>3} │ {'Spin':^9} │ {'<S²>':^10} │ "
        f"{'E_tot (Ha)':^18} │ {'σ(E) (Ha)':^12} │ "
        f"{'ΔE (eV)':^12} │ {'σ(ΔE) (eV)':^12} │"
    )
    
    print(header)
    print("├" + "─" * (width - 2) + "┤")
    
    for i in range(len(energies)):
        if np.isnan(energies[i]):
            continue
        
        s2 = s2vals.get(i, np.nan)
        spin = spin_label(s2)
        
        if i == 0:
            exc = 0.0
            exc_err = 0.0
        else:
            exc = (energies[i] - E0) * HA_TO_EV
            exc_err = np.sqrt(errors[i]**2 + dE0**2) * HA_TO_EV
        
        print(
            f"│ {i:3d} │ "
            f"{spin:^9} │ "
            f"{s2:10.4f} │ "
            f"{energies[i]:18.10f} │ "
            f"{errors[i]:12.3e} │ "
            f"{exc:12.3f} │ "
            f"{exc_err:12.3f} │"
        )
    
    print("═" * width)


# ============================================================
# Main Function
# ============================================================

def compute(filename: str):
    """
    Main computation workflow.
    
    Args:
        filename: Path to JSON data file
    """
    # Load and parse data
    iters = load_data(filename)
    
    # Track states across iterations
    trajectories, s2_values = track_states(iters)
    
    # Extrapolate each state
    energies = []
    errors = []
    for i in range(len(trajectories)):
        result = extrapolate(trajectories[i])
        if result is None:
            energies.append(np.nan)
            errors.append(np.nan)
        else:
            e, err = result
            energies.append(e)
            errors.append(err)
    
    # Print results
    print_summary(filename, iters, len(trajectories))
    print_table(energies, errors, s2_values)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extrapolate SCI energies to the complete basis set limit.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python3 sci_extrapolator.py HF/aug-cc-pvdz/MoreStates/HF/json/00003.json
        """
    )
    parser.add_argument("filename", help="Path to the JSON data file.")
    args = parser.parse_args()
    compute(args.filename)
