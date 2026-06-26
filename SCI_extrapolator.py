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

# ANSI color codes for terminal output (optional)
COLORS = {
    "HEADER": '\033[95m',
    "BLUE": '\033[94m',
    "CYAN": '\033[96m',
    "GREEN": '\033[92m',
    "YELLOW": '\033[93m',
    "RED": '\033[91m',
    "BOLD": '\033[1m',
    "UNDERLINE": '\033[4m',
    "END": '\033[0m',
    "DIM": '\033[2m'
}

# Try to determine if we're in a terminal that supports colors
try:
    USE_COLORS = sys.stdout.isatty()
except:
    USE_COLORS = False


def colorize(text: str, color: str) -> str:
    """Apply color to text if terminal supports it."""
    if USE_COLORS and color in COLORS:
        return f"{COLORS[color]}{text}{COLORS['END']}"
    return text


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
        spin_states: Dictionary mapping state index to S2 value from largest iteration
    """
    # Start with the largest wave function (most accurate)
    ref = [fingerprint(s) for s in iters[-1]["states"]]
    n_states = len(ref)
    
    trajectories = {i: [] for i in range(n_states)}
    
    # Store S2 values from the largest (most accurate) wave function
    spin_states = {}
    for i, state in enumerate(iters[-1]["states"]):
        spin_states[i] = state.get("s2", 0.0)
    
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
        
        # Update reference to the current (smaller) wave function
        prev = [curr[j] for _, j in zip(row_indices, col_indices)]
    
    # Each trajectory now has points from largest to smallest rPT2
    # Reverse to get smallest to largest for extrapolation
    for i in range(n_states):
        trajectories[i].reverse()
    
    return trajectories, spin_states


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
    
    print("\n" + colorize("═" * 78, "CYAN"))
    print(colorize(" Selected CI Extrapolation Summary ".center(78), "BOLD"))
    print(colorize("═" * 78, "CYAN"))
    
    print(f" {colorize('File', 'BOLD'):<30} : {filename}")
    print(f" {colorize('Iterations', 'BOLD'):<30} : {len(iters)}")
    print(f" {colorize('States tracked', 'BOLD'):<30} : {n_states}")
    print(f" {colorize('Determinants', 'BOLD'):<30} : {min(ndets):>10d} → {max(ndets):>10d}")
    print(f" {colorize('Extrapolation', 'BOLD'):<30} : linear fit of E_var vs rPT2")
    print(f" {colorize('Fit points', 'BOLD'):<30} : 3 to 6 best (minimum residual)")
    print(f" {colorize('Weights', 'BOLD'):<30} : 1 / rPT2²")
    print(f" {colorize('State tracking', 'BOLD'):<30} : Hungarian matching (largest → smallest)")
    print(f" {colorize('Matching metric', 'BOLD'):<30} : E + S² + variance + ex_energy fingerprint")
    print(f" {colorize('Spin assignment', 'BOLD'):<30} : from largest wave function")
    
    print(colorize("═" * 78, "CYAN"))


def print_table(energies: list[float], errors: list[float], s2vals: dict):
    """
    Print a beautifully formatted table of extrapolated energies and excitations.
    
    Args:
        energies: List of extrapolated total energies
        errors: List of extrapolation errors
        s2vals: Dictionary mapping state index to S2 value (from largest wave function)
    """
    E0 = energies[0]
    dE0 = errors[0]
    
    # Define column widths (including padding)
    col_widths = {
        'state': 6,      # width for state number
        'spin': 10,      # width for spin label
        's2': 10,        # width for S2 value
        'energy': 16,    # width for total energy
        'error': 14,     # width for energy error
        'excitation': 14, # width for excitation energy
        'exc_error': 14   # width for excitation error
    }
    
    # Calculate total width
    total_width = sum(col_widths.values()) + 3 * 6  # 6 separators (│) with spaces
    
    # Print header with decoration
    print("\n" + colorize("═" * total_width, "CYAN"))
    print(colorize(" Extrapolated Excitation Energies ".center(total_width), "BOLD"))
    print(colorize("═" * total_width, "CYAN"))
    
    # Print column headers with proper alignment
    header = (
        f"{colorize('State', 'BOLD'):>{col_widths['state']}} │ "
        f"{colorize('Spin', 'BOLD'):^{col_widths['spin']}} │ "
        f"{colorize('<S²>', 'BOLD'):^{col_widths['s2']}} │ "
        f"{colorize('E_total (Ha)', 'BOLD'):^{col_widths['energy']}} │ "
        f"{colorize('σ(E) (Ha)', 'BOLD'):^{col_widths['error']}} │ "
        f"{colorize('ΔE (eV)', 'BOLD'):^{col_widths['excitation']}} │ "
        f"{colorize('σ(ΔE) (eV)', 'BOLD'):^{col_widths['exc_error']}}"
    )
    print(header)
    
    # Print separator with proper widths
    separator = (
        "─" * col_widths['state'] + "┼" +
        "─" * col_widths['spin'] + "┼" +
        "─" * col_widths['s2'] + "┼" +
        "─" * col_widths['energy'] + "┼" +
        "─" * col_widths['error'] + "┼" +
        "─" * col_widths['excitation'] + "┼" +
        "─" * col_widths['exc_error']
    )
    print(colorize(separator, "DIM"))
    
    # Print each state
    for i in range(len(energies)):
        if np.isnan(energies[i]):
            continue
        
        s2 = s2vals.get(i, np.nan)
        spin = spin_label(s2)
        
        # State label
        state_label = f"{i:>3d}"
        
        exc = (energies[i] - E0) * HA_TO_EV
        exc_err = np.sqrt(errors[i]**2 + dE0**2) * HA_TO_EV
        
        # Format values
        energy_str = f"{energies[i]:14.8f}" if not np.isnan(energies[i]) else "N/A"
        error_str = f"{errors[i]:10.3e}" if not np.isnan(errors[i]) else "N/A"
        
        # Color only the error columns based on magnitude
        if not np.isnan(errors[i]):
            if errors[i] < 1e-4:
                error_str = colorize(error_str, "GREEN")
            elif errors[i] < 1e-3:
                error_str = colorize(error_str, "YELLOW")
            else:
                error_str = colorize(error_str, "RED")
        
        # Color excitation error
        if not np.isnan(exc_err):
            if exc_err < 0.01:
                exc_err_str = colorize(f"{exc_err:10.3f}", "GREEN")
            elif exc_err < 0.1:
                exc_err_str = colorize(f"{exc_err:10.3f}", "YELLOW")
            else:
                exc_err_str = colorize(f"{exc_err:10.3f}", "RED")
        else:
            exc_err_str = "N/A"
        
        # Format excitation
        exc_str = f"{exc:10.3f}" if i >= 0 else "N/A"
        
        # Print row with proper alignment
        row = (
            f"{state_label:>{col_widths['state']}} │ "
            f"{spin:^{col_widths['spin']}} │ "
            f"{s2:^{col_widths['s2']}.4f} │ "
            f"{energy_str:^{col_widths['energy']}} │ "
            f"{error_str:^{col_widths['error']}} │ "
            f"{exc_str:^{col_widths['excitation']}} │ "
            f"{exc_err_str:^{col_widths['exc_error']}}"
        )
        print(row)
    
    # Print footer
    print(colorize("═" * total_width, "CYAN"))
    
    # Print legend for colored entries (only errors)
    if USE_COLORS:
        print(f"\n{colorize('Error legend:', 'BOLD')} "
              f"{colorize('σ(E) < 1e-4', 'GREEN')} │ "
              f"{colorize('1e-4 < σ(E) < 1e-3', 'YELLOW')} │ "
              f"{colorize('σ(E) > 1e-3', 'RED')} │ "
              f"{colorize('σ(ΔE) < 0.01 eV', 'GREEN')} │ "
              f"{colorize('0.01 < σ(ΔE) < 0.1 eV', 'YELLOW')} │ "
              f"{colorize('σ(ΔE) > 0.1 eV', 'RED')}")


def print_compact_table(energies: list[float], errors: list[float], s2vals: dict):
    """
    Print a compact table format suitable for inclusion in papers.
    
    Args:
        energies: List of extrapolated total energies
        errors: List of extrapolation errors
        s2vals: Dictionary mapping state index to S2 value (from largest wave function)
    """
    E0 = energies[0]
    dE0 = errors[0]
    
    # Define column widths
    col_widths = {
        'state': 6,
        'spin': 10,
        's2': 8,
        'exc': 12,
        'exc_err': 10,
        'energy': 14
    }
    
    total_width = sum(col_widths.values()) + 5 * 3  # separators
    
    print("\n" + colorize("═" * total_width, "CYAN"))
    print(colorize(" Extrapolated Excitation Energies (Compact) ".center(total_width), "BOLD"))
    print(colorize("═" * total_width, "CYAN"))
    
    # Header
    header = (
        f"{colorize('State', 'BOLD'):>{col_widths['state']}} │ "
        f"{colorize('Spin', 'BOLD'):^{col_widths['spin']}} │ "
        f"{colorize('<S²>', 'BOLD'):^{col_widths['s2']}} │ "
        f"{colorize('ΔE (eV)', 'BOLD'):^{col_widths['exc']}} │ "
        f"{colorize('σ(ΔE)', 'BOLD'):^{col_widths['exc_err']}} │ "
        f"{colorize('E_tot (Ha)', 'BOLD'):^{col_widths['energy']}}"
    )
    print(header)
    
    # Separator
    separator = (
        "─" * col_widths['state'] + "┼" +
        "─" * col_widths['spin'] + "┼" +
        "─" * col_widths['s2'] + "┼" +
        "─" * col_widths['exc'] + "┼" +
        "─" * col_widths['exc_err'] + "┼" +
        "─" * col_widths['energy']
    )
    print(colorize(separator, "DIM"))
    
    # Print each state
    for i in range(len(energies)):
        if np.isnan(energies[i]):
            continue
        
        s2 = s2vals.get(i, np.nan)
        spin = spin_label(s2)
        
        exc = (energies[i] - E0) * HA_TO_EV
        exc_err = np.sqrt(errors[i]**2 + dE0**2) * HA_TO_EV
        
        energy_str = f"{energies[i]:14.8f}" if not np.isnan(energies[i]) else "N/A"
        
        # Color only errors
        error_str = f"{exc_err:10.4f}" if i > 0 else f"{0.0:10.4f}"
        if USE_COLORS:
            if exc_err < 0.01:
                error_str = colorize(error_str, "GREEN")
            elif exc_err < 0.1:
                error_str = colorize(error_str, "YELLOW")
            else:
                error_str = colorize(error_str, "RED")
        
        # Print row
        row = (
            f"{i:>{col_widths['state']}} │ "
            f"{spin:^{col_widths['spin']}} │ "
            f"{s2:^{col_widths['s2']}.4f} │ "
            f"{exc:^{col_widths['exc']}.4f} │ "
            f"{error_str:^{col_widths['exc_err']}} │ "
            f"{energy_str:^{col_widths['energy']}}"
        )
        print(row)
    
    print(colorize("═" * total_width, "CYAN"))


# ============================================================
# Main Function
# ============================================================

def compute(filename: str, compact: bool = False):
    """
    Main computation workflow.
    
    Args:
        filename: Path to JSON data file
        compact: If True, use compact table format
    """
    # Load and parse data
    iters = load_data(filename)
    
    # Track states across iterations
    trajectories, spin_states = track_states(iters)
    
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
    
    if compact:
        print_compact_table(energies, errors, spin_states)
    else:
        print_table(energies, errors, spin_states)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extrapolate SCI energies to the complete basis set limit.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python3 sci_extrapolator.py HF/aug-cc-pvdz/MoreStates/HF/json/00003.json
  python3 sci_extrapolator.py --compact HF/aug-cc-pvdz/MoreStates/HF/json/00003.json
        """
    )
    parser.add_argument("filename", help="Path to the JSON data file.")
    parser.add_argument("--compact", "-c", action="store_true", 
                       help="Use compact table format (good for papers)")
    parser.add_argument("--no-colors", action="store_true",
                       help="Disable colored output")
    args = parser.parse_args()
    
    # Override color setting if requested
    if args.no_colors:
        USE_COLORS = False
    
    compute(args.filename, compact=args.compact)
