#!/usr/bin/env python3
"""Extrapolate selected-CI energies to the FCI limit of a fixed orbital basis.

The reported ``u`` values are heuristic extrapolation-stability indicators,
not calibrated statistical standard deviations or confidence intervals.
"""

import argparse
import json
import math
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

HA_TO_EV = 27.211386245988
SPIN_REFERENCES = {
    "Singlet": 0.0,
    "Doublet": 0.75,
    "Triplet": 2.0,
    "Quartet": 3.75,
    "Quintet": 6.0,
}
COLORS = {
    "HEADER": "\033[95m", "BLUE": "\033[94m", "CYAN": "\033[96m",
    "GREEN": "\033[92m", "YELLOW": "\033[93m", "RED": "\033[91m",
    "BOLD": "\033[1m", "END": "\033[0m", "DIM": "\033[2m",
}
USE_COLORS = sys.stdout.isatty()


def colorize(text: str, color: str) -> str:
    if USE_COLORS and color in COLORS:
        return f"{COLORS[color]}{text}{COLORS['END']}"
    return text


def _close_truncated_json(content: str) -> str:
    """Repair only an incomplete JSON suffix, without changing completed data.

    Brackets inside strings are ignored. Malformed inner structure, dangling
    strings and mismatched closing brackets are deliberately not repaired.
    """
    stack = []
    in_string = False
    escape = False
    for char in content:
        if in_string:
            if escape:
                escape = False
            elif char == "\\":
                escape = True
            elif char == '"':
                in_string = False
        elif char == '"':
            in_string = True
        elif char in "{[":
            stack.append(char)
        elif char in "}]":
            if not stack or stack.pop() != {"}": "{", "]": "["}[char]:
                raise ValueError("Mismatched JSON closing bracket")
    if in_string:
        raise ValueError("Truncated JSON ends inside a string; refusing repair")
    if not stack:
        raise ValueError("Invalid JSON with no unclosed containers; refusing repair")

    # A trailing comma after the last *complete* item is safe to remove.
    text = content.rstrip()
    if text.endswith(","):
        text = text[:-1].rstrip()
    return text + "".join("}" if opening == "{" else "]" for opening in reversed(stack))


def fix_json_file(filename: str, verbose: bool = True) -> str:
    """Return a valid JSON path, optionally writing a conservatively repaired copy.

    The original input is never modified. An unrecoverable file raises ValueError
    rather than silently returning the broken original.
    """
    with open(filename, "r", encoding="utf-8-sig") as handle:
        content = handle.read()
    try:
        json.loads(content)
        if verbose:
            print(f"JSON validated: {filename}")
        return filename
    except json.JSONDecodeError as original_error:
        try:
            fixed = _close_truncated_json(content)
            json.loads(fixed)
        except (json.JSONDecodeError, ValueError) as error:
            raise ValueError(
                f"Invalid JSON ({original_error}). Automatic suffix repair failed: {error}"
            ) from original_error

    target = filename + ".fixed"
    with open(target, "w", encoding="utf-8") as handle:
        handle.write(fixed)
    if verbose:
        print(f"WARNING: recovered incomplete JSON suffix; check {target}")
    return target


def extract_iterations(obj: dict | list) -> list[dict]:
    """Recursively collect records containing both 'states' and 'n_det'."""
    if isinstance(obj, dict):
        result = [obj] if "states" in obj and "n_det" in obj else []
        for value in obj.values():
            result.extend(extract_iterations(value))
        return result
    if isinstance(obj, list):
        return [entry for value in obj for entry in extract_iterations(value)]
    return []


def _validate_iterations(iterations: list[dict]) -> None:
    if not iterations:
        raise ValueError("No SCI iterations with 'n_det' and 'states' found")
    for step, record in enumerate(iterations):
        ndet = record["n_det"]
        if isinstance(ndet, bool) or not isinstance(ndet, int) or ndet <= 0:
            raise ValueError(f"Iteration {step}: n_det must be a positive integer")
        states = record["states"]
        if not isinstance(states, list) or not states:
            raise ValueError(f"Iteration {step}: 'states' must be a nonempty list")
        for root, state in enumerate(states):
            if not isinstance(state, dict):
                raise ValueError(f"Iteration {step}, state {root}: expected an object")
            for key in ("energy", "rpt2", "s2", "variance"):
                try:
                    number = float(state[key])
                except (KeyError, TypeError, ValueError) as error:
                    raise ValueError(
                        f"Iteration {step}, state {root}: missing/invalid {key}"
                    ) from error
                if not math.isfinite(number):
                    raise ValueError(f"Iteration {step}, state {root}: nonfinite {key}")
            excitation_fingerprint = state.get("ex_energy", [])
            if not isinstance(excitation_fingerprint, list) or not all(
                isinstance(value, (int, float)) and math.isfinite(value)
                for value in excitation_fingerprint
            ):
                raise ValueError(f"Iteration {step}, state {root}: invalid ex_energy")


def load_data(filename: str) -> list[dict]:
    """Load validated SCI records sorted by number of determinants."""
    try:
        path = fix_json_file(filename)
        with open(path, "r", encoding="utf-8-sig") as handle:
            data = json.load(handle)
        iterations = extract_iterations(data)
        _validate_iterations(iterations)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise ValueError(f"Cannot load SCI data from {filename}: {error}") from error
    iterations.sort(key=lambda item: item["n_det"])
    return iterations


def fingerprint(state: dict) -> dict:
    return {
        "energy": float(state["energy"]),
        "s2": float(state["s2"]),
        "variance": float(state["variance"]),
        "ex": np.asarray(state.get("ex_energy", []), dtype=float),
    }


def distance(a: dict, b: dict) -> float:
    """Original heuristic matching metric, kept for compatibility."""
    energy_diff = abs(a["energy"] - b["energy"])
    s2_diff = abs(a["s2"] - b["s2"])
    variance_diff = abs(a["variance"] - b["variance"])
    size = min(len(a["ex"]), len(b["ex"]))
    ex_diff = np.linalg.norm(a["ex"][:size] - b["ex"][:size]) if size else 0.0
    return float(energy_diff + 0.1 * s2_diff + 0.1 * variance_diff + 0.5 * ex_diff)


def spin_label(s2: float) -> str:
    if not math.isfinite(s2):
        return "Unknown"
    return min(SPIN_REFERENCES, key=lambda spin: abs(s2 - SPIN_REFERENCES[spin]))


def track_states(iters: list[dict]) -> tuple[dict, dict]:
    """Match roots backwards, preserving reference IDs when root counts change.

    Each trajectory includes the largest SCI wave function. Unmatched roots
    cease at the iteration where they disappear, rather than being renumbered.
    """
    final_states = iters[-1]["states"]
    trajectories = {
        root: [(float(state["rpt2"]), float(state["energy"]))]
        for root, state in enumerate(final_states)
    }
    spin_states = {
        root: float(state["s2"]) for root, state in enumerate(final_states)
    }
    active = [(root, fingerprint(state)) for root, state in enumerate(final_states)]

    for iteration in reversed(iters[:-1]):
        current_states = iteration["states"]
        current = [fingerprint(state) for state in current_states]
        if not active or not current:
            break
        costs = np.array(
            [[distance(previous, state) for state in current] for _, previous in active],
            dtype=float,
        )
        rows, columns = linear_sum_assignment(costs)
        matched = []
        for row, column in zip(rows, columns):
            original_id = active[row][0]
            state = current_states[column]
            trajectories[original_id].append((float(state["rpt2"]), float(state["energy"])))
            matched.append((original_id, current[column]))
        active = matched

    for points in trajectories.values():
        points.reverse()
    return trajectories, spin_states


def _fit_window(points: np.ndarray) -> dict | None:
    """WLS with objective sum((residual / abs(rPT2))**2).

    The weight is a heuristic emphasis on the converged points, not a measured
    error variance. The intercept covariance is therefore only conditional.
    """
    x = points[:, 0]
    y = points[:, 1]
    if len(x) < 3 or np.any(x == 0):
        return None
    scaled_x = x / np.max(np.abs(x))
    design = np.column_stack((np.ones(len(x)), scaled_x))
    w = np.min(np.abs(x)) / np.abs(x)  # same optimum as w=1/abs(x), stable scale
    weighted_design = design * w[:, None]
    weighted_y = y * w
    coeff, _, rank, _ = np.linalg.lstsq(weighted_design, weighted_y, rcond=None)
    if rank != 2 or np.linalg.cond(weighted_design) > 1e12:
        return None
    residuals = y - design @ coeff
    dof = len(x) - 2
    weighted_residual_variance = float(np.sum((w * residuals) ** 2) / dof)
    gram_inverse = np.linalg.inv(weighted_design.T @ weighted_design)
    intercept_variance = max(0.0, weighted_residual_variance * gram_inverse[0, 0])
    return {
        "energy": float(coeff[0]),
        "fit_uncertainty": float(math.sqrt(intercept_variance)),
        "rms_residual": float(np.sqrt(np.mean(residuals ** 2))),
    }


def extrapolate(points: list[tuple], details: bool = False) -> tuple | dict | None:
    """Weighted linear extrapolation using a fixed window and window sensitivity.

    Primary result: 5 points with the smallest absolute rPT2 (or 3/4 if fewer).
    Diagnostic fits: 3 through 6 points. The reported uncertainty indicator is
    max(primary fit's residual-based intercept error, maximum difference between
    the primary intercept and any diagnostic-window intercept). It is *not* a
    calibrated 1-sigma statistical uncertainty.

    By default returns (energy, indicator), preserving the prior public API.
    Set details=True to also obtain all fit diagnostics.
    """
    finite_points = []
    for x, y in points:
        if not (math.isfinite(x) and math.isfinite(y)):
            raise ValueError("Nonfinite energy or rPT2 in extrapolation input")
        finite_points.append((float(x), float(y)))
    if not finite_points:
        return None
    exact = [y for x, y in finite_points if x == 0.0]
    if exact:
        if np.ptp(exact) > 1e-9:
            raise ValueError("Inconsistent energies at zero rPT2")
        result = {
            "energy": float(np.mean(exact)), "uncertainty": 0.0,
            "fit_uncertainty": 0.0, "window_sensitivity": 0.0,
            "n_points": len(finite_points), "fit_points": 1,
            "windows": {}, "direct_zero": True,
        }
        return result if details else (result["energy"], result["uncertainty"])
    if len(finite_points) < 3:
        return None

    sorted_points = np.asarray(
        sorted(finite_points, key=lambda pair: abs(pair[0])), dtype=float
    )
    windows = {}
    for k in range(3, min(6, len(sorted_points)) + 1):
        result = _fit_window(sorted_points[:k])
        if result is not None:
            windows[k] = result
    if not windows:
        return None
    primary_k = min(windows, key=lambda k: (abs(k - 5), -k))
    primary = windows[primary_k]
    window_sensitivity = max(
        abs(candidate["energy"] - primary["energy"]) for candidate in windows.values()
    )
    uncertainty = max(primary["fit_uncertainty"], window_sensitivity)
    result = {
        "energy": primary["energy"], "uncertainty": uncertainty,
        "fit_uncertainty": primary["fit_uncertainty"],
        "window_sensitivity": window_sensitivity,
        "n_points": len(finite_points), "fit_points": primary_k,
        "windows": windows, "direct_zero": False,
    }
    return result if details else (result["energy"], result["uncertainty"])


def print_summary(filename: str, iters: list[dict], n_states: int) -> None:
    ndets = [record["n_det"] for record in iters]
    print("\n" + colorize("═" * 78, "CYAN"))
    print(colorize(" Selected CI Extrapolation Summary ".center(78), "BOLD"))
    print(colorize("═" * 78, "CYAN"))
    for label, value in (
        ("File", filename),
        ("Iterations", str(len(iters))),
        ("States tracked", str(n_states)),
        ("Determinants", f"{min(ndets):d} → {max(ndets):d}"),
        ("Extrapolation", "linear fit of E_var vs rPT2"),
        ("Primary fit", "5 most converged points (or fewer if unavailable)"),
        ("Window sensitivity", "3–6 most converged points"),
        ("Weights", "objective: sum[(residual / abs(rPT2))²]"),
        ("State tracking", "Hungarian matching (largest → smallest)"),
        ("Matching metric", "E + S² + variance + ex_energy fingerprint"),
    ):
        print(f" {label:<29} : {value}")
    print(colorize("═" * 78, "CYAN"))


def _fmt(value: float, fmt: str, width: int) -> str:
    return format(value, fmt) if math.isfinite(value) else "N/A".rjust(width)


def print_table(energies: list[float], errors: list[float], s2vals: dict) -> None:
    e0 = energies[0]
    u0 = errors[0]
    width = 102
    print("\n" + colorize("═" * width, "CYAN"))
    print(colorize(" Extrapolated Excitation Energies ".center(width), "BOLD"))
    print(colorize("═" * width, "CYAN"))
    print(
        f"{'State':>6} │ {'Spin':^10} │ {'<S²>':>9} │ {'E_tot (Ha)':>16} │ "
        f"{'u(E) (Ha)':>11} │ {'ΔE (eV)':>12} │ {'u(ΔE) (eV)':>12}"
    )
    print(colorize("─" * width, "DIM"))
    for root, energy in enumerate(energies):
        spin_value = s2vals.get(root, np.nan)
        spin = spin_label(spin_value)
        exc = (energy - e0) * HA_TO_EV if math.isfinite(energy) and math.isfinite(e0) else np.nan
        if root == 0:
            exc_err = 0.0 if math.isfinite(e0) else np.nan
        else:
            exc_err = math.hypot(errors[root], u0) * HA_TO_EV
        print(
            f"{root:6d} │ {spin:^10} │ {_fmt(spin_value, '9.4f', 9)} │ "
            f"{_fmt(energy, '16.8f', 16)} │ {_fmt(errors[root], '11.3e', 11)} │ "
            f"{_fmt(exc, '12.3f', 12)} │ {_fmt(exc_err, '12.3f', 12)}"
        )
    print(colorize("═" * width, "CYAN"))
    print("u: heuristic fit/window-stability indicator (not a calibrated 1σ uncertainty).")
    print("Excitation-energy indicators assume independent state errors; correlations are not estimated.")


def print_compact_table(energies: list[float], errors: list[float], s2vals: dict) -> None:
    e0 = energies[0]
    u0 = errors[0]
    width = 75
    print("\n" + colorize("═" * width, "CYAN"))
    print(colorize(" Extrapolated Excitation Energies (Compact) ".center(width), "BOLD"))
    print(colorize("═" * width, "CYAN"))
    print(f"{'State':>6} {'Spin':^10} {'<S²>':^9} {'ΔE (eV)':>12} {'u(ΔE)':>11} {'E_tot (Ha)':>16}")
    print("─" * width)
    for root, energy in enumerate(energies):
        spin_value = s2vals.get(root, np.nan)
        exc = (energy - e0) * HA_TO_EV if math.isfinite(energy) and math.isfinite(e0) else np.nan
        exc_err = 0.0 if root == 0 and math.isfinite(e0) else math.hypot(errors[root], u0) * HA_TO_EV
        print(
            f"{root:6d} {spin_label(spin_value):^10} {_fmt(spin_value, '9.4f', 9)} "
            f"{_fmt(exc, '12.4f', 12)} {_fmt(exc_err, '11.4f', 11)} "
            f"{_fmt(energy, '16.8f', 16)}"
        )
    print(colorize("═" * width, "CYAN"))
    print("u: heuristic extrapolation-stability indicator, not a calibrated 1σ error bar.")


def print_fit_diagnostics(fits: list[dict | None]) -> None:
    """Display individual 3-, 4-, 5-, and 6-point extrapolated intercepts."""
    print("\nFit-window diagnostics (energies and indicators in Hartree):")
    print(
        f"{'State':>5} {'N':>3} {'Fit':>3} {'u_fit':>10} {'u_window':>10} "
        f"{'b(3)':>15} {'b(4)':>15} {'b(5)':>15} {'b(6)':>15}"
    )
    for root, fit in enumerate(fits):
        if fit is None:
            print(f"{root:5d}  no stable fit (fewer than three usable points or singular windows)")
            continue
        intercepts = [
            f"{fit['windows'][size]['energy']:15.8f}"
            if size in fit["windows"] else f"{'N/A':>15}"
            for size in (3, 4, 5, 6)
        ]
        print(
            f"{root:5d} {fit['n_points']:3d} {fit['fit_points']:3d} "
            f"{fit['fit_uncertainty']:10.3e} {fit['window_sensitivity']:10.3e} "
            + " ".join(intercepts)
        )


def compute(filename: str, compact: bool = False, diagnostics: bool = False):
    iterations = load_data(filename)
    trajectories, spins = track_states(iterations)
    fits = [extrapolate(trajectories[root], details=True) for root in range(len(trajectories))]
    energies = [fit["energy"] if fit is not None else float("nan") for fit in fits]
    errors = [fit["uncertainty"] if fit is not None else float("nan") for fit in fits]
    print_summary(filename, iterations, len(trajectories))
    if compact:
        print_compact_table(energies, errors, spins)
    else:
        print_table(energies, errors, spins)
    if diagnostics or not compact:
        print_fit_diagnostics(fits)
    return energies, errors, fits


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Extrapolate SCI energies to the FCI limit in a fixed one-electron basis.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python3 SCI_extrapolator.py example/00003.json\n"
            "  python3 SCI_extrapolator.py --compact example/00003.json"
        ),
    )
    parser.add_argument("filename", help="Path to a Quantum Package SCI JSON file")
    parser.add_argument("--compact", "-c", action="store_true", help="Compact paper-style table")
    parser.add_argument("--diagnostics", action="store_true", help="Show individual fit diagnostics")
    parser.add_argument("--no-colors", action="store_true", help="Disable terminal colors")
    args = parser.parse_args(argv)
    global USE_COLORS
    if args.no_colors:
        USE_COLORS = False
    try:
        compute(args.filename, compact=args.compact, diagnostics=args.diagnostics)
    except (ValueError, np.linalg.LinAlgError) as error:
        print(f"SCI extrapolation error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
