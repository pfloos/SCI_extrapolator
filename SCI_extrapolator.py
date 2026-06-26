import sys
import json
import numpy as np
from scipy.optimize import linear_sum_assignment

HaToEv = 27.211386245988


# ============================================================
# Extract all CI iterations recursively
# ============================================================
def extract_iterations(obj):
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


def load_data(filename):
    with open(filename) as f:
        data = json.load(f)

    iters = extract_iterations(data)
    iters.sort(key=lambda x: x["n_det"])

    return iters


# ============================================================
# Fingerprint for state matching
# ============================================================
def fingerprint(st):
    return {
        "energy": st["energy"],
        "s2": st.get("s2", 0.0),
        "variance": st.get("variance", 0.0),
        "ex": np.array(st.get("ex_energy", []), dtype=float),
    }


def distance(a, b):

    dE = abs(a["energy"] - b["energy"])
    dS = abs(a["s2"] - b["s2"])
    dV = abs(a["variance"] - b["variance"])

    if len(a["ex"]) and len(b["ex"]):
        n = min(len(a["ex"]), len(b["ex"]))
        dX = np.linalg.norm(a["ex"][:n] - b["ex"][:n])
    else:
        dX = 0.0

    return dE + 0.1*dS + 0.1*dV + 0.5*dX


# ============================================================
# Spin label
# ============================================================
def spin_label(s2):

    refs = {
        "Singlet": 0.0,
        "Doublet": 0.75,
        "Triplet": 2.0,
        "Quartet": 3.75,
        "Quintet": 6.0
    }

    return min(refs, key=lambda x: abs(s2 - refs[x]))


# ============================================================
# State tracking
# ============================================================
def track_states(iters):

    ref = [fingerprint(s) for s in iters[-1]["states"]]
    n_states = len(ref)

    trajectories = {i: [] for i in range(n_states)}
    final_s2 = {}

    prev = ref

    for idx, it in enumerate(reversed(iters)):

        curr = [fingerprint(s) for s in it["states"]]

        C = np.zeros((len(prev), len(curr)))

        for i in range(len(prev)):
            for j in range(len(curr)):
                C[i, j] = distance(prev[i], curr[j])

        row, col = linear_sum_assignment(C)

        new_prev = []

        for i, j in zip(row, col):

            trajectories[i].append(
                (it["states"][j]["rpt2"], it["states"][j]["energy"])
            )

            new_prev.append(curr[j])

            if idx == 0:
                final_s2[i] = curr[j]["s2"]

        prev = new_prev

    return trajectories, final_s2


# ============================================================
# Extrapolation
# ============================================================
def extrapolate(points):

    points = sorted(points, key=lambda x: abs(x[0]))

    if len(points) < 3:
        return None

    best = None

    for k in range(3, min(6, len(points)) + 1):

        sub = np.array(points[:k])

        x = sub[:, 0]
        y = sub[:, 1]

        w = 1.0 / (np.abs(x) ** 2)

        p = np.polyfit(x, y, 1, w=np.sqrt(w))

        yfit = np.polyval(p, x)

        err = np.std(y - yfit)

        E0 = p[1]

        if best is None or err < best[1]:
            best = (E0, err)

    return best


# ============================================================
# Pretty print table
# ============================================================
def print_table(energies, errors, s2vals):

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
    print("├" + "─"*(width-2) + "┤")

    for i in range(len(energies)):

        if np.isnan(energies[i]):
            continue

        s2 = s2vals.get(i, np.nan)
        spin = spin_label(s2)

        if i == 0:
            exc = 0.0
            exc_err = 0.0
        else:
            exc = (energies[i] - E0) * HaToEv
            exc_err = np.sqrt(errors[i]**2 + dE0**2) * HaToEv

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

def print_summary(filename, iters, n_states):

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
    print(f" State tracking    : Hungarian matching")
    print(f" Matching metric   : E + S² + variance + ex_energy fingerprint")

    print("═" * 78)

# ============================================================
# Main
# ============================================================
def compute(filename):

    iters = load_data(filename)

    traj, s2vals = track_states(iters)

    energies = []
    errors = []

    for i in range(len(traj)):

        res = extrapolate(traj[i])

        if res is None:
            energies.append(np.nan)
            errors.append(np.nan)
        else:
            E, err = res
            energies.append(E)
            errors.append(err)

    print_summary(filename, iters, len(traj))

    print_table(energies, errors, s2vals)


if __name__ == "__main__":
    compute(sys.argv[1])

