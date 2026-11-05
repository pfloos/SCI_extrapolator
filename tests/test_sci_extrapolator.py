"""Regression tests for SCI extrapolation and input validation (stdlib unittest)."""

import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import SCI_extrapolator as sci


def state(energy, rpt2, s2=0.0, variance=0.1, ex=None):
    return {
        "energy": energy, "rpt2": rpt2, "s2": s2,
        "variance": variance, "ex_energy": [] if ex is None else ex,
    }


class Tests(unittest.TestCase):
    def test_tracking_includes_last_iteration(self):
        iters = [
            {"n_det": 100, "states": [state(-10, -0.3), state(-9.0, -0.25)]},
            {"n_det": 200, "states": [state(-10.1, -0.15), state(-9.1, -0.1)]},
            {"n_det": 400, "states": [state(-10.2, -0.06), state(-9.2, -0.05)]},
        ]
        trajectories, _ = sci.track_states(iters)
        self.assertEqual(trajectories[0][-1], (-0.06, -10.2))
        self.assertEqual(trajectories[1][-1], (-0.05, -9.2))
        self.assertEqual([len(points) for points in trajectories.values()], [3, 3])

    def test_tracking_preserves_ids_with_varying_root_counts(self):
        # At the smallest iteration, the middle root (original ID 1) disappears.
        # Root ID 2 must not be accidentally appended to ID 1's trajectory.
        iters = [
            {"n_det": 100, "states": [state(-10.0, -0.3), state(-8.0, -0.3)]},
            {"n_det": 200, "states": [state(-10.1, -0.2), state(-9.1, -0.2), state(-8.1, -0.2)]},
            {"n_det": 300, "states": [state(-10.2, -0.1), state(-9.2, -0.1), state(-8.2, -0.1)]},
        ]
        trajectories, _ = sci.track_states(iters)
        self.assertEqual(len(trajectories[0]), 3)
        self.assertEqual(len(trajectories[1]), 2)
        self.assertEqual(len(trajectories[2]), 3)
        self.assertEqual(trajectories[2][0][1], -8.0)

    def test_exact_linear_extrapolation(self):
        xs = [-0.15, -0.1, -0.07, -0.05, -0.03, -0.02]
        data = [(x, -100.5 - 0.93 * x) for x in xs]
        fit = sci.extrapolate(data, details=True)
        self.assertAlmostEqual(fit["energy"], -100.5, places=10)
        self.assertEqual(fit["fit_points"], 5)
        self.assertEqual(sorted(fit["windows"]), [3, 4, 5, 6])
        self.assertLess(fit["uncertainty"], 1e-10)

    def test_regression_weights_match_numpy(self):
        data = [
            (-0.1, -30.0), (-0.08, -30.008), (-0.06, -30.02),
            (-0.04, -30.012), (-0.025, -30.026),
        ]
        # Five points: primary and external polyfit should agree exactly.
        fit = sci.extrapolate(data, details=True)
        x = np.array([t[0] for t in data])
        y = np.array([t[1] for t in data])
        expected = np.polyfit(x, y, 1, w=1.0 / np.abs(x))[1]
        self.assertAlmostEqual(fit["energy"], expected, places=10)

    def test_window_dependence_is_visible_and_not_three_point_selection(self):
        xs = [-0.12, -0.10, -0.08, -0.06, -0.04, -0.02]
        data = [(x, -50.0 - 1.0 * x + 0.6 * x * x) for x in xs]
        fit = sci.extrapolate(data, details=True)
        self.assertEqual(fit["fit_points"], 5)
        self.assertGreater(fit["window_sensitivity"], 0)
        self.assertGreaterEqual(fit["uncertainty"], fit["window_sensitivity"])

    def test_hf_example_trajectories(self):
        # These are the matching root trajectories extracted from the published
        # example/00003.json, including the final 299786-determinant iteration.
        fixture = Path(__file__).parent / "fixtures" / "hf_trajectories.json"
        trajectories = json.loads(fixture.read_text())["trajectories"]
        fits = [sci.extrapolate(t, details=True) for t in trajectories]
        expected_energies = [
            -100.34935819, -99.97919020, -99.97915906, -99.96500885,
            -99.96507457, -99.85388342, -99.84827297, -99.84835780,
            -99.84342696, -99.84340936, -99.83248652, -99.82047402,
        ]
        self.assertEqual(len(fits), 12)
        for fit, expected in zip(fits, expected_energies):
            self.assertAlmostEqual(fit["energy"], expected, places=7)
            self.assertEqual(fit["fit_points"], 5)
        self.assertEqual(len(trajectories[0]), 9)
        self.assertAlmostEqual(trajectories[0][-1][0], -0.009520987857783)

    def test_end_to_end_from_json_to_energies(self):
        fixture = Path(__file__).parent / "fixtures" / "hf_trajectories.json"
        data = json.loads(fixture.read_text())
        iterations = []
        for step in range(9):
            states = [
                state(
                    root_traj[step][1], root_traj[step][0],
                    s2=data["referenceSpin"][root], variance=100.0 * root,
                )
                for root, root_traj in enumerate(data["trajectories"])
            ]
            iterations.append({"n_det": 10 ** (step + 1), "states": states})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "synthetic_hf.json"
            path.write_text(json.dumps({"fci": iterations}))
            with patch("builtins.print"):
                energies, errors, fits = sci.compute(str(path))
        self.assertEqual(len(energies), 12)
        self.assertAlmostEqual(energies[0], -100.34935819, places=7)
        self.assertAlmostEqual(energies[1], -99.97919020, places=7)
        self.assertTrue(np.all(np.isfinite(errors)))
        self.assertTrue(all(fit["n_points"] == 9 for fit in fits))

    def test_zero_rpt2_does_not_crash(self):
        data = [(-0.08, -8.0), (-0.04, -8.1), (0.0, -8.15)]
        self.assertEqual(sci.extrapolate(data), (-8.15, 0.0))
        with self.assertRaises(ValueError):
            sci.extrapolate(data + [(0.0, -8.5)])

    def test_identical_rpt2_is_not_fitted(self):
        self.assertIsNone(sci.extrapolate([(-0.1, -5.0)] * 5))

    def test_one_root_excitation_is_exact_zero(self):
        with patch("builtins.print") as print_mock:
            sci.print_table([-4.2], [0.01], {0: 0.0})
        output = "\n".join(" ".join(str(arg) for arg in call.args) for call in print_mock.call_args_list)
        self.assertIn("0.000", output)
        self.assertIn("not a calibrated 1σ", output)

    def test_json_recovery_adds_only_missing_closers(self):
        payload = '{"fci": [{"n_det": 100, "states": [' + json.dumps(state(-10, -0.1))
        repaired = sci._close_truncated_json(payload)
        data = json.loads(repaired)
        self.assertEqual(data["fci"][0]["n_det"], 100)

    def test_json_recovery_ignores_brackets_inside_strings(self):
        payload = '{"description": "literal {[}]", "fci": []'
        self.assertEqual(json.loads(sci._close_truncated_json(payload))["description"], "literal {[}]")
        with self.assertRaises(ValueError):
            sci._close_truncated_json('{"a": "unfinished')

    def test_json_with_bom_is_accepted(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "bom.json"
            path.write_text(json.dumps({"fci": [
                {"n_det": 10, "states": [state(-10, -0.1)]},
            ]}), encoding="utf-8-sig")
            self.assertEqual(len(sci.load_data(str(path))), 1)

    def test_file_loading_and_missing_field_validation(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "sample.json"
            source.write_text(json.dumps({"fci": [
                {"n_det": 20, "states": [state(-3, -0.1)]},
                {"n_det": 10, "states": [state(-2, -0.2)]},
            ]}))
            self.assertEqual([record["n_det"] for record in sci.load_data(str(source))], [10, 20])
            source.write_text(json.dumps({"fci": [
                {"n_det": 10, "states": [{"energy": -2, "s2": 0, "variance": 0.2}]},
            ]}))
            with self.assertRaisesRegex(ValueError, "rpt2"):
                sci.load_data(str(source))

    def test_json_repaired_file_keeps_original(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "truncated.json"
            original = '{"fci": [{"n_det": 10, "states": [' + json.dumps(state(-10, -0.1))
            source.write_text(original)
            fixed = sci.fix_json_file(str(source), verbose=False)
            self.assertEqual(source.read_text(), original)
            self.assertTrue(Path(fixed).exists())
            self.assertEqual(len(sci.load_data(str(source))), 1)


if __name__ == "__main__":
    unittest.main()
