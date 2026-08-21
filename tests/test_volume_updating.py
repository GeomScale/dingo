from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.io import loadmat, savemat


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "dingo" / "dynamic_volume.py"
if not MODULE_PATH.is_file():
    MODULE_PATH = ROOT / "dingo" / "dynamic_volume.py"
SPEC = importlib.util.spec_from_file_location("volume_updating_test", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

HalfspaceCut = MODULE.HalfspaceCut
StaticVolumeUpdater = MODULE.StaticVolumeUpdater
cuts_from_bounds = MODULE.cuts_from_bounds


class FakeUpdater(StaticVolumeUpdater):
    last_request = None

    def _run_matlab(self, request_path, output_path, timeout):
        self.last_request = loadmat(request_path, squeeze_me=True, struct_as_record=False)
        n_cuts = np.atleast_2d(self.last_request["cut_normals"]).shape[0]
        ratios = np.array([0.5, 0.25], dtype=float)[:n_cuts]
        logs = np.log(ratios)
        savemat(
            output_path,
            {
                "rho": float(np.prod(ratios)),
                "log_rho": float(np.sum(logs)),
                "cut_ratios": ratios,
                "cut_log_ratios": logs,
                "cut_n_phases": np.arange(1, n_cuts + 1),
                "cut_total_steps": 10 * np.arange(1, n_cuts + 1),
                "cut_converged": np.ones(n_cuts, dtype=np.uint8),
                "algorithm": "static_relative_volume_sequential",
                "objective_floor": np.uint8(0),
            },
        )


class VolumeUpdatingWrapperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.matlab_dir = MODULE_PATH.parent / "matlab" / "volume_updating"

    def updater(self):
        return FakeUpdater(
            matlab_executable="/bin/true",
            matlab_source_dir=self.matlab_dir,
        )

    def test_sequential_static_result_contract(self):
        updater = self.updater()
        cuts = [
            HalfspaceCut(np.array([1.0, 0.0]), 0.5, "x"),
            HalfspaceCut(np.array([0.0, 1.0]), 0.25, "y"),
        ]
        result = updater.estimate(
            np.zeros((0, 2)),
            np.zeros(2),
            np.ones(2),
            cuts,
            options={"epsilon": 0.1, "algorithm": "static"},
        )
        self.assertAlmostEqual(result.rho, 0.125)
        self.assertAlmostEqual(result.log_rho, np.log(0.125))
        np.testing.assert_allclose(result.cut_ratios, [0.5, 0.25])
        self.assertEqual(result.cut_labels, ("x", "y"))
        self.assertFalse(result.objective_floor)
        np.testing.assert_allclose(
            np.atleast_2d(updater.last_request["cut_normals"]),
            [[1.0, 0.0], [0.0, 1.0]],
        )

    def test_empty_cut_stream_is_identity(self):
        result = self.updater().estimate(
            np.zeros((0, 2)), np.zeros(2), np.ones(2), []
        )
        self.assertEqual(result.rho, 1.0)
        self.assertEqual(result.cut_ratios.size, 0)

    def test_reuse_algorithm_is_not_faked(self):
        with self.assertRaisesRegex(ValueError, "does not provide a separate reuse"):
            self.updater().estimate(
                np.zeros((0, 2)),
                np.zeros(2),
                np.ones(2),
                [HalfspaceCut(np.ones(2), 1.0)],
                options={"algorithm": "reuse"},
            )

    def test_unknown_non_source_option_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "not part of the supplied static API"):
            self.updater().estimate(
                np.zeros((0, 2)),
                np.zeros(2),
                np.ones(2),
                [HalfspaceCut(np.ones(2), 1.0)],
                options={"ratio_ess_min": 3000},
            )

    def test_invalid_static_option_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "positive integer"):
            self.updater().estimate(
                np.zeros((0, 2)),
                np.zeros(2),
                np.ones(2),
                [HalfspaceCut(np.ones(2), 1.0)],
                options={"N_utest": 12.5},
            )
        with self.assertRaisesRegex(ValueError, r"r \+ delta"):
            self.updater().estimate(
                np.zeros((0, 2)),
                np.zeros(2),
                np.ones(2),
                [HalfspaceCut(np.ones(2), 1.0)],
                options={"r": 0.8, "delta": 0.3},
            )
        with self.assertRaisesRegex(ValueError, r"r \+ delta"):
            self.updater().estimate(
                np.zeros((0, 2)),
                np.zeros(2),
                np.ones(2),
                [HalfspaceCut(np.ones(2), 1.0)],
                options={"r": 0.98},
            )

    def test_cut_dimension_is_checked(self):
        with self.assertRaisesRegex(ValueError, "expected 2"):
            self.updater().estimate(
                np.zeros((0, 2)),
                np.zeros(2),
                np.ones(2),
                [HalfspaceCut(np.ones(3), 1.0)],
            )

    def test_bound_update_compiles_to_two_sided_cuts(self):
        cuts = cuts_from_bounds(
            np.array([-2.0, 0.0]),
            np.array([2.0, 4.0]),
            np.array([-1.0, 0.0]),
            np.array([1.5, 3.0]),
        )
        self.assertEqual([cut.label for cut in cuts], ["ub[0]", "lb[0]", "ub[1]"])
        np.testing.assert_allclose(cuts[0].normal, [1.0, 0.0])
        np.testing.assert_allclose(cuts[1].normal, [-1.0, 0.0])
        self.assertEqual(cuts[1].threshold, 1.0)

    def test_bound_relaxation_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "subset"):
            cuts_from_bounds([0.0], [1.0], [-0.1], [1.0])

    def test_sampler_interface_initialization_is_present(self):
        source = (
            MODULE_PATH.parent
            / "matlab"
            / "volume_updating"
            / "relative_volume.m"
        ).read_text(encoding="utf-8")
        self.assertIn("sampler_opts.startTime = start_time;", source)
        self.assertIn("rng(sampler_opts.seed, 'simdTwister');", source)


if __name__ == "__main__":
    unittest.main()
