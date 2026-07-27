import math
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(__file__))
PYTHON_DIR = os.path.join(ROOT, "python")
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

from rss_controller import build_reference_path, compute_control_step, normalize_angle, run_closed_loop


class RssControllerTests(unittest.TestCase):
    def test_build_reference_path_shape(self):
        path = build_reference_path(num_points=6)
        self.assertEqual(path.shape[0], 3)
        self.assertEqual(path.shape[1], 6)

    def test_normalize_angle_wraps(self):
        self.assertAlmostEqual(normalize_angle(3 * math.pi), -math.pi)

    def test_compute_control_step_returns_body_velocity(self):
        path = build_reference_path(num_points=8)
        state = [0.0, 0.0, 0.0]
        state_dot = [0.1, 0.0, 0.0]
        control = compute_control_step(path, 1, state_dot, state)
        self.assertEqual(len(control), 3)
        self.assertTrue(all(math.isfinite(v) for v in control))

    def test_run_closed_loop_returns_history(self):
        path = build_reference_path(num_points=8)
        history = run_closed_loop(path, num_steps=3)
        self.assertEqual(len(history["states"]), 3)
        self.assertEqual(len(history["controls"]), 3)
        self.assertTrue(all(len(state) == 3 for state in history["states"]))


if __name__ == "__main__":
    unittest.main()
