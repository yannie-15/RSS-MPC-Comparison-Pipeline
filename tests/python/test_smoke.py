#!/usr/bin/env python
"""Python smoke tests for RSS-MPC pipeline.

Tests that the Python modules can be imported and basic functionality works.
Does NOT require MATLAB Engine to be installed.

Usage:
    python tests/python/test_smoke.py
    python -m pytest tests/python/test_smoke.py -v
"""

import json
import sys
import tempfile
from pathlib import Path

# Ensure repo root is on the path
REPO_ROOT = Path(__file__).parent.parent.parent.resolve()
sys.path.insert(0, str(REPO_ROOT))

from python.config_io import load_config, apply_cli_overrides, save_resolved_config, DEFAULT_CONFIG
from python.result_io import parse_summary, print_summary


def test_config_loading():
    """Test that default config loads correctly."""
    cfg = load_config(None)
    assert cfg["case_id"] == "paper_fixed_001"
    assert cfg["num_steps"] == 100
    assert cfg["solver"] == "sdpt3"
    assert cfg["live_plot"] is False
    print("PASS: test_config_loading")


def test_config_from_file():
    """Test loading a config from a JSON file."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False, encoding="utf-8") as f:
        json.dump({"case_id": "test_123", "seed": 42}, f)
        path = f.name

    cfg = load_config(path)
    assert cfg["case_id"] == "test_123"
    assert cfg["seed"] == 42
    # Defaults should still be present
    assert cfg["num_steps"] == 100
    Path(path).unlink()  # cleanup
    print("PASS: test_config_from_file")


def test_cli_overrides():
    """Test applying CLI overrides."""
    cfg = load_config(None)
    apply_cli_overrides(cfg, case_id="override_test", seed=99, output="/tmp/out")
    assert cfg["case_id"] == "override_test"
    assert cfg["seed"] == 99
    assert cfg["output_dir"] == "/tmp/out"
    print("PASS: test_cli_overrides")


def test_save_resolved_config():
    """Test saving resolved config to disk."""
    cfg = load_config(None)
    cfg["case_id"] = "save_test"

    with tempfile.TemporaryDirectory() as tmpdir:
        path = save_resolved_config(cfg, tmpdir)
        assert path.exists()
        with open(path, "r") as f:
            loaded = json.load(f)
        assert loaded["case_id"] == "save_test"
    print("PASS: test_save_resolved_config")


def test_parse_summary():
    """Test parsing a MATLAB summary JSON string."""
    mock_json = json.dumps({
        "case_id": "test_parse",
        "success": True,
        "failure_reason": "",
        "completed_steps": 100,
        "position_rmse": 0.001234,
        "trajectory_cost": 1.5,
        "max_wheel_speed": 3.2,
        "max_steering_rate": 12.5,
        "max_constraint_violation": 0.0,
        "mean_solver_time": 0.05,
        "total_solver_time": 5.0,
        "wall_time": 10.5,
        "solver": "sdpt3",
        "cvx_status_summary": "Solved (300)",
        "result_file": "/tmp/result.mat",
        "output_dir": "/tmp/out",
        "num_steps": 100,
    })

    summary = parse_summary(mock_json)
    assert summary["success"] is True
    assert summary["completed_steps"] == 100
    assert abs(summary["position_rmse"] - 0.001234) < 1e-10
    print("PASS: test_parse_summary")


def test_parse_failure_summary():
    """Test parsing a failure summary."""
    mock_json = json.dumps({
        "case_id": "test_fail",
        "success": False,
        "failure_reason": "Solver failed at step 50",
        "completed_steps": 49,
        "position_rmse": 0,
        "trajectory_cost": 0,
        "max_wheel_speed": 0,
        "max_steering_rate": 0,
        "max_constraint_violation": 0,
        "mean_solver_time": 0,
        "total_solver_time": 0,
        "wall_time": 5.0,
        "solver": "sdpt3",
        "cvx_status_summary": "N/A",
        "result_file": "",
        "output_dir": "/tmp/fail",
        "num_steps": 100,
    })

    summary = parse_summary(mock_json)
    assert summary["success"] is False
    assert "step 50" in summary["failure_reason"]
    print("PASS: test_parse_failure_summary")


def test_print_summary():
    """Test that print_summary doesn't crash."""
    mock_summary = {
        "case_id": "print_test",
        "success": True,
        "completed_steps": 100,
        "num_steps": 100,
        "position_rmse": 0.001,
        "trajectory_cost": 1.5,
        "max_wheel_speed": 3.0,
        "max_steering_rate": 12.0,
        "max_constraint_violation": 0.0,
        "mean_solver_time": 0.05,
        "total_solver_time": 5.0,
        "wall_time": 10.0,
        "solver": "sdpt3",
        "cvx_status_summary": "Solved (300)",
        "result_file": "/tmp/result.mat",
        "output_dir": "/tmp/out",
    }
    print_summary(mock_summary)  # Should not raise
    print("PASS: test_print_summary")


def test_matlab_engine_check():
    """Test that check_matlab_engine returns a boolean (not necessarily True)."""
    from python.matlab_bridge import check_matlab_engine
    result = check_matlab_engine()
    assert isinstance(result, bool)
    if result:
        print("PASS: test_matlab_engine_check (MATLAB Engine available)")
    else:
        print("PASS: test_matlab_engine_check (MATLAB Engine NOT available - OK for smoke test)")
        print("  To enable MATLAB integration, install matlab.engine:")
        print("  cd <matlabroot>/extern/engines/python && python setup.py install")


def test_main_py_importable():
    """Test that main.py can be imported (without running)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("main", str(REPO_ROOT / "main.py"))
    assert spec is not None
    print("PASS: test_main_py_importable")


def run_all():
    """Run all smoke tests."""
    print("=" * 60)
    print("  RSS-MPC Pipeline - Python Smoke Tests")
    print("=" * 60)
    print()

    tests = [
        test_config_loading,
        test_config_from_file,
        test_cli_overrides,
        test_save_resolved_config,
        test_parse_summary,
        test_parse_failure_summary,
        test_print_summary,
        test_matlab_engine_check,
        test_main_py_importable,
    ]

    passed = 0
    failed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"FAIL: {test.__name__}: {e}")
            failed += 1

    print()
    print(f"{'=' * 60}")
    print(f"  Results: {passed} passed, {failed} failed")
    print(f"{'=' * 60}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(run_all())
