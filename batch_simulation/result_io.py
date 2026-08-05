"""Result I/O utilities for RSS-MPC pipeline.

Handles parsing and displaying MATLAB summary results.
"""

import json
from pathlib import Path
from typing import Any


def parse_summary(raw_json: str) -> dict[str, Any]:
    """Parse a JSON string returned by MATLAB's run_one_case.

    Args:
        raw_json: JSON string from MATLAB.

    Returns:
        Parsed summary dict.

    Raises:
        json.JSONDecodeError: If the string is not valid JSON.
    """
    return json.loads(str(raw_json))


def load_summary(output_dir: str | Path) -> dict[str, Any]:
    """Load summary.json from an output directory.

    Args:
        output_dir: Directory containing summary.json.

    Returns:
        Parsed summary dict.
    """
    summary_path = Path(output_dir) / "summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"summary.json not found in: {output_dir}")
    with open(summary_path, "r", encoding="utf-8") as f:
        return json.load(f)


def print_summary(summary: dict[str, Any]) -> None:
    """Print a human-readable summary to stdout.

    Args:
        summary: Parsed summary dict from MATLAB.
    """
    print("=" * 60)
    print(f"  Case ID:        {summary.get('case_id', 'N/A')}")
    print(f"  Success:        {summary.get('success', False)}")
    print(f"  Completed Steps:{summary.get('completed_steps', 0)}/{summary.get('num_steps', 0)}")
    print(f"  Position RMSE:  {summary.get('position_rmse', 0):.10f}")
    print(f"  Trajectory Cost:{summary.get('trajectory_cost', 0):.10f}")
    print(f"  Max Wheel Speed:{summary.get('max_wheel_speed', 0):.6f} m/s")
    print(f"  Max Steer Rate: {summary.get('max_steering_rate', 0):.6f} rad/s")
    print(f"  Max Constraint: {summary.get('max_constraint_violation', 0):.6e}")
    print(f"  Mean Solve Time:{summary.get('mean_solver_time', 0):.6f} s")
    print(f"  Total Solve:    {summary.get('total_solver_time', 0):.6f} s")
    print(f"  Wall Time:      {summary.get('wall_time', 0):.2f} s")
    print(f"  Solver:         {summary.get('solver', 'N/A')}")
    print(f"  CVX Status:     {summary.get('cvx_status_summary', 'N/A')}")
    print(f"  Result File:    {summary.get('result_file', 'N/A')}")
    print(f"  Output Dir:     {summary.get('output_dir', 'N/A')}")

    if not summary.get("success", False) and summary.get("failure_reason"):
        print(f"  Failure Reason: {summary['failure_reason'][:200]}...")
    print("=" * 60)
