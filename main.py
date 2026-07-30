#!/usr/bin/env python
"""RSS-MPC Comparison Pipeline - Python Entry Point.

This is the recommended entry point for running RSS proposed algorithm cases.
It uses MATLAB Engine to call the MATLAB closed-loop simulation exactly once
per case (not per MPC step), preserving all CVX/SDPT3 solver behavior.

Usage:
    python main.py --config default_python_config.json
    python main.py --seed 3 --output results/seed_0003 --no-plot
    python main.py --case-id test_001 --output results/test_001

Exit codes:
    0  - Success (algorithm solved all steps)
    2  - Algorithm failure (solver failed, constraint violation, etc.)
    3  - Infrastructure error (MATLAB Engine not installed, config error, etc.)
"""

import argparse
import json
import sys
from pathlib import Path

# Ensure the python/ package is importable
REPO_ROOT = Path(__file__).parent.resolve()
sys.path.insert(0, str(REPO_ROOT))

from python.config_io import load_config, apply_cli_overrides, save_resolved_config
from python.result_io import print_summary
from python.matlab_bridge import MatlabBridge, check_matlab_engine


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="RSS-MPC Comparison Pipeline - Python entry point",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exit codes:
  0  Success
  2  Algorithm failure
  3  Infrastructure error
        """,
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to JSON config file (default: default_python_config.json)",
    )
    parser.add_argument(
        "--case-id",
        type=str,
        default=None,
        help="Case identifier (overrides config)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed (overrides config)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output directory (overrides config)",
    )
    parser.add_argument(
        "--live-plot",
        action="store_true",
        default=None,
        help="Enable live plotting in MATLAB",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        default=False,
        help="Disable live plotting (default)",
    )
    parser.add_argument(
        "--save-figures",
        action="store_true",
        default=None,
        help="Save figures to output directory",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing output directory",
    )
    return parser.parse_args()


def main() -> int:
    """Main entry point.

    Returns:
        Exit code: 0 (success), 2 (algorithm failure), 3 (infrastructure error).
    """
    args = parse_args()

    # Determine config file path
    if args.config is None:
        default_config = REPO_ROOT / "default_python_config.json"
        if default_config.exists():
            config_path = default_config
        else:
            print("ERROR: No config file specified and default_python_config.json not found.")
            return 3
    else:
        config_path = Path(args.config)
        if not config_path.is_absolute():
            config_path = REPO_ROOT / args.config

    # Load and merge config
    try:
        cfg = load_config(config_path)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        return 3
    except json.JSONDecodeError as e:
        print(f"ERROR: Invalid JSON in config file: {e}")
        return 3

    # Apply CLI overrides
    live_plot = args.live_plot if args.live_plot is not None else (not args.no_plot)
    apply_cli_overrides(
        cfg,
        case_id=args.case_id,
        seed=args.seed,
        output=args.output,
        live_plot=live_plot,
        save_figures=args.save_figures,
        force=args.force,
    )

    # Resolve output directory
    if not cfg.get("output_dir"):
        cfg["output_dir"] = str(REPO_ROOT / "results" / cfg.get("case_id", "default"))
    output_dir = Path(cfg["output_dir"])

    # Check for existing output
    if output_dir.exists() and any(output_dir.iterdir()):
        if not cfg.get("force", False):
            print(f"WARNING: Output directory already exists and is non-empty: {output_dir}")
            print("Use --force to overwrite.")
            # Continue anyway - files will be overwritten

    # Save resolved config
    try:
        resolved_path = save_resolved_config(cfg, output_dir)
        print(f"Resolved config saved to: {resolved_path}")
    except Exception as e:
        print(f"ERROR: Failed to save resolved config: {e}")
        return 3

    # Check MATLAB Engine availability
    if not check_matlab_engine():
        print("ERROR: matlab.engine is not installed.")
        print("Install MATLAB Engine for Python:")
        print(f"  1. Find your MATLAB root (run 'matlabroot' in MATLAB)")
        print(f"  2. cd <matlabroot>/extern/engines/python")
        print(f"  3. {sys.executable} setup.py install")
        return 3

    # Run the case via MATLAB Engine
    try:
        with MatlabBridge(REPO_ROOT) as bridge:
            summary = bridge.run_one_case(resolved_path)
    except ImportError as e:
        print(f"ERROR: {e}")
        return 3
    except Exception as e:
        print(f"ERROR: Infrastructure failure: {e}")
        return 3

    # Set exit code based on summary
    if summary.get("success", False):
        print(f"\nSUCCESS: Case '{summary.get('case_id', 'N/A')}' completed.")
        print(f"  Steps: {summary.get('completed_steps', 0)}/{summary.get('num_steps', 0)}")
        print(f"  RMSE:  {summary.get('position_rmse', 0):.10f}")
        print(f"  Cost:  {summary.get('trajectory_cost', 0):.10f}")
        return 0
    else:
        print(f"\nALGORITHM FAILURE: Case '{summary.get('case_id', 'N/A')}' failed.")
        reason = summary.get("failure_reason", "Unknown")
        if reason:
            print(f"  Reason: {reason[:300]}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
