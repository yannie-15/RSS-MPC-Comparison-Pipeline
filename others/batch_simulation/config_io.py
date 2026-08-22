"""Config I/O utilities for RSS-MPC pipeline.

Handles loading, merging, and saving JSON configuration files.
"""

import json
from pathlib import Path
from copy import deepcopy
from typing import Any


DEFAULT_CONFIG: dict[str, Any] = {
    "case_id": "paper_fixed_001",
    "seed": 1,
    "trajectory_mode": "paper_fixed",
    "solver": "sdpt3",
    "num_steps": 100,
    "live_plot": False,
    "save_figures": False,
    "save_full_log": False,
    "output_dir": "results/cases/paper_fixed_001",
}


def load_config(config_path: str | Path | None = None) -> dict[str, Any]:
    """Load a JSON config file, falling back to defaults.

    Args:
        config_path: Path to JSON config file. If None, returns defaults.

    Returns:
        Merged configuration dictionary.
    """
    cfg = deepcopy(DEFAULT_CONFIG)
    if config_path is None:
        return cfg

    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with open(path, "r", encoding="utf-8") as f:
        user_cfg = json.load(f)

    cfg.update(user_cfg)
    return cfg


def apply_cli_overrides(cfg: dict[str, Any], **overrides: Any) -> dict[str, Any]:
    """Apply CLI argument overrides to a config dict.

    Only non-None values are applied.

    Args:
        cfg: Base configuration.
        **overrides: CLI override values (case_id, seed, output, etc.)

    Returns:
        Updated configuration dict (mutates and returns input).
    """
    cli_map = {
        "case_id": "case_id",
        "seed": "seed",
        "output": "output_dir",
        "live_plot": "live_plot",
        "save_figures": "save_figures",
        "force": "force",
    }

    for cli_key, cfg_key in cli_map.items():
        val = overrides.get(cli_key)
        if val is not None:
            cfg[cfg_key] = val

    return cfg


def save_resolved_config(cfg: dict[str, Any], output_dir: str | Path) -> Path:
    """Save the resolved configuration to the output directory.

    Args:
        cfg: Configuration dict.
        output_dir: Output directory path.

    Returns:
        Path to the saved config file.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    config_path = out / "resolved_config.json"
    with open(config_path, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2, ensure_ascii=False)
    return config_path
