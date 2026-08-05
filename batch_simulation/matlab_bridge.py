"""MATLAB Engine bridge for RSS-MPC pipeline.

Starts a single MATLAB Engine instance, adds the repository to the MATLAB path,
and calls run_one_case once per case (not per MPC step).
"""

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from .config_io import save_resolved_config
from .result_io import parse_summary, print_summary


class MatlabBridge:
    """Manages a single MATLAB Engine session for running RSS cases.

    Usage:
        bridge = MatlabBridge(repo_root)
        bridge.start()
        summary = bridge.run_one_case(config_path)
        bridge.close()

    Or as a context manager:
        with MatlabBridge(repo_root) as bridge:
            summary = bridge.run_one_case(config_path)
    """

    def __init__(self, repo_root: str | Path):
        self.repo_root = Path(repo_root).resolve()
        self.engine: Any = None
        self._started = False

    def start(self) -> None:
        """Start the MATLAB Engine and add repo to path.

        Raises:
            ImportError: If matlab.engine is not installed.
            RuntimeError: If the engine fails to start.
        """
        try:
            import matlab.engine
        except ImportError as e:
            raise ImportError(
                "matlab.engine not found. Install MATLAB Engine for Python:\n"
                f"  cd <matlabroot>/extern/engines/python\n"
                f"  python setup.py install\n"
                f"Original error: {e}"
            ) from e

        print("[matlab_bridge] Starting MATLAB Engine...")
        self.engine = matlab.engine.start_matlab()

        # Windows 下 MATLAB 默认用系统 ANSI 编码 (GBK/cp936) 输出中文,
        # 传给 Python 时按 UTF-8 解码会乱码; 强制设为 UTF-8
        if sys.platform.startswith('win'):
            try:
                self.engine.feature('DefaultCharacterSet', 'UTF-8', nargout=0)
            except Exception:
                pass

        # Add repo root to MATLAB path (highest priority)
        # genpath adds root first, then subdirs, so root files take precedence
        self.engine.addpath(
            str(self.repo_root),
            nargout=0,
        )
        # 设置 MATLAB 工作目录为 batch_simulation/ (main.m 所在目录)
        batch_dir = str(self.repo_root / "batch_simulation")
        if Path(batch_dir).exists():
            self.engine.cd(batch_dir, nargout=0)
            self.engine.addpath(batch_dir, nargout=0)
        # 添加 core/ (共享工具) 和 algorithms/ (算法包)
        core_dir = str(self.repo_root / "core")
        if Path(core_dir).exists():
            self.engine.addpath(core_dir, nargout=0)
        alg_dir = str(self.repo_root / "algorithms")
        if Path(alg_dir).exists():
            self.engine.addpath(alg_dir, nargout=0)

        self._started = True
        print("[matlab_bridge] MATLAB Engine ready.")

    def run_one_case(self, config_path: str | Path) -> dict[str, Any]:
        """Run a single case via MATLAB run_one_case.

        Args:
            config_path: Path to the JSON config file.

        Returns:
            Parsed summary dict.

        Raises:
            RuntimeError: If the engine is not started.
            json.JSONDecodeError: If MATLAB returns invalid JSON.
        """
        if not self._started or self.engine is None:
            raise RuntimeError("MATLAB Engine not started. Call bridge.start() first.")

        config_path = Path(config_path).resolve()
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        print(f"[matlab_bridge] Running case: {config_path}")
        raw = self.engine.run_one_case(str(config_path), nargout=1)
        summary = parse_summary(raw)
        print_summary(summary)
        return summary

    def run_main(self, seeds: list[int], algorithms: list[str], force_regen: bool = False) -> bool:
        """调用 MATLAB main.m 批量入口 (Step 0-7 全流程)。

        等价于在 MATLAB 中执行:
            main('seeds', 1:10, 'algorithms', {...}, 'forceRegen', false)

        Args:
            seeds: 种子列表, 如 [1, 2, ..., 10]。
            algorithms: 算法名列表, 如 ['proposed-3iter', 'e-lmpc']。
            force_regen: 是否强制重新生成场景文件。

        Returns:
            True 如果 MATLAB main.m 正常执行完成, False 如果抛出异常。

        Raises:
            RuntimeError: 如果 MATLAB Engine 未启动。
        """
        if not self._started or self.engine is None:
            raise RuntimeError("MATLAB Engine not started. Call bridge.start() first.")

        import matlab

        # Python list[int] → matlab.double (对应 MATLAB 1×N double 数组)
        seeds_matlab = matlab.double(seeds)
        # Python list[str] → MATLAB 1×N cell array of char (main.m 的 parse_args 兼容此格式)
        algorithms_matlab = algorithms

        print(f"[matlab_bridge] 调用 main.m: seeds={seeds[0]}:{seeds[-1]} "
              f"(共 {len(seeds)} 个), algorithms={algorithms}, forceRegen={force_regen}")

        try:
            # main.m 返回 comparison 结构体; 这里不取返回值 (Step 6/7 内部已完成画图/打印)
            # nargout=0 让 MATLAB 不强制要求接收返回值
            self.engine.main(
                'seeds', seeds_matlab,
                'algorithms', algorithms_matlab,
                'forceRegen', bool(force_regen),
                nargout=0,
            )
            print("[matlab_bridge] main.m 执行完成。")
            return True
        except Exception as e:
            err_msg = str(e)
            # MATLAB Engine 在 Windows 上可能仍用 GBK 编码错误信息 (即便已设置
            # DefaultCharacterSet), 尝试 latin-1 → gbk 重新解码还原中文
            if sys.platform.startswith('win'):
                try:
                    err_msg = err_msg.encode('latin-1').decode('gbk')
                except (UnicodeDecodeError, UnicodeEncodeError):
                    pass
            print(f"[matlab_bridge] main.m 执行失败: {err_msg}")
            return False

    def close(self) -> None:
        """Shut down the MATLAB Engine."""
        if self.engine is not None:
            try:
                self.engine.quit()
            except Exception:
                pass
            self.engine = None
            self._started = False
            print("[matlab_bridge] MATLAB Engine closed.")

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


def check_matlab_engine() -> bool:
    """Check if matlab.engine is available.

    Returns:
        True if matlab.engine can be imported.
    """
    try:
        import matlab.engine  # noqa: F401
        return True
    except ImportError:
        return False


def check_matlab_install() -> dict[str, str]:
    """Detect MATLAB installation path and version.

    Returns:
        Dict with 'matlabroot' and 'version' keys, or empty strings if not found.
    """
    info = {"matlabroot": "", "version": ""}

    # Try common Windows locations
    if sys.platform == "win32":
        import winreg
        try:
            with winreg.OpenKey(winreg.HKEY_CLASSES_ROOT, r"MATLAB.Application.1\\shell\\open\\command") as key:
                cmd = winreg.QueryValue(key, None)
                if cmd:
                    info["matlabroot"] = str(Path(cmd.split('"')[1]).parent.parent)
        except (OSError, IndexError):
            pass

    return info


def install_matlab_engine(matlabroot: str | None = None, python_exe: str | None = None) -> int:
    """Guide the user through installing MATLAB Engine for Python.

    Args:
        matlabroot: Path to MATLAB root. If None, attempts auto-detection.
        python_exe: Python executable to use. If None, uses sys.executable.

    Returns:
        Exit code from the install command (0 = success).
    """
    if matlabroot is None:
        info = check_matlab_install()
        matlabroot = info["matlabroot"]

    if not matlabroot or not Path(matlabroot).exists():
        print("ERROR: Could not find MATLAB installation.")
        print("Please provide the MATLAB root path manually.")
        return 1

    if python_exe is None:
        python_exe = sys.executable

    engine_dir = str(Path(matlabroot) / "extern" / "engines" / "python")
    print(f"Installing MATLAB Engine from: {engine_dir}")
    print(f"Using Python: {python_exe}")

    result = subprocess.run(
        [python_exe, "setup.py", "install"],
        cwd=engine_dir,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(f"Install failed (exit {result.returncode}):")
        print(result.stderr)
    else:
        print("MATLAB Engine installed successfully.")

    return result.returncode
