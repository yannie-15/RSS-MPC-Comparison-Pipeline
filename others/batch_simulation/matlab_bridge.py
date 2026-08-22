"""MATLAB CLI bridge for RSS-MPC pipeline.

通过 `matlab -batch` 调用 MATLAB 脚本 (main.m / run_one_case), 替代
MATLAB Engine for Python。MATLAB Engine 在部分 Windows 环境下启动卡死,
而 `matlab -batch` 方式更稳定且兼容性更好。
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from .config_io import save_resolved_config
from .result_io import parse_summary, print_summary


class MatlabBridge:
    """通过 `matlab -batch` 调用 MATLAB 脚本。

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
        self.matlab_exe: str | None = None
        self._started = False

    def start(self) -> None:
        """检测 matlab 可执行文件是否在 PATH 中。

        Raises:
            RuntimeError: 如果找不到 matlab 命令。
        """
        self.matlab_exe = shutil.which("matlab")
        if not self.matlab_exe:
            raise RuntimeError(
                "matlab 命令未找到。请确保 MATLAB 已安装且其 bin 目录在 PATH 中。\n"
                "Windows: 通常为 C:\\Program Files\\MATLAB\\R20XXx\\bin"
            )

        self._started = True
        print("[matlab_bridge] MATLAB CLI ready.")

    def _build_matlab_cmd(self, body: str) -> str:
        """构造 MATLAB -batch 命令字符串。

        Args:
            body: MATLAB 代码主体 (不含 exit, 由本函数添加)。

        Returns:
            完整的 MATLAB 命令字符串。
        """
        # 将仓库根目录、本目录 (main.m 所在) 和 batch_simulation/ (run_one_case 等) 加入 path
        repo_str = str(self.repo_root).replace("\\", "/")
        self_dir_str = str(Path(__file__).parent).replace("\\", "/")
        batch_str = str(self.repo_root / "batch_simulation").replace("\\", "/")
        return (
            f"addpath('{repo_str}', '{self_dir_str}', '{batch_str}'); "
            f"{body}; exit"
        )

    def run_one_case(self, config_path: str | Path) -> dict[str, Any]:
        """通过 matlab -batch 运行单个 case。

        Args:
            config_path: JSON 配置文件路径。

        Returns:
            解析后的 summary 字典。

        Raises:
            RuntimeError: 如果 bridge 未启动。
            FileNotFoundError: 配置文件不存在。
        """
        if not self._started or self.matlab_exe is None:
            raise RuntimeError("MATLAB bridge not started. Call bridge.start() first.")

        config_path = Path(config_path).resolve()
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        print(f"[matlab_bridge] Running case: {config_path}")

        # run_one_case 返回 summary 结构体, 通过 JSON 文件传递结果
        result_json = config_path.parent / f"_result_{config_path.stem}.json"
        if result_json.exists():
            result_json.unlink()

        config_str = str(config_path).replace("\\", "/")
        result_str = str(result_json).replace("\\", "/")
        batch_str = str(self.repo_root / "batch_simulation").replace("\\", "/")
        body = (
            f"addpath('{batch_str}'); "
            f"summary = run_one_case('{config_str}'); "
            f"fid = fopen('{result_str}', 'w'); fwrite(fid, jsonencode(summary)); fclose(fid)"
        )
        matlab_cmd = self._build_matlab_cmd(body)

        result = subprocess.run(
            [self.matlab_exe, "-batch", matlab_cmd],
            cwd=str(self.repo_root),
        )

        if result.returncode != 0:
            raise RuntimeError(f"MATLAB run_one_case 失败 (exit code: {result.returncode})")

        if not result_json.exists():
            raise RuntimeError("MATLAB 未生成结果 JSON 文件")

        with open(result_json, "r", encoding="utf-8") as f:
            raw = json.load(f)
        result_json.unlink()

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
            RuntimeError: 如果 bridge 未启动。
        """
        if not self._started or self.matlab_exe is None:
            raise RuntimeError("MATLAB bridge not started. Call bridge.start() first.")

        # 构造 MATLAB 参数
        seeds_matlab = "[" + ",".join(str(s) for s in seeds) + "]"
        alg_matlab = "{" + ",".join(f"'{a}'" for a in algorithms) + "}"
        force_regen_matlab = "true" if force_regen else "false"

        body = (
            f"main('seeds', {seeds_matlab}, "
            f"'algorithms', {alg_matlab}, "
            f"'forceRegen', {force_regen_matlab})"
        )
        matlab_cmd = self._build_matlab_cmd(body)

        print(f"[matlab_bridge] 调用 matlab -batch: seeds={seeds[0]}:{seeds[-1]} "
              f"(共 {len(seeds)} 个), algorithms={algorithms}, forceRegen={force_regen}")

        try:
            # matlab -batch 直接输出到终端 stdout/stderr, 实时可见
            result = subprocess.run(
                [self.matlab_exe, "-batch", matlab_cmd],
                cwd=str(self.repo_root),
            )
            if result.returncode == 0:
                print("[matlab_bridge] main.m 执行完成。")
                return True
            else:
                print(f"[matlab_bridge] main.m 执行失败 (exit code: {result.returncode})")
                return False
        except Exception as e:
            print(f"[matlab_bridge] main.m 执行异常: {e}")
            return False

    def close(self) -> None:
        """清理资源 (subprocess 方式无需清理)。"""
        self._started = False
        print("[matlab_bridge] MATLAB bridge closed.")

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


def check_matlab_engine() -> bool:
    """检测 MATLAB 是否可用 (matlab CLI 在 PATH 中)。

    Returns:
        True 如果 matlab 命令可用。
    """
    return shutil.which("matlab") is not None


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
