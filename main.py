#!/usr/bin/env python
"""RSS-MPC 批量仿真 Python 入口 (对齐 main.m)

通过 MATLAB Engine 调用 main.m, 实现 Step 0-7 全流程:
  路径设置 → seed 列表 → 算法选择 → 断点续跑 → 批量仿真 → 保存 → 对比图 → Table II

用法:
    python main.py                                              # 默认: 1:10 seeds, 全部 4 种算法
    python main.py --seeds 1:100                                # 自定义 seed 范围
    python main.py --algorithms proposed-3iter                  # 只跑一种算法
    python main.py --seeds 1:50 --algorithms proposed-3iter,e-lmpc --force-regen

退出码:
    0  成功
    2  仿真失败
    3  基础设施错误 (MATLAB Engine 未安装等)
"""

import argparse
import sys
from pathlib import Path

# Windows 下强制 stdout/stderr 为 UTF-8, 避免 MATLAB Engine 中文输出乱码
if sys.platform.startswith('win'):
    try:
        sys.stdout.reconfigure(encoding='utf-8')
        sys.stderr.reconfigure(encoding='utf-8')
    except (AttributeError, OSError):
        pass

# 确保 python/ 包可导入
REPO_ROOT = Path(__file__).parent.resolve()
sys.path.insert(0, str(REPO_ROOT))

from batch_simulation.matlab_bridge import MatlabBridge, check_matlab_engine


DEFAULT_ALGORITHMS = ['e-lmpc', 'active-set', 'interior-point', 'proposed-3iter']


def parse_args() -> argparse.Namespace:
    """解析命令行参数 (对齐 main.m 的 name-value 参数)。"""
    parser = argparse.ArgumentParser(
        description="RSS-MPC 批量仿真 Python 入口 (对齐 main.m)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
退出码:
  0  成功
  2  仿真失败
  3  基础设施错误
        """,
    )
    parser.add_argument(
        '--seeds',
        type=str,
        default='1:10',
        help='种子范围, 格式 start:end (默认 1:10), 也支持单个整数',
    )
    parser.add_argument(
        '--algorithms',
        type=str,
        default=None,
        help='算法列表, 逗号分隔 (默认全部 4 种: e-lmpc,active-set,interior-point,proposed-3iter)',
    )
    parser.add_argument(
        '--force-regen',
        action='store_true',
        default=False,
        help='强制重新生成场景文件',
    )
    return parser.parse_args()


def parse_seeds(seed_str: str) -> list[int]:
    """解析 seed 范围字符串。

    支持两种格式:
      '1:10'  → [1, 2, ..., 10]
      '5'     → [5]
    """
    seed_str = seed_str.strip()
    if ':' not in seed_str:
        return [int(seed_str)]
    parts = seed_str.split(':')
    if len(parts) != 2:
        raise ValueError(f"无效的 seed 范围格式: {seed_str} (应为 start:end)")
    start = int(parts[0])
    end = int(parts[1])
    if end < start:
        raise ValueError(f"seed 范围结束值 {end} 小于起始值 {start}")
    return list(range(start, end + 1))


def main() -> int:
    """主入口: 调用 MATLAB main.m 完成批量仿真。

    Returns:
        退出码: 0 (成功), 2 (仿真失败), 3 (基础设施错误)。
    """
    args = parse_args()

    # 解析 seeds
    try:
        seeds = parse_seeds(args.seeds)
    except (ValueError, TypeError) as e:
        print(f"ERROR: 解析 seed 范围失败: {e}")
        return 3

    # 解析 algorithms
    if args.algorithms:
        algorithms = [a.strip() for a in args.algorithms.split(',') if a.strip()]
        if not algorithms:
            print("ERROR: --algorithms 参数为空")
            return 3
    else:
        algorithms = DEFAULT_ALGORITHMS

    print("=" * 60)
    print("RSS-MPC 批量仿真 (Python 入口 → MATLAB main.m)")
    print("=" * 60)
    print(f"  种子范围: {seeds[0]}:{seeds[-1]} (共 {len(seeds)} 个)")
    print(f"  算法:     {', '.join(algorithms)}")
    print(f"  总组合数: {len(seeds) * len(algorithms)}")
    print(f"  强制重新生成场景: {args.force_regen}")
    print(f"  工作目录: {REPO_ROOT}")
    print("=" * 60)

    # 检查 MATLAB Engine 是否可用
    if not check_matlab_engine():
        print("ERROR: matlab.engine 未安装。")
        print("安装 MATLAB Engine for Python:")
        print(f"  1. 在 MATLAB 中运行 'matlabroot' 获取 MATLAB 根目录")
        print(f"  2. cd <matlabroot>/extern/engines/python")
        print(f"  3. {sys.executable} setup.py install")
        return 3

    # 通过 MATLAB Engine 调用 main.m (Step 0-7 全流程)
    try:
        with MatlabBridge(REPO_ROOT) as bridge:
            success = bridge.run_main(seeds, algorithms, args.force_regen)
    except ImportError as e:
        print(f"ERROR: {e}")
        return 3
    except Exception as e:
        print(f"ERROR: 基础设施失败: {e}")
        return 3

    if success:
        print("\n全部步骤完成! 结果保存在 results/batch/ 目录。")
        return 0
    else:
        print("\n仿真执行失败, 详见上方 MATLAB 日志。")
        return 2


if __name__ == '__main__':
    sys.exit(main())
