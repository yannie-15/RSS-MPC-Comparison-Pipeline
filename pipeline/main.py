"""Pipeline 入口: 传入 seed_id / 算法 / 步长 K / SCP 迭代数, 生成结果三件套.

用法:
    python pipeline/main.py --seed 0 --algorithm proposed-3iter --K 6
    python pipeline/main.py --seed 0 --algorithm proposed-3iter --K 6 --iters 5
    python pipeline/main.py --seed 0 --algorithm e-lmpc --K 6        (MATLAB Engine)
    python -m pipeline.main --seed 0 --algorithm proposed-3iter --K 6

流程 (4 算法同一条 Python 主干, MATLAB 只做算法求解):
    main(seed_id, algorithm, K, iters)             [dt (tau) 写死 0.01 s]
      -> 轨迹生成: make_trajectory (.py, seed 与 K 传入)
      -> 闭环仿真: simulate (.py, 原始动力学循环调用算法)
           proposed-3iter             -> algorithms/RSS_proposed/control_rss_ocpqcqp.py
           e-lmpc/active-set/interior-point -> matlab_algorithm.py
             (常驻 MATLAB Engine 会话, 每步调 submodule control_RSS 求解)
      -> compute_metrics(result, params)           评估 (含求解时长统计)
      -> 落盘结果:
           run_config.json    1) 记录传入参数 (seed_id/算法/K/iters/dt/场景/权重/求解器)
           metrics.json       2) 评估指标 (RMSE/J_total/求解时长/约束违反...)
           figures/*.png      3) 画图 (轨迹/误差/轮速/转向速率/求解时长/控制)
           simulation_data.npz 原始数组 (复画图/复核用)

算法选择 (4 选 1):
    proposed-3iter    纯 Python (HPIPM OCP QCQP + SCP, 实现在
                      algorithms/RSS_proposed/), --iters 可指定
                      SCP 外层迭代次数 (默认 3 = 论文基准 proposed-3iter)
    e-lmpc / active-set / interior-point
                      MATLAB submodule 算法 (fmincon 系), 经 MATLAB Engine
                      每步求解 (需 matlabengine 包 + MATLAB 在 PATH);
                      K 由 submodule 硬编码为 6, --iters 对其无效
"""

import argparse
import json
import sys
import time
from pathlib import Path

# 支持 `python pipeline/main.py` 与 `python -m pipeline.main` 两种启动方式
if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline.params import DT, AlgorithmParams, RunConfig
from pipeline.trajectory_generator import make_trajectory
from pipeline.simulator import simulate, ALGORITHMS, MATLAB_ONLY_ALGORITHMS
from pipeline.metrics import compute_metrics
from pipeline.plotting import plot_results


def run(seed_id: int, algorithm: str, K: int, iters: int = 3,
        out_root=None, verbose: bool = True) -> dict:
    """单次完整 pipeline 运行, 返回 {'run_dir', 'metrics', 'success', ...}.

    iters: proposed 算法的 SCP 外层迭代数 (每步 HPIPM 求解次数), 默认 3.
    """
    algorithm = algorithm.lower()
    is_matlab_algo = algorithm in MATLAB_ONLY_ALGORITHMS

    if out_root is None:
        out_root = Path(__file__).resolve().parent.parent / 'results' / 'pipeline'
    out_root = Path(out_root)
    run_name = f'seed{seed_id}_K{K}'
    if iters != 3 and not is_matlab_algo:
        run_name += f'_iters{iters}'
    run_dir = out_root / algorithm / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    t_sim0 = time.perf_counter()

    # ============ 统一 Python 主干 (4 算法) ============
    if is_matlab_algo:
        if K != 6:
            print(f'警告: {algorithm} 为 MATLAB submodule 算法, 预测时域硬编码 '
                  f'K=6, --K {K} 被忽略')
        if iters != 3:
            print(f'警告: --iters 仅对 proposed-3iter 生效, {algorithm} 忽略')

    # 1. 轨迹生成 (K 与 seed_id 传入, 4 算法同源)
    trajectory = make_trajectory(seed_id, K, DT)
    # 2. 闭环仿真 (轨迹 + 算法 + 参数 传入仿真器; MATLAB 算法经 Engine 桥每步求解)
    params = AlgorithmParams(K=K, dt=DT, max_iter=iters)
    result = simulate(trajectory, algorithm, params, verbose=verbose)
    traj_num_steps = trajectory.num_steps
    traj_source = trajectory.source
    traj_name = trajectory.scenario_name

    wall_time = time.perf_counter() - t_sim0

    # ================= 3. 评估 =================
    metrics = compute_metrics(result, params)

    # ================= 4. 结果落盘 =================
    # 4.1 运行参数记录 (seed_id / 算法 / K / dt / 场景 / 全部算法参数)
    run_config = RunConfig(
        seed_id=seed_id,
        algorithm=algorithm,
        K=K,
        dt=DT,
        num_steps=traj_num_steps,
        trajectory_source=traj_source,
        scenario_name=traj_name,
        algorithm_params=params,
    )
    config_payload = run_config.to_dict()
    config_payload['run_info'] = {
        'backend': 'matlab-engine' if is_matlab_algo else 'python',
        'wall_time_s': wall_time,
        'solved_count': result.solvedCount,
        'success': result.success,
        'failure_reason': result.failureReason,
        'first_fail_step': result.firstFailStep,
    }
    with open(run_dir / 'run_config.json', 'w', encoding='utf-8') as f:
        json.dump(config_payload, f, ensure_ascii=False, indent=2)

    # 4.2 评估指标 (含求解时长)
    with open(run_dir / 'metrics.json', 'w', encoding='utf-8') as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    # 4.3 原始数组 (复画图 / golden 对比用)
    np_save(run_dir / 'simulation_data.npz', result)

    # 4.4 画图
    figures = plot_results(result, params, run_dir / 'figures', metrics)

    # ================= 5. 汇总打印 (与 MATLAB 汇总行风格一致) =================
    n_valid = int((~result.stepFailed).sum()) if result.solvedCount > 0 else 0
    if result.success:
        print(f'[{result.algorithm}: {result.scenario_name}] '
              f'RMSE={metrics["rmse"]:.6f}, '
              f'medianSolveTime={metrics["medianSolveTime"]:.4f}s, '
              f'J_total={metrics["J_total"]:.4f}, '
              f'validSteps={n_valid}/{result.solvedCount}')
    else:
        print(f'[{result.algorithm}: {result.scenario_name}] '
              f'FAILED at step {result.firstFailStep}/{result.num_steps}: '
              f'{result.failureReason}')
    print(f'结果目录: {run_dir}')

    return {
        'run_dir': str(run_dir),
        'metrics': metrics,
        'figures': [str(p) for p in figures],
        'success': result.success,
        'wall_time_s': wall_time,
    }


def np_save(path, result):
    """SimResult 数组字段存 npz (诊断 dict 不序列化, 只存数值数组)."""
    import numpy as np
    payload = {
        'states': result.states,
        'worldVelocities': result.worldVelocities,
        'bodyVelocities': result.bodyVelocities,
        'executedU': result.executedU,
        'wheelSpeeds': result.wheelSpeeds,
        'wheelAngles': result.wheelAngles,
        'path': result.path,
        'solveTimes': result.solveTimes,
        'iterations': result.iterations,
        'stepFailed': result.stepFailed,
        'stepFinite': result.stepFinite,
        'wheelViolPerStep': result.wheelViolPerStep,
        'coneViolPerStep': result.coneViolPerStep,
        'incumbentTypes': np.array(result.incumbentTypes, dtype=object),
    }
    np.savez(path, **payload)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='RSS MPC pipeline: main(seed_id, algorithm, K, iters), '
                    'dt=0.01s 固定')
    parser.add_argument('--seed', type=int, default=0,
                        help='场景 seed_id (0=paper_fixed 固定场景, >=1=scenario_bank)')
    parser.add_argument('--algorithm', default='proposed-3iter',
                        choices=sorted(ALGORITHMS.keys()),
                        help='算法名 (4 选 1); proposed-3iter 为纯 Python, '
                             'e-lmpc/active-set/interior-point 经 MATLAB Engine '
                             '每步求解 (需 matlabengine 包 + MATLAB 在 PATH, K 固定 6)')
    parser.add_argument('--K', type=int, default=6,
                        help='MPC 预测时域步长 K (默认 6, 论文 IV-A)')
    parser.add_argument('--iters', type=int, default=3,
                        help='proposed 算法 SCP 外层迭代数 (每步 HPIPM 求解次数, '
                             '默认 3 = 论文基准 proposed-3iter)')
    parser.add_argument('--out', default=None,
                        help='结果输出根目录 (默认 results/pipeline)')
    parser.add_argument('--quiet', action='store_true',
                        help='关闭算法 per-iteration 打印')
    args = parser.parse_args(argv)

    # dt (tau) 写死为 0.01 s, 不作为命令行参数
    run(args.seed, args.algorithm, args.K, iters=args.iters,
        out_root=args.out, verbose=not args.quiet)


if __name__ == '__main__':
    main()
