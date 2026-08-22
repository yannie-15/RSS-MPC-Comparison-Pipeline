"""结果画图: 轨迹跟踪 / 跟踪误差 / 轮速 / 转向速率 / 求解时长 / 控制输入.

对应 MATLAB plot_results.m / plot_one_algorithm.m 的核心图, 全部落盘 PNG
(matplotlib Agg 后端, 无需显示环境).
"""

from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402

from pipeline.params import AlgorithmParams   # noqa: E402
from pipeline.simulator import SimResult      # noqa: E402
from pipeline.metrics import wrap_angle       # noqa: E402


def plot_results(result: SimResult, params: AlgorithmParams,
                 out_dir, metrics: dict = None) -> list:
    """生成全部结果图并保存到 out_dir, 返回生成的文件路径列表."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    vimax = float(params.vehicle.vimax)
    phidotmax = float(params.vehicle.phidotmax)
    dt = float(params.dt)

    M = result.solvedCount
    t = np.arange(M) * dt          # 第 k 步控制作用时刻 t=(k-1)*dt
    files = []

    tag = f"{result.algorithm}_{result.scenario_name}_K{params.K}"

    # ================= 1. 轨迹跟踪 (XY) =================
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(result.path[0, :], result.path[1, :], '-', color='tab:blue',
            linewidth=2, label='Reference Trajectory')
    ax.plot(result.states[0, :M], result.states[1, :M], 'x-',
            color='tab:red', linewidth=1.0, markersize=4,
            label='Simulation Results')
    ax.plot(result.states[0, 0], result.states[1, 0], 'o', color='tab:green',
            markersize=7, label='Start')
    ax.set_xlabel('x_w (m)')
    ax.set_ylabel('y_w (m)')
    ax.set_title('Trajectory Tracking Performance')
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_aspect('equal', adjustable='datalim')
    f = out_dir / f'trajectory_{tag}.png'
    fig.tight_layout()
    fig.savefig(f, dpi=150)
    plt.close(fig)
    files.append(f)

    if M == 0:
        return files

    # ================= 2. 跟踪误差 =================
    num_ref = result.path.shape[1]
    pos_err = np.zeros(M)
    ori_err = np.zeros(M)
    for i in range(M):
        ref_idx = min(i + 1, num_ref)
        pos_err[i] = np.linalg.norm(
            result.states[0:2, i] - result.path[0:2, ref_idx - 1])
        ori_err[i] = abs(wrap_angle(result.states[2, i] - result.path[2, ref_idx - 1]))

    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    axes[0].plot(t, pos_err, '-', color='tab:blue', linewidth=1.5)
    axes[0].set_ylabel('Position error (m)')
    axes[0].set_title('Tracking Errors')
    axes[0].grid(True)
    axes[1].plot(t, ori_err, '-', color='tab:orange', linewidth=1.5)
    axes[1].set_ylabel('Orientation error (rad)')
    axes[1].set_xlabel('Time (s)')
    axes[1].grid(True)
    f = out_dir / f'tracking_errors_{tag}.png'
    fig.tight_layout()
    fig.savefig(f, dpi=150)
    plt.close(fig)
    files.append(f)

    # ================= 3. 轮速 (含 vimax 约束线) =================
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for n in range(result.numWheels):
        ax.plot(t, result.wheelSpeeds[n, :], '-', linewidth=1.5,
                label=f'Wheel {n + 1}')
    ax.axhline(vimax, linestyle='--', color='k', linewidth=1.5,
               label='Constraints')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Output Velocity (m/s)')
    ax.set_title('Wheel Speeds')
    ax.legend(fontsize=7)
    ax.grid(True)
    f = out_dir / f'wheel_speeds_{tag}.png'
    fig.tight_layout()
    fig.savefig(f, dpi=150)
    plt.close(fig)
    files.append(f)

    # ================= 4. 转向速率 (含 ±phidotmax 约束线) =================
    if result.wheelAngles.shape[1] > 1:
        d_angles = np.diff(result.wheelAngles, axis=1)
        d_angles = np.mod(d_angles + np.pi, 2 * np.pi) - np.pi
        phidot = d_angles / dt
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for n in range(result.numWheels):
            ax.plot(t[1:], phidot[n, :], '-', linewidth=1.5,
                    label=f'Wheel {n + 1}')
        ax.axhline(phidotmax, linestyle='--', color='k', linewidth=1.5,
                   label='Constraints')
        ax.axhline(-phidotmax, linestyle='--', color='k', linewidth=1.5)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Steering rate (rad/s)')
        ax.set_title('Steering Rate of Each Wheel')
        ax.legend(fontsize=7)
        ax.grid(True)
        f = out_dir / f'steering_rates_{tag}.png'
        fig.tight_layout()
        fig.savefig(f, dpi=150)
        plt.close(fig)
        files.append(f)

    # ================= 5. 求解时长 =================
    if result.solveTimes.size > 0:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(t, result.solveTimes, '.-', color='tab:purple',
                linewidth=1.0, markersize=4, label='Per-step solve time')
        if metrics is not None and metrics.get('medianSolveTime') is not None:
            med = metrics.get('medianSolveTime')
            if med == med:   # not NaN
                ax.axhline(med, linestyle='--', color='tab:red', linewidth=1.5,
                           label=f'Median (post warm-up) = {med:.4f} s')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Solve time (s)')
        ax.set_title('Per-step Solver Time')
        ax.legend(fontsize=8)
        ax.grid(True)
        f = out_dir / f'solve_times_{tag}.png'
        fig.tight_layout()
        fig.savefig(f, dpi=150)
        plt.close(fig)
        files.append(f)

    # ================= 6. 控制输入 =================
    fig, ax = plt.subplots(figsize=(7, 4.5))
    labels = ['$u_{v_x}$', '$u_{v_y}$', '$u_{\\omega}$']
    for i in range(3):
        ax.plot(t, result.executedU[i, :], '-', linewidth=1.5, label=labels[i])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Control increment')
    ax.set_title('Executed Control Inputs')
    ax.legend(fontsize=8)
    ax.grid(True)
    f = out_dir / f'controls_{tag}.png'
    fig.tight_layout()
    fig.savefig(f, dpi=150)
    plt.close(fig)
    files.append(f)

    return files
