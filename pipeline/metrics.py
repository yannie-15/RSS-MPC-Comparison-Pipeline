"""评估指标计算.

复刻两处 MATLAB 逻辑 (合并):
    - core/computeMetrics.m            : 跟踪精度/约束违反率/平滑度/路径信息
    - others/paper_reproduction/run_paper_baseline_case.m 指标段:
        J 分解 (J_position/J_heading/J_control/J_total),
        warm-up 排除后的耗时统计, successRate/validStepRate,
        约束违反量汇总 (原始双线性约束)

J_total = sum_k [ w_pos*(ex^2+ey^2) + w_psi*e_psi^2 + w_control*||u_k||^2 ]
权重取 params.weights (与 MATLAB 硬编码 [30,30,1]/0.3 一致, 因 config 默认相同).
"""

import numpy as np

from pipeline.params import AlgorithmParams
from pipeline.simulator import SimResult


def wrap_angle(a):
    """角度归一化到 (-pi, pi] (MATLAB mod 语义)."""
    return np.mod(a + np.pi, 2.0 * np.pi) - np.pi


def compute_metrics(result: SimResult, params: AlgorithmParams) -> dict:
    """闭环仿真结果 -> 完整指标 dict (含求解时长统计)."""
    path = result.path
    states = result.states
    num_ref = path.shape[1]
    num_steps = result.solvedCount
    planned_steps = result.num_steps
    wheelSpeeds = result.wheelSpeeds
    wheelAngles = result.wheelAngles
    executedU = result.executedU
    solveTimes = result.solveTimes

    vimax = float(params.vehicle.vimax)
    phidotmax = float(params.vehicle.phidotmax)
    dt = float(params.dt)
    w_pos = float(params.weights.w_pos)
    w_psi = float(params.weights.w_psi)
    w_control = float(params.weights.w_control)

    metrics = {}

    # ================= 空仿真 (第一步即失败) =================
    if num_steps == 0:
        metrics.update({
            'rmse': float('nan'), 'maxPositionError': float('nan'),
            'meanPositionError': float('nan'), 'finalPositionError': float('nan'),
            'rmseOrientation': float('nan'), 'maxOrientationError': float('nan'),
            'meanSolveTime': float('nan'), 'maxSolveTime': float('nan'),
            'totalSolveTime': float('nan'),
            'wheelSpeedViolationRatio': float('nan'),
            'maxWheelSpeedViolation': float('nan'),
            'steeringRateViolationRatio': float('nan'),
            'maxSteeringRateViolation': float('nan'),
            'meanWheelSpeedChange': float('nan'),
            'maxWheelSpeedChange': float('nan'),
            'pathLength': _path_length(path), 'numSteps': 0,
            'trajectoryCost': float('nan'),
            'J_position': float('nan'), 'J_heading': float('nan'),
            'J_control': float('nan'), 'J_total': float('nan'),
            'successRate': 0.0, 'validStepRate': 0.0,
            'medianSolveTime': float('nan'), 'q1SolveTime': float('nan'),
            'q3SolveTime': float('nan'), 'meanSolveTimePostWarmup': float('nan'),
            'warmupExcluded': 0,
            'maxWheelViol': 0.0, 'maxConeViol': 0.0,
            'cntWheelViol': 0, 'cntConeViol': 0,
            'cntStrictIncumbent': 0, 'cntApproximateIncumbent': 0,
            'cntNoneIncumbent': 0,
        })
        return metrics

    # ================= 1. 位置跟踪精度 =================
    # states[:, i] 是第 i 步控制前的状态, 对应参考 path[:, i] (1-based 语义)
    position_errors = np.zeros(num_steps)
    for i in range(num_steps):
        ref_idx = min(i + 1, num_ref)
        position_errors[i] = np.linalg.norm(
            states[0:2, i] - path[0:2, ref_idx - 1])

    metrics['rmse'] = float(np.sqrt(np.mean(position_errors ** 2)))
    metrics['maxPositionError'] = float(np.max(position_errors))
    metrics['meanPositionError'] = float(np.mean(position_errors))
    metrics['finalPositionError'] = float(position_errors[-1])

    # ================= 2. 姿态跟踪精度 =================
    orientation_errors = np.zeros(num_steps)
    for i in range(num_steps):
        ref_idx = min(i + 1, num_ref)
        orientation_errors[i] = abs(wrap_angle(states[2, i] - path[2, ref_idx - 1]))

    metrics['rmseOrientation'] = float(np.sqrt(np.mean(orientation_errors ** 2)))
    metrics['maxOrientationError'] = float(np.max(orientation_errors))

    # ================= 3. 求解时间 (computeMetrics 语义) =================
    valid_times = solveTimes[solveTimes > 0] if solveTimes.size else solveTimes
    if valid_times.size > 0:
        metrics['meanSolveTime'] = float(np.mean(valid_times))
        metrics['maxSolveTime'] = float(np.max(valid_times))
        metrics['totalSolveTime'] = float(np.sum(valid_times))
    else:
        metrics['meanSolveTime'] = float('nan')
        metrics['maxSolveTime'] = float('nan')
        metrics['totalSolveTime'] = float('nan')

    # ================= 4. 约束违反率 =================
    if wheelSpeeds is not None and wheelSpeeds.size > 0:
        speed_violations = int(np.sum(wheelSpeeds > vimax * (1 + 1e-6)))
        metrics['wheelSpeedViolationRatio'] = speed_violations / wheelSpeeds.size
        metrics['maxWheelSpeedViolation'] = float(
            max(0.0, np.max(wheelSpeeds) - vimax))
    else:
        metrics['wheelSpeedViolationRatio'] = float('nan')
        metrics['maxWheelSpeedViolation'] = float('nan')

    if wheelAngles is not None and wheelAngles.shape[1] > 1:
        d_angles = np.diff(wheelAngles, axis=1)
        d_angles = np.mod(d_angles + np.pi, 2 * np.pi) - np.pi
        phidot = d_angles / dt
        phidot_violations = int(np.sum(np.abs(phidot) > phidotmax * (1 + 1e-6)))
        metrics['steeringRateViolationRatio'] = phidot_violations / phidot.size
        metrics['maxSteeringRateViolation'] = float(
            max(0.0, np.max(np.abs(phidot)) - phidotmax))
    else:
        metrics['steeringRateViolationRatio'] = float('nan')
        metrics['maxSteeringRateViolation'] = float('nan')

    # ================= 5. 控制平滑度 =================
    if wheelSpeeds is not None and wheelSpeeds.shape[1] > 1:
        speed_diff = np.diff(wheelSpeeds, axis=1)
        metrics['meanWheelSpeedChange'] = float(np.mean(np.abs(speed_diff)))
        metrics['maxWheelSpeedChange'] = float(np.max(np.abs(speed_diff)))
    else:
        metrics['meanWheelSpeedChange'] = float('nan')
        metrics['maxWheelSpeedChange'] = float('nan')

    # ================= 6. 路径信息 =================
    metrics['pathLength'] = _path_length(path)
    metrics['numSteps'] = num_steps

    # ================= 7. J 分解 (式 21 轨迹代价) =================
    # J_total = sum_k [ w_pos*(ex^2+ey^2) + w_psi*e_psi^2 + w_control*||u_k||^2 ]
    J_position = 0.0
    J_heading = 0.0
    J_control = 0.0
    for k in range(num_steps):
        ref_idx = min(k + 1, num_ref)
        e = states[:, k] - path[:, ref_idx - 1]
        e[2] = wrap_angle(e[2])
        u_k = executedU[:, k]
        J_position += w_pos * (e[0] ** 2 + e[1] ** 2)
        J_heading += w_psi * e[2] ** 2
        J_control += w_control * float(np.sum(u_k ** 2))
    metrics['J_position'] = J_position
    metrics['J_heading'] = J_heading
    metrics['J_control'] = J_control
    metrics['J_total'] = J_position + J_heading + J_control
    metrics['trajectoryCost'] = metrics['J_total']

    # ================= 8. 成功率 =================
    metrics['successRate'] = num_steps / planned_steps if planned_steps > 0 else 0.0
    n_valid_steps = int(np.sum(~result.stepFailed))
    metrics['validStepRate'] = n_valid_steps / max(num_steps, 1)

    # ================= 9. 耗时统计 (warm-up 排除, P1-4) =================
    n_warmup = min(5, max(1, int(np.floor(num_steps * 0.05))))
    if num_steps > n_warmup:
        post_warmup_times = solveTimes[n_warmup:]
        metrics['meanSolveTime'] = float(np.mean(solveTimes))
        metrics['maxSolveTime'] = float(np.max(solveTimes))
        metrics['totalSolveTime'] = float(np.sum(solveTimes))
        metrics['medianSolveTime'] = float(np.median(post_warmup_times))
        metrics['q1SolveTime'] = float(np.percentile(post_warmup_times, 25))
        metrics['q3SolveTime'] = float(np.percentile(post_warmup_times, 75))
        metrics['meanSolveTimePostWarmup'] = float(np.mean(post_warmup_times))
    else:
        metrics['meanSolveTime'] = float(np.mean(solveTimes))
        metrics['maxSolveTime'] = float(np.max(solveTimes))
        metrics['totalSolveTime'] = float(np.sum(solveTimes))
        metrics['medianSolveTime'] = float(np.median(solveTimes))
        metrics['q1SolveTime'] = float(np.percentile(solveTimes, 25))
        metrics['q3SolveTime'] = float(np.percentile(solveTimes, 75))
        metrics['meanSolveTimePostWarmup'] = metrics['meanSolveTime']
    metrics['warmupExcluded'] = n_warmup

    # ================= 10. 原始约束违反量汇总 =================
    wheel_viol = result.wheelViolPerStep
    cone_viol = result.coneViolPerStep
    if wheel_viol.size > 0:
        metrics['maxWheelViol'] = float(np.max(wheel_viol))
        metrics['stepWheelMax'] = int(np.argmax(wheel_viol)) + 1
        metrics['maxConeViol'] = float(np.max(cone_viol))
        metrics['stepConeMax'] = int(np.argmax(cone_viol)) + 1
        metrics['cntWheelViol'] = int(np.sum(wheel_viol > 1e-6))
        metrics['cntConeViol'] = int(np.sum(cone_viol > 1e-6))
        metrics['cntStrictIncumbent'] = sum(
            1 for s in result.incumbentTypes if s == 'strict')
        metrics['cntApproximateIncumbent'] = sum(
            1 for s in result.incumbentTypes if s == 'approximate')
        metrics['cntNoneIncumbent'] = sum(
            1 for s in result.incumbentTypes if s == 'none')
    else:
        metrics['maxWheelViol'] = 0.0
        metrics['stepWheelMax'] = 0
        metrics['maxConeViol'] = 0.0
        metrics['stepConeMax'] = 0
        metrics['cntWheelViol'] = 0
        metrics['cntConeViol'] = 0
        metrics['cntStrictIncumbent'] = 0
        metrics['cntApproximateIncumbent'] = 0
        metrics['cntNoneIncumbent'] = 0

    return metrics


def _path_length(path: np.ndarray) -> float:
    """参考路径长度 (XY 平面折线段长度和)."""
    d = np.diff(path[0:2, :], axis=1)
    return float(np.sum(np.sqrt(np.sum(d ** 2, axis=0))))
