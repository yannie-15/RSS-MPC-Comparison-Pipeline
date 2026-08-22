"""闭环仿真器: 原始动力学 + 循环调用算法 (4 算法统一入口).

为 others/paper_reproduction/run_paper_baseline_case.m 闭环部分的 Python 等价移植:
    - 原始动力学推进 (core/propagateState.m) 与轮子正运动学
      (core/computeWheelOutputs.m), 不使用算法内部的展开/线性化误差动力学
    - 每步: control(path, k, lastBodyVelocity, state, params) -> u/世界速度/车体速度/诊断
    - step_failed / NaN-Inf 检查失败即终止 (与 MATLAB 一致)
    - per-step 记录: 状态/速度/控制/轮速/轮角/求解时长/迭代数/约束违反量/诊断

算法通过注册表分发, 4 种算法全部在本仿真器内闭环:
    - proposed-3iter : 纯 Python (HPIPM OCP QCQP + SCP), 可用 --iters 指定
                       SCP 外层迭代次数 (默认 3, 即论文基准 proposed-3iter)
    - e-lmpc / active-set / interior-point : MATLAB 子模块算法 (fmincon 系列),
                       经 matlab_algorithm.MatlabAlgorithmBridge (常驻 MATLAB
                       Engine 会话) 每步调用 control_RSS 求解 —— 轨迹/参数/闭环
                       推进/评估全在本仿真器 (Python), MATLAB 只做算法求解
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional
import sys
import numpy as np

from pipeline.params import AlgorithmParams, DT
from pipeline.dynamics import propagate_state, compute_wheel_outputs
from pipeline.trajectory_generator import Trajectory

# proposed 控制器实现归位算法包 algorithms/RSS_proposed/ (与 hpipm_qp_solver/
# construct_ocp_qcqp 同目录); pipeline/ 只留轨迹/仿真/评估主干
_RSS_PROPOSED_DIR = str(Path(__file__).resolve().parent.parent
                        / 'algorithms' / 'RSS_proposed')
if _RSS_PROPOSED_DIR not in sys.path:
    sys.path.insert(0, _RSS_PROPOSED_DIR)
from control_rss_ocpqcqp import control_rss_ocpqcqp

# 算法注册表: 算法名 -> 控制函数 (None 表示该算法未在 Python pipeline 移植)
# 统一签名: control(path, step, state_dot, state, params, verbose)
#           -> (u(3,K), world_velocity(3,), body_velocity(3,), diagnostics dict)
MATLAB_ONLY_ALGORITHMS = ('e-lmpc', 'active-set', 'interior-point')
ALGORITHMS: Dict[str, Optional[Callable]] = {
    'proposed-3iter': control_rss_ocpqcqp,
    'e-lmpc': None,
    'active-set': None,
    'interior-point': None,
}


@dataclass
class SimResult:
    """闭环仿真输出 (与 MATLAB summary struct 对应)."""
    algorithm: str
    scenario_name: str
    seed_id: int
    path: np.ndarray                    # (3, N) 参考轨迹
    num_steps: int                      # 计划步数
    dt: float
    states: np.ndarray = None           # (3, M+1) 世界系位姿序列
    worldVelocities: np.ndarray = None  # (3, M)
    bodyVelocities: np.ndarray = None   # (3, M)
    executedU: np.ndarray = None        # (3, M)
    wheelSpeeds: np.ndarray = None      # (num_wheels, M)
    wheelAngles: np.ndarray = None      # (num_wheels, M)
    solveTimes: np.ndarray = None       # (M,) 每步 SCP 总求解时长
    iterations: np.ndarray = None       # (M,) 真实 solver 调用数
    stepFailed: np.ndarray = None       # (M,) bool
    stepFinite: np.ndarray = None       # (M,) bool
    stepApproximate: np.ndarray = None  # (M,) bool
    incumbentTypes: List[str] = field(default_factory=list)
    wheelViolPerStep: np.ndarray = None  # (M,) 原始轮速约束违反
    coneViolPerStep: np.ndarray = None   # (M,) 原始转向锥约束违反
    stepDiagnostics: List[dict] = field(default_factory=list)  # 每步完整诊断
    success: bool = True
    failureReason: str = ''
    firstFailStep: int = 0
    solvedCount: int = 0

    @property
    def numWheels(self) -> int:
        return int(self.wheelSpeeds.shape[0])


def simulate(trajectory: Trajectory, algorithm: str,
             params: Optional[AlgorithmParams] = None,
             verbose: bool = True) -> SimResult:
    """闭环仿真入口.

    输入:
        trajectory : trajectory_generator.make_trajectory 输出 (含场景/轮子参数)
        algorithm  : 算法名 (4 选 1, 见模块 docstring)
        params     : AlgorithmParams (None 则默认; vehicle 会被场景参数覆盖)
        verbose    : 透传给控制函数的打印开关

    仿真器使用原始动力学 (propagate_state / compute_wheel_outputs) 推进闭环,
    算法只提供每步控制增量; 场景的轮子参数 (wheel_pos/vimax/phidotmax)
    覆盖 params.vehicle, 保证与参考场景一致.
    MATLAB 算法 (e-lmpc/active-set/interior-point) 经 MatlabAlgorithmBridge
    每步求解: 引擎会话随本函数开启/关闭 (启动约 30s, 之后每步毫秒级调用).
    """
    algorithm = algorithm.lower()
    if algorithm not in ALGORITHMS:
        raise ValueError(
            f'不支持算法: {algorithm} (可用: {list(ALGORITHMS)})')
    control_fn = ALGORITHMS[algorithm]
    bridge = None
    # Ensure params defaulted before possibly using it to initialize MATLAB bridge
    if params is None:
        params = AlgorithmParams(K=6)

    # 场景参数覆盖: 轨迹生成器给出的车辆 (轮子) 参数优先 (CLI --vimax/--phidotmax
    # 的覆盖已写入 trajectory.vehicle; 必须在 MATLAB 桥初始化前生效, 使
    # run_config 携带覆盖后的约束值; rho_schedule 一并透传, 供逐步 rho 调度)
    params = AlgorithmParams(
        K=params.K, dt=params.dt, vehicle=trajectory.vehicle,
        weights=params.weights, solver=params.solver, max_iter=params.max_iter,
        rho_schedule=params.rho_schedule)

    if control_fn is None:
        if algorithm not in MATLAB_ONLY_ALGORITHMS:
            raise ValueError(
                f'算法 {algorithm} 注册表项为空且非 MATLAB 算法, 请检查 ALGORITHMS')
        from pipeline.matlab_algorithm import MatlabAlgorithmBridge
        # 将完整 run config 传给 MATLAB 桥，使 MATLAB 端可使用与 Python 一致的参数
        run_config = {
            'seed_id': int(trajectory.seed_id),
            'algorithm': algorithm,
            'K': int(params.K),
            'dt': float(params.dt),
            'num_steps': int(trajectory.num_steps),
            'trajectory_source': trajectory.source,
            'scenario_name': trajectory.scenario_name,
            'algorithm_params': params.to_dict(),
        }
        bridge = MatlabAlgorithmBridge(algorithm, trajectory, verbose=verbose, K=params.K, run_config=run_config)
        control_fn = bridge.control

    dt = float(params.dt)
    num_steps = int(trajectory.num_steps)
    num_wheels = trajectory.vehicle.num_wheels
    wheel_pos = np.asarray(trajectory.vehicle.wheel_pos, dtype=np.float64)

    path = np.asarray(trajectory.path, dtype=np.float64)

    state = np.asarray(trajectory.initial_state, dtype=np.float64).flatten()
    last_body_velocity = np.asarray(trajectory.initial_velocity,
                                    dtype=np.float64).flatten()

    # ================= 预分配 =================
    states = np.zeros((3, num_steps + 1))
    worldVelocities = np.zeros((3, num_steps))
    bodyVelocities = np.zeros((3, num_steps))
    wheelSpeeds = np.zeros((num_wheels, num_steps))
    wheelAngles = np.zeros((num_wheels, num_steps))
    executedU = np.zeros((3, num_steps))

    solveTimes = np.zeros(num_steps)
    stepIterations = np.zeros(num_steps, dtype=np.int64)
    stepFinite = np.ones(num_steps, dtype=bool)
    stepFailed = np.zeros(num_steps, dtype=bool)
    stepApproximate = np.zeros(num_steps, dtype=bool)
    stepSolverCallCount = np.zeros(num_steps, dtype=np.int64)
    stepWheelViol = np.zeros(num_steps)
    stepConeViol = np.zeros(num_steps)
    stepIncumbentType: List[str] = [''] * num_steps
    stepDiagnostics: List[dict] = []

    success = True
    failure_reason = ''
    solved_count = 0
    first_fail_step = 0

    states[:, 0] = state

    # ================= 闭环仿真 =================
    for k in range(1, num_steps + 1):
        try:
            # ---- 调用算法 (传入 path/步号/车体速度/位姿/参数) ----
            u_full, world_velocity, body_velocity, diagnostics = control_fn(
                path, k, last_body_velocity, state, params, verbose)

            u = u_full[:, 0]
            solve_time = float(diagnostics.get('total_solve_time', float('nan')))
            iter_num = int(diagnostics.get('solver_call_count', 0))

            solveTimes[k - 1] = solve_time
            stepIterations[k - 1] = iter_num
            stepSolverCallCount[k - 1] = iter_num
            stepDiagnostics.append(diagnostics)

            # ---- 记录 per-step 约束违反量 ----
            if 'orig_wheel_viol_final' in diagnostics:
                stepWheelViol[k - 1] = diagnostics['orig_wheel_viol_final']
                stepConeViol[k - 1] = diagnostics['orig_cone_viol_final']
            if 'incumbent_type' in diagnostics:
                stepIncumbentType[k - 1] = diagnostics['incumbent_type']
            if diagnostics.get('step_approximate', False):
                stepApproximate[k - 1] = True

            # ---- 检查 step_failed: 不得忽略 ----
            if diagnostics.get('step_failed', False):
                stepFailed[k - 1] = True
                stepFinite[k - 1] = False
                if first_fail_step == 0:
                    first_fail_step = k
                success = False
                failure_reason = ('diagnostics.step_failed=true '
                                  '(三次 outer 后无 incumbent)')
                print(f'[{algorithm}: step {k}] STEP FAILED: {failure_reason}')
                break

            # ---- 检查 solver_call_count == max_iter (仅 proposed SCP; MATLAB 算法无此口径) ----
            if algorithm not in MATLAB_ONLY_ALGORITHMS and \
                    iter_num != params.max_iter:
                print(f'[{algorithm}: step {k}] 警告: solver_call_count='
                      f'{iter_num} (预期 {params.max_iter})')

            # ---- 解有限性检查 ----
            step_ok = True
            fail_reason = ''
            if not np.all(np.isfinite(world_velocity)):
                step_ok = False
                fail_reason = 'worldVelocity contains NaN/Inf'
            if not np.all(np.isfinite(u)):
                step_ok = False
                if not fail_reason:
                    fail_reason = 'u contains NaN/Inf'
            stepFinite[k - 1] = step_ok

            if not step_ok:
                stepFailed[k - 1] = True
                print(f'[Step {k}: {algorithm}] 求解器诊断: {fail_reason}')
                # 仍继续执行, 但记录为失败步

            # ---- 执行控制: 轮子正运动学 (原始动力学) ----
            wheel_speed, wheel_angle = compute_wheel_outputs(body_velocity,
                                                             wheel_pos)

            executedU[:, k - 1] = u
            states[:, k - 1] = state
            worldVelocities[:, k - 1] = world_velocity
            bodyVelocities[:, k - 1] = body_velocity
            wheelSpeeds[:, k - 1] = wheel_speed
            wheelAngles[:, k - 1] = wheel_angle
            solved_count += 1

            # ---- 推进状态 (原始动力学) ----
            state = propagate_state(state, world_velocity, dt)
            last_body_velocity = body_velocity
            states[:, k] = state

        except Exception as exc:   # noqa: BLE001 - 与 MATLAB catch ME 一致
            import traceback
            stepFailed[k - 1] = True
            stepFinite[k - 1] = False
            if first_fail_step == 0:
                first_fail_step = k
            success = False
            failure_reason = ''.join(
                traceback.format_exception_only(type(exc), exc)).strip()
            print(f'[{algorithm}: step {k}] EXCEPTION: {exc}')
            traceback.print_exc()
            break

    # ================= 关闭 MATLAB Engine 会话 (MATLAB 算法) =================
    if bridge is not None:
        bridge.close()

    # ================= 截断到实际完成步数 =================
    if solved_count == 0:
        states = states[:, :1]
        worldVelocities = np.zeros((3, 0))
        bodyVelocities = np.zeros((3, 0))
        executedU = np.zeros((3, 0))
        wheelSpeeds = np.zeros((num_wheels, 0))
        wheelAngles = np.zeros((num_wheels, 0))
        solveTimes = np.zeros(0)
        stepIterations = stepIterations[:0]
        stepFinite = stepFinite[:0]
        stepFailed = stepFailed[:0]
        stepApproximate = stepApproximate[:0]
        stepWheelViol = stepWheelViol[:0]
        stepConeViol = stepConeViol[:0]
        stepIncumbentType = []
    else:
        states = states[:, :solved_count + 1]
        worldVelocities = worldVelocities[:, :solved_count]
        bodyVelocities = bodyVelocities[:, :solved_count]
        executedU = executedU[:, :solved_count]
        wheelSpeeds = wheelSpeeds[:, :solved_count]
        wheelAngles = wheelAngles[:, :solved_count]
        solveTimes = solveTimes[:solved_count]
        stepIterations = stepIterations[:solved_count]
        stepFinite = stepFinite[:solved_count]
        stepFailed = stepFailed[:solved_count]
        stepApproximate = stepApproximate[:solved_count]
        stepWheelViol = stepWheelViol[:solved_count]
        stepConeViol = stepConeViol[:solved_count]
        stepIncumbentType = stepIncumbentType[:solved_count]

    n_valid_steps = int(np.sum(~stepFailed)) if solved_count > 0 else 0

    return SimResult(
        algorithm=algorithm,
        scenario_name=trajectory.scenario_name,
        seed_id=trajectory.seed_id,
        path=path,
        num_steps=num_steps,
        dt=dt,
        states=states,
        worldVelocities=worldVelocities,
        bodyVelocities=bodyVelocities,
        executedU=executedU,
        wheelSpeeds=wheelSpeeds,
        wheelAngles=wheelAngles,
        solveTimes=solveTimes,
        iterations=stepIterations,
        stepFailed=stepFailed,
        stepFinite=stepFinite,
        stepApproximate=stepApproximate,
        incumbentTypes=stepIncumbentType,
        wheelViolPerStep=stepWheelViol,
        coneViolPerStep=stepConeViol,
        stepDiagnostics=stepDiagnostics,
        success=success,
        failureReason=failure_reason,
        firstFailStep=first_fail_step,
        solvedCount=solved_count,
    )
