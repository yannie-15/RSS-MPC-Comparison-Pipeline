"""RSS proposed 控制器 (HPIPM OCP QCQP, 精确凸二次约束).

RSS proposed 控制器 (OCP QCQP + SCP) 的唯一实现 (原 MATLAB 版 control_RSS_ocpqcqp.m 已删除).
Golden 对齐语义:
    - u_hat = zeros(3, K) (零初始化, 无 warm start, 无 cut_points)
    - max_iter 次 outer 迭代 (默认 3 = 论文基准 proposed-3iter; 可由
      pipeline/main.py --iters 指定), 每次恰好 1 次 solve_ocp_qcqp
      (solver_call_count=1)
    - status==0 且解有限: u_hat = u_sol
    - 否则: 保留上一轮 u_hat, 不推进
    - 所有成功步标记为 strict (无 approximate)

求解调用同目录 hpipm_qp_solver.solve_ocp_qcqp (不复制).

本模块与 construct_ocp_qcqp.py / hpipm_qp_solver.py 同在 algorithms/RSS_proposed/
(算法实现归位算法包; pipeline/ 只留轨迹/仿真/评估主干)。

输入参数 params 为 pipeline.params.AlgorithmParams (由仿真器组装传入,
替代 MATLAB 版内部 config() 的角色), 包含:
    K, dt, vehicle(wheel_pos/vimax/phidotmax), weights(w_pos/w_psi/w_control/rho),
    max_iter (SCP 外层迭代数)
"""

import os
import sys
import numpy as np

from pipeline.params import AlgorithmParams
from pipeline.dynamics import rotation_matrix, wheel_matrices

# 同目录模块自举 (本目录由加载方或自身加入 sys.path)
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)
from construct_ocp_qcqp import construct_ocp_qcqp_from_rss

# hpipm_qp_solver 惰性导入缓存 (模块 import 时不触发 DLL 加载)
_solver_module = None


def _get_solver_module():
    """导入同目录 hpipm_qp_solver."""
    global _solver_module
    if _solver_module is None:
        import hpipm_qp_solver
        _solver_module = hpipm_qp_solver
    return _solver_module


def control_rss_ocpqcqp(path: np.ndarray, step: int, state_dot: np.ndarray,
                        state: np.ndarray, params: AlgorithmParams,
                        verbose: bool = True):
    """单个 MPC 步的 RSS 控制求解 (OCP QCQP + SCP).

    输入:
        path      : (3, N) 参考轨迹 [x; y; theta]
        step      : 当前 MPC 步号 (1-based, 与 MATLAB 一致)
        state_dot : (3,) 世界系速度 xdot_k
        state     : (3,) 世界系位姿 [x, y, psi]
        params    : AlgorithmParams
        verbose   : 是否打印 per-iteration 日志

    返回 (u, new_state_dot, velocity, diagnostics):
        u              : (3, K) 控制增量序列 (3 次 outer 后的 incumbent)
        new_state_dot  : (3,) 下一时刻世界系速度 = R(psi)*(state_dot + u_1)
        velocity       : (3,) 下一时刻车体系速度 = v0 + u_1
        diagnostics    : dict, 与 MATLAB diagnostics struct 对应
    """
    path = np.asarray(path, dtype=np.float64)
    state = np.asarray(state, dtype=np.float64).flatten()
    state_dot = np.asarray(state_dot, dtype=np.float64).flatten()

    # ================= 参数提取 =================
    K = int(params.K)
    current_xy = state[0:2].copy()
    psi0 = float(state[2])
    v0 = state_dot.copy()   # MATLAB: v0 = state_dot (闭环中实际传入车体系速度 nu_0)
    rho = float(params.weights.rho)

    # ================= 迭代 Setup =================
    max_iter = int(params.max_iter)   # 论文 IV-B: 固定 3 次外层迭代
    u_hat = np.zeros((3, K))          # u^(0) = 0 (static init, 与 Dense QCQP 一致)

    # 诊断结构体 (记录每次 outer 迭代)
    diagnostics = {
        'iterations': {
            'status': [''] * max_iter,
            'optval': [float('nan')] * max_iter,
            'solve_time': [float('nan')] * max_iter,
            'solver_name': [''] * max_iter,
            'hpipm_iters': [0] * max_iter,
            'nq': [0] * max_iter,
            'max_viol': [float('nan')] * max_iter,
            'u_diff_to_prev': [float('nan')] * max_iter,
            'qk_wheel_viol': [float('nan')] * max_iter,
            'qk_cone_viol': [float('nan')] * max_iter,
            'orig_wheel_viol': [float('nan')] * max_iter,
            'orig_cone_viol': [float('nan')] * max_iter,
        },
        'step': int(step),
        'max_iter': max_iter,
        'step_failed': False,
        'solver_call_count': 0,
    }

    total_solve_time = 0.0

    # ================= Outer 循环 (严格 3 次, 无 inner loop) =================
    for outer in range(max_iter):
        u_anchor = u_hat.copy()   # RSS 凸化锚点

        # 诊断默认值 (失败时保留)
        outer_status = 'Failed'
        outer_optval = float('nan')
        outer_solve_time = float('nan')
        outer_solver_name = 'HPIPM-OCPQCQP'
        outer_hpipm_iters = 0
        outer_nq = 0
        outer_max_viol = float('nan')
        outer_u_diff_val = float('nan')
        solver_calls_this = 0

        try:
            # ===== 论文 Alg.1 line 4: 构造精确凸子问题 Q_K(u_anchor) =====
            ocp = construct_ocp_qcqp_from_rss(path, step, v0, state, u_anchor, params)
            outer_nq = int(np.sum(ocp['nq_per_stage']))

            # ===== 论文 Alg.1 line 5: 求解 u^(m+1) = S(u^(m)) =====
            solver = _get_solver_module()
            result = solver.solve_ocp_qcqp(
                ocp['A'], ocp['B'], ocp['b_stack'],
                ocp['Q_stack'], ocp['S_eff'], ocp['R_eff'],
                ocp['q_stack'], ocp['r_stack'],
                ocp['nx'], ocp['nu'], ocp['nq'], ocp['nbx'], ocp['nq_per_stage'],
                ocp['Qq_stack'], ocp['Sq_stack'], ocp['Rq_stack'],
                ocp['qq_stack'], ocp['rq_stack'], ocp['uq_stack'],
                ocp['x0'], ocp['idxbx'], float(ocp['const']),
                False,                               # verbose
            )

            # ===== 提取结果 =====
            x = np.asarray(result['x'], dtype=np.float64).flatten()
            status_code = int(result['status'])
            outer_optval = float(result['obj_value_full'])
            outer_solve_time = float(result['solve_time'])
            outer_hpipm_iters = int(result['iters'])
            solver_calls_this = int(result['solver_call_count'])
            outer_solver_name = 'HPIPM-OCPQCQP'

            # x_out = [u(:); nu(:)] (与 Dense QCQP 格式一致)
            # MATLAB: u_sol = reshape(x(1:3K), 3, K) 列主序
            u_sol = x[0:3 * K].reshape((3, K), order='F')

        except Exception as exc:   # noqa: BLE001 - 与 MATLAB catch ME 一致
            u_sol = np.zeros((3, K))
            status_code = -1
            outer_optval = float('nan')
            outer_solve_time = float('nan')
            outer_hpipm_iters = 0
            solver_calls_this = 0
            outer_solver_name = 'HPIPM-Error'
            if verbose:
                print(f'第{step}步 outer={outer + 1} - Python 调用异常: {exc}')

        # 累计真实 solver 调用次数
        diagnostics['solver_call_count'] += solver_calls_this
        total_solve_time += 0.0 if np.isnan(outer_solve_time) else outer_solve_time

        if status_code == 0:
            outer_status = 'Solved'
        else:
            outer_status = f'Failed({status_code})'

        if verbose:
            print(f'第{step}步第{outer + 1}次迭代 - HPIPM-OCPQCQP内部求解时间：'
                  f'{outer_solve_time:.6f}秒 (status={status_code})')
            print(f'最优代价: {outer_optval} | 求解状态: {outer_status}')

        # ===== 严格接受条件: status==0 且解有限 =====
        # status=1 (MAX_ITER) 不视为成功
        if status_code == 0 and np.all(np.isfinite(u_sol)):
            outer_u_diff = float(np.max(np.abs(u_sol - u_hat)))
            outer_u_diff_val = outer_u_diff

            # 精确约束检查 (固定 Q_K(u_anchor))
            max_viol, wheel_viol, cone_viol = check_qk_violation(
                u_sol, u_anchor, v0, params)
            outer_max_viol = max_viol

            # 原始物理约束违反量 (双线性约束, 与 Q_K 凸化约束不同)
            orig_wheel_viol, orig_cone_viol = compute_original_violation(
                u_sol, v0, params)
            diagnostics['iterations']['orig_wheel_viol'][outer] = orig_wheel_viol
            diagnostics['iterations']['orig_cone_viol'][outer] = orig_cone_viol
            diagnostics['iterations']['qk_wheel_viol'][outer] = wheel_viol
            diagnostics['iterations']['qk_cone_viol'][outer] = cone_viol

            # 更新 u_hat (与 Dense QCQP 一致: status==0 即更新)
            u_hat = u_sol.copy()

        # 记录 outer 诊断
        diagnostics['iterations']['status'][outer] = outer_status
        diagnostics['iterations']['optval'][outer] = outer_optval
        diagnostics['iterations']['solve_time'][outer] = outer_solve_time
        diagnostics['iterations']['solver_name'][outer] = outer_solver_name
        diagnostics['iterations']['hpipm_iters'][outer] = outer_hpipm_iters
        diagnostics['iterations']['nq'][outer] = outer_nq
        if not np.isnan(outer_max_viol):
            diagnostics['iterations']['max_viol'][outer] = outer_max_viol
        if not np.isnan(outer_u_diff_val):
            diagnostics['iterations']['u_diff_to_prev'][outer] = outer_u_diff_val

    # ================= 输出: 3 次 outer 全部成功才标记 strict =================
    u = u_hat
    all_outer_solved = all(
        s == 'Solved' for s in diagnostics['iterations']['status'])
    if all_outer_solved:
        diagnostics['incumbent_type'] = 'strict'
        diagnostics['step_approximate'] = False
        diagnostics['has_strict_incumbent'] = True
        diagnostics['step_failed'] = False
    else:
        diagnostics['incumbent_type'] = 'none'
        diagnostics['step_approximate'] = False
        diagnostics['has_strict_incumbent'] = False
        diagnostics['step_failed'] = True

    # ================= 论文 Alg.1 line 11: 输出 nu_1 = nu_0 + u_1 =================
    new_state_dot = rotation_matrix(psi0) @ (state_dot + 1.00 * u[:, 0])
    velocity = v0 + u[:, 0]

    # ================= [DIAGNOSTIC] 原始约束违反量检查 (per-step) =================
    orig_wheel_viol_final, orig_cone_viol_final = compute_original_violation(
        u, v0, params)
    diagnostics['orig_wheel_viol_final'] = orig_wheel_viol_final
    diagnostics['orig_cone_viol_final'] = orig_cone_viol_final
    if verbose and (step % 25 == 0 or step == 1):
        print(f'[DIAG step={step}] wheel_viol_cur={orig_wheel_viol_final:.6e}, '
              f'cone_viol_cur={orig_cone_viol_final:.6e}')

    # 汇总诊断
    diagnostics['total_solve_time'] = total_solve_time

    return u, new_state_dot, velocity, diagnostics


def check_qk_violation(u_sol: np.ndarray, u_anchor: np.ndarray,
                       v0: np.ndarray, params: AlgorithmParams):
    """检查 u_sol 在固定 Q_K(u_anchor) 下的精确约束违反量.

    约束定义 (固定 Q_K(u_anchor)):
        1. 轮速 SOC (论文 20b): ||H_n * nu_k||^2 - vimax^2 <= 0
        2. 转向锥凸化 (论文 15-16): C = A - B(u_anchor) - L(u, u_anchor) <= 0

    返回 (max_viol, wheel_viol, cone_viol).
    """
    u_sol = np.asarray(u_sol, dtype=np.float64)
    u_anchor = np.asarray(u_anchor, dtype=np.float64)
    v0 = np.asarray(v0, dtype=np.float64).flatten()

    K = int(params.K)
    dt = float(params.dt)
    phidotmax = float(params.vehicle.phidotmax)
    vimax = float(params.vehicle.vimax)
    wheel_pos = np.asarray(params.vehicle.wheel_pos, dtype=np.float64)
    num_wheels = wheel_pos.shape[0]

    Hn = wheel_matrices(wheel_pos)

    delta_theta = dt * phidotmax
    # 论文 (12): R1 = R(pi/2 - delta_theta), R2 = R1^T
    R1 = np.array([
        [np.sin(delta_theta), -np.cos(delta_theta)],
        [np.cos(delta_theta), np.sin(delta_theta)],
    ], dtype=np.float64)
    R2 = np.array([
        [np.sin(delta_theta), np.cos(delta_theta)],
        [-np.cos(delta_theta), np.sin(delta_theta)],
    ], dtype=np.float64)

    # nu_sol 序列 (from u_sol)
    nu_sol = np.zeros((3, K + 1))
    nu_sol[:, 0] = v0
    for k in range(K):
        nu_sol[:, k + 1] = nu_sol[:, k] + u_sol[:, k]

    # nu_hat_anchor 序列 (from u_anchor)
    nu_hat_anchor = np.zeros((3, K + 1))
    nu_hat_anchor[:, 0] = v0
    for k in range(K):
        nu_hat_anchor[:, k + 1] = nu_hat_anchor[:, k] + u_anchor[:, k]

    # ---- 1. 轮速 SOC 约束 (k=1..K) ----
    wheel_viol = 0.0
    for k in range(1, K + 1):
        nu_k = nu_sol[:, k]
        for n in range(num_wheels):
            Mn = Hn[n].T @ Hn[n]
            val = float(nu_k @ Mn @ nu_k - vimax ** 2)   # ||H*nu_k||^2 - vimax^2
            if val > wheel_viol:
                wheel_viol = val

    # ---- 2. 转向锥凸化约束 C = A - B - L <= 0 (k=1..K) ----
    cone_viol = 0.0
    for k in range(1, K + 1):
        v_km1 = nu_sol[:, k - 1]           # nu_{k-1} (from u_sol)
        u_k = u_sol[:, k - 1]              # u_k (from u_sol)
        v_hat_km1 = nu_hat_anchor[:, k - 1]  # nu_hat_{k-1} (from u_anchor)
        u_hat_k = u_anchor[:, k - 1]       # u_hat_k (from u_anchor)

        for n in range(num_wheels):
            Hi = Hn[n]
            Mn = Hi.T @ Hi

            for gg in range(2):
                Rg = R1 if gg == 0 else R2
                T = (np.eye(2) + Rg) @ Hi    # (2,3)
                U_mat = Rg @ Hi              # (2,3)

                # B 项 (常数, 固定 u_anchor)
                ell_anchor = T @ v_hat_km1 + U_mat @ u_hat_k
                B_const = 0.5 * float(ell_anchor @ ell_anchor)

                # A 项 (凸二次, 在 u_sol 处求值)
                A_val = (float(v_km1 @ Mn @ v_km1)
                         + float(v_km1 @ Mn @ u_k)
                         + 0.5 * float(u_k @ Mn @ u_k))

                # L 项 (线性, 锚点 u_anchor, 求值点 u_sol)
                L_val = float(ell_anchor @ (T @ (v_km1 - v_hat_km1)
                                            + U_mat @ (u_k - u_hat_k)))

                # C = A - B - L
                C_val = A_val - B_const - L_val
                if C_val > cone_viol:
                    cone_viol = C_val

    max_viol = max(wheel_viol, cone_viol)
    return max_viol, wheel_viol, cone_viol


def compute_original_violation(u: np.ndarray, v0: np.ndarray,
                               params: AlgorithmParams):
    """原始物理约束违反量 (双线性, 非 Q_K 凸化) — compute_original_violation 等价.

    轮速: ||H_n * nu_k|| <= vimax
    转向锥: a'Mb >= 0 (a = nu_{k-1}, b = nu_k)

    返回 (wheel_viol, cone_viol).
    """
    u = np.asarray(u, dtype=np.float64)
    v0 = np.asarray(v0, dtype=np.float64).flatten()

    K = int(params.K)
    dt = float(params.dt)
    phidotmax = float(params.vehicle.phidotmax)
    vimax = float(params.vehicle.vimax)
    wheel_pos = np.asarray(params.vehicle.wheel_pos, dtype=np.float64)
    num_wheels = wheel_pos.shape[0]

    Hn = wheel_matrices(wheel_pos)
    delta_th = dt * phidotmax
    R1d = np.array([
        [np.sin(delta_th), -np.cos(delta_th)],
        [np.cos(delta_th), np.sin(delta_th)],
    ], dtype=np.float64)
    R2d = R1d.T

    nu = np.zeros((3, K + 1))
    nu[:, 0] = v0
    for kk in range(K):
        nu[:, kk + 1] = nu[:, kk] + u[:, kk]

    # 轮速约束 (k=1..K)
    wheel_viol = 0.0
    for kk in range(1, K + 1):
        nkk = nu[:, kk]
        for n in range(num_wheels):
            vn = float(np.linalg.norm(Hn[n] @ nkk, 2))
            viol = vn - vimax
            if viol > wheel_viol:
                wheel_viol = viol

    # 转向锥 (k=1..K)
    cone_viol = 0.0
    for kk in range(1, K + 1):
        a = nu[:, kk - 1]
        b = nu[:, kk]
        for n in range(num_wheels):
            Hn1 = Hn[n]
            for gg in range(2):
                Rgg = R1d if gg == 0 else R2d
                Mgg = Hn1.T @ Rgg @ Hn1
                term = float(a @ Mgg @ b)
                viol = -term
                if viol > cone_viol:
                    cone_viol = viol

    return wheel_viol, cone_viol
