"""velocity-only OCP QCQP 实验 (HPIPM ocp_qcqp, nx=3 速度态).

实验性降维构造 (与原版 construct_ocp_qcqp.py 并存, 不替代):
    状态只剩车体系速度 v_n ∈ R^3 (nx=3, 原版为 [e; v] ∈ R^6),
    控制仍为 u_n ∈ R^3 (u_n = v_{n+1} - v_n),
    动力学退化为 v_{n+1} = v_n + u_n (A=I, B=I, b=0).

与原版的关系:
    - 72 条二次约束 (轮速 SOC + 转向锥凸化) 逐字不变 —— 它们本来就只
      作用在 (v, u) 上, 位置误差 e 从不出现, 打包维度 6→3
    - 代价不再含位置/航向误差项, 改为速度跟踪 ‖v_n - v_ref_n‖^2_W:
        v_ref_n = 前馈 (路径差分转车体系) + P 反馈 (当前位姿误差)
      反馈在 QCQP 之外 (外环 P, 增益 VEL_K_POS/VEL_K_PSI), MPC 内无
      位置状态, 每步只重置 v_0
    - 数学上与原版不等价: 原版位置代价在 u 坐标下稠密, OCP stage-wise
      结构在 nx=3 下无法表达 (跨 stage 块依赖 i), 故闭环结果不重合
      golden 基准 (RMSE=0.036793 / J_total=13.3838), 仅作复杂度/性能对比

调参结论 (seed0_K6 euler, VEL_W_SCALE=0.03 / VEL_K_POS=VEL_K_PSI=20):
    RMSE=0.033529 (基线 0.036793), J_total=9.9711 (基线 13.3838,
    J_psi 7.36→1.94 / J_u 1.96→4.66, 以控制能量换航向精度),
    validSteps=100/100, 约束违反 ~0 (maxWheelViol=1.3e-10);
    seed1 复核: RMSE=0.018443 (基线 0.021528), 182/182.
    medianSolveTime≈0.0019s vs 基线 0.0021s —— 72 条二次约束主导
    IPM 耗时, nx 减半 (primal 60→39) 只带来 ~8% 提速 (墙钟噪声级).

求解仍走同目录 hpipm_qp_solver.solve_ocp_qcqp (HPIPM OCP QCQP IPM,
同一 DLL 与设置); 求解器接口维度无关 (nx_arr/nu_arr 驱动, 零改动).
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
from control_rss_ocpqcqp import check_qk_violation, compute_original_violation

# P 反馈增益 (1/s): 位姿误差 → 参考速度修正 (QCQP 之外的外环)
# v_ref_xy += k_pos * R(psi)^T (p_ref - p_cur);  omega_ref += k_psi * e_psi
# seed0 扫描 (k=5/10/15/20/25): RMSE 0.064/0.051/0.048/0.034/0.034,
# J 25.1/24.8/14.8/9.97/11.0 — k=20 最优 (RMSE 0.0335 vs 基线 0.0368)
VEL_K_POS = 20.0  # 位置 → 速度 (时间常数 1/k = 0.05 s)
VEL_K_PSI = 20.0  # 航向 → 角速度

# 速度跟踪权重缩放 (相对论文 w_pos/w_psi). 原版刚度作用在位置误差 (量级
# ~1e-1) 上, 消掉位置后若把同样的权重直接压在速度 (量级 ~1e0) 上, QP 线性
# 项比约束曲率大 ~2 个量级, SCP 锚点激进时触发 IPM MIN_STEP (step 49 复
# 现). 缩放后 W=diag(30s, 30s, s), 与控制正则化 R=0.3 平衡.
# seed0 扫描 (s=1/0.1/0.03/0.01 @k=20): 0.1 仍 step 55 失败, 0.03 最优.
VEL_W_SCALE = 0.03

# 环境变量覆盖 (参数扫描实验用; 正式跑不设即用上面的默认值)
if 'VEL_W_SCALE' in os.environ:
    VEL_W_SCALE = float(os.environ['VEL_W_SCALE'])
if 'VEL_K_POS' in os.environ:
    VEL_K_POS = float(os.environ['VEL_K_POS'])
if 'VEL_K_PSI' in os.environ:
    VEL_K_PSI = float(os.environ['VEL_K_PSI'])


# ===== 可选: 增益调度 + 反馈 preview 衰减 (默认关闭, 保持已发布行为) =====
# 贴向基准的两个结构差异对策 (均环境变量门控, 不设即完全关闭):
#   VEL_K_SCHED  : P 增益按 MPC 步号调度, '40:15,20' = 第1-15步 k=40, 之后 20
#   VEL_W_SCHED  : 速度跟踪权重缩放按步号调度, '0.05:15,0.01' 同上格式
#   VEL_FB_PREVIEW='1' : v_ref 的 P 反馈项按闭环误差衰减律 (1-k*dt)^n 逐
#       stage 衰减 —— 模拟基准 MPC 位置代价的 preview 效应 (基准惩罚
#       horizon 内每一步的预测位置误差, 修正自然前移; 常数反馈则把修正
#       平摊到整个 horizon, 起步收敛形状与基准不同)
def _parse_sched(spec):
    """'40:15,20' -> [(15, 40.0), (inf, 20.0)]; 空 -> None (不调度)."""
    if not spec:
        return None
    out = []
    for part in spec.split(','):
        seg = part.split(':')
        if len(seg) == 2:
            out.append((float(seg[0]), float(seg[1])))
        else:
            out.append((float('inf'), float(seg[0])))
    return out


_K_SCHED = _parse_sched(os.environ.get('VEL_K_SCHED'))
_W_SCHED = _parse_sched(os.environ.get('VEL_W_SCHED'))
_FB_PREVIEW = os.environ.get('VEL_FB_PREVIEW', '0') == '1'


def _sched_at(sched, default, step):
    if sched is None:
        return default
    for upto, val in sched:
        if step <= upto:
            return val
    return sched[-1][1]

# hpipm_qp_solver 惰性导入缓存 (模块 import 时不触发 DLL 加载)
_solver_module = None


def _get_solver_module():
    """导入同目录 hpipm_qp_solver."""
    global _solver_module
    if _solver_module is None:
        import hpipm_qp_solver
        _solver_module = hpipm_qp_solver
    return _solver_module


def _wrap_angle(a: float) -> float:
    return (a + np.pi) % (2.0 * np.pi) - np.pi


def construct_ocp_qcqp_vel(path: np.ndarray, step: int, v0: np.ndarray,
                           state: np.ndarray, u_anchor: np.ndarray,
                           params: AlgorithmParams) -> dict:
    """构造 velocity-only OCP QCQP 全部数据 (实验性, nx=3).

    输入 (与 construct_ocp_qcqp_from_rss 相同):
        path      : (3, N) 参考轨迹 [x; y; theta]
        step      : 当前 MPC 步号 (1-based)
        v0        : (3,) 当前车体系速度 nu_0
        state     : (3,) 世界系位姿 [x, y, psi] (供 v_ref 反馈项使用)
        u_anchor  : (3, K) RSS 凸化锚点
        params    : AlgorithmParams

    返回 dict (键与原版一致):
        A (3,3)=I, B (3,3)=I, b_stack (3,K)=0, r_stack (3,K),
        Q_stack (3, 3*(K+1)), S_eff (3,3)=0, R_eff (3,3), q_stack (3, K+1),
        const, K, N_stages, nx, nu, nbx, nq, nq_per_stage,
        Qq_stack (3,3,3NK), Sq_stack (3,3,3NK), Rq_stack (3,3,3NK),
        qq_stack (3,3NK), rq_stack (3,3NK), uq_stack (3NK,),
        x0 (3,)=v0, idxbx (3,), v_ref (3, K+1), metadata
    """
    path = np.asarray(path, dtype=np.float64)
    v0 = np.asarray(v0, dtype=np.float64).flatten()
    state = np.asarray(state, dtype=np.float64).flatten()
    u_anchor = np.asarray(u_anchor, dtype=np.float64)

    # ===== 参数提取 =====
    K = int(params.K)
    dt = float(params.dt)
    phidotmax = float(params.vehicle.phidotmax)
    vimax = float(params.vehicle.vimax)
    wheel_pos = np.asarray(params.vehicle.wheel_pos, dtype=np.float64)
    num_wheels = wheel_pos.shape[0]

    w_scale = _sched_at(_W_SCHED, VEL_W_SCALE, step)
    w_track = w_scale * float(params.weights.w_pos)    # xy 速度跟踪权重
    w_om = w_scale * float(params.weights.w_psi)       # 角速度跟踪权重
    w_control = float(params.weights.w_control)
    rho = float(params.weights.rho)

    current_xy = state[0:2].copy()
    psi0 = float(state[2])

    Hn = wheel_matrices(wheel_pos)           # 每轮 (2,3)
    delta_theta = dt * phidotmax             # 论文 (11)
    # 论文 (12): R1 = R(pi/2 - delta_theta), R2 = R1^T
    R1 = np.array([
        [np.sin(delta_theta), -np.cos(delta_theta)],
        [np.cos(delta_theta), np.sin(delta_theta)],
    ], dtype=np.float64)
    R2 = np.array([
        [np.sin(delta_theta), np.cos(delta_theta)],
        [-np.cos(delta_theta), np.sin(delta_theta)],
    ], dtype=np.float64)

    num_path_pts = path.shape[1]

    # ===== 1. 动力学 (LTI): v_{n+1} = v_n + u_n =====
    A = np.eye(3)
    B = np.eye(3)
    b_stack = np.zeros((3, K))

    # ===== 2. dim 设置 =====
    N_stages = K + 1
    nx = 3
    nu = 3

    nq_per_stage = np.zeros(N_stages, dtype=np.int32)
    nq_per_stage[0] = 2 * num_wheels                 # stage 0: steering
    for s in range(1, K):                            # stage 1..K-1
        nq_per_stage[s] = num_wheels + 2 * num_wheels
    nq_per_stage[K] = num_wheels                     # stage K: terminal wheel
    total_nq = int(nq_per_stage.sum())

    # ===== 3. 速度参考 v_ref = 前馈 + P 反馈 =====
    # 反馈 (整个 horizon 共享, 每步 MPC 重算): 当前位姿误差
    # VEL_FB_PREVIEW=1 时逐 stage 按闭环衰减律 (1-k*dt)^n 缩放 (preview 模拟)
    ref_idx_0 = min(num_path_pts, step)
    ref_0 = path[:, ref_idx_0 - 1]
    e_pos_world = ref_0[0:2] - current_xy            # (2,) 世界系位置误差
    e_psi = _wrap_angle(ref_0[2] - psi0)
    _c, _s = np.cos(psi0), np.sin(psi0)
    R_psi0_2 = np.array([[_c, -_s], [_s, _c]])       # R(psi0) 平面部分
    k_pos = _sched_at(_K_SCHED, VEL_K_POS, step)
    k_psi = _sched_at(_K_SCHED, VEL_K_PSI, step)
    fb_body_xy = k_pos * (R_psi0_2.T @ e_pos_world)   # 车体系
    fb_omega = k_psi * e_psi
    decay_xy = (1.0 - k_pos * dt) if _FB_PREVIEW else 1.0
    decay_om = (1.0 - k_psi * dt) if _FB_PREVIEW else 1.0

    # 前馈: 逐 stage 路径差分转车体系 (参考航向处)
    v_ref = np.zeros((3, N_stages))
    v_ref[:, 0] = v0                                 # stage 0 代价为零
    for n in range(K):
        ref_k = path[:, min(num_path_pts, step + n) - 1]
        ref_kp1 = path[:, min(num_path_pts, step + n + 1) - 1]
        dp = ref_kp1[0:2] - ref_k[0:2]
        dth = _wrap_angle(ref_kp1[2] - ref_k[2])
        ck, sk = np.cos(ref_k[2]), np.sin(ref_k[2])
        ff_body = np.array([[ck, sk], [-sk, ck]]) @ (dp / dt)  # R(theta)^T dp/dt
        v_ref[:, n + 1] = np.array([
            ff_body[0] + (decay_xy ** (n + 1)) * fb_body_xy[0],
            ff_body[1] + (decay_xy ** (n + 1)) * fb_body_xy[1],
            dth / dt + (decay_om ** (n + 1)) * fb_omega,
        ])

    # ===== 4. 代价矩阵 (速度跟踪 + 控制正则化) =====
    W = np.diag([w_track, w_track, w_om])
    Q_stack = np.zeros((nx, nx * N_stages))
    q_stack = np.zeros((nx, N_stages))
    const = 0.0
    for n in range(1, N_stages):                     # stage 0 无代价 (v_0 固定)
        Q_stack[:, n * nx:(n + 1) * nx] = 2.0 * W
        q_stack[:, n] = -2.0 * (W @ v_ref[:, n])
        const += float(v_ref[:, n] @ W @ v_ref[:, n])

    R_eff = 2.0 * (w_control * np.eye(nu) + rho * np.eye(nu))
    S_eff = np.zeros((nu, nx))
    r_stack = np.zeros((nu, K))
    for n in range(K):
        r_stack[:, n] = -2.0 * rho * u_anchor[:, n]
        const += rho * float(u_anchor[:, n] @ u_anchor[:, n])

    # ===== 5. 初始状态 (box bounds 固定 v0) =====
    x0 = v0.copy()

    # ===== 6. 锚点速度序列 v_hat (转向锥凸化 B/L 项用) =====
    nu_hat_anchor = np.zeros((3, K + 1))
    nu_hat_anchor[:, 0] = v0
    for k in range(1, K + 1):
        nu_hat_anchor[:, k] = nu_hat_anchor[:, k - 1] + u_anchor[:, k - 1]

    # ===== 7. 构造全部精确二次约束 (与原版逐字一致, 打包 nx=3) =====
    Qq_stack = np.zeros((nx, nx, total_nq))
    Sq_stack = np.zeros((nu, nx, total_nq))
    Rq_stack = np.zeros((nu, nu, total_nq))
    qq_stack = np.zeros((nx, total_nq))
    rq_stack = np.zeros((nu, total_nq))
    uq_stack = np.zeros(total_nq)

    metadata = {
        'kind': [], 'k': [], 'wheel': [], 'rotation': [],
        'stage': [], 'local_index': [],
    }

    idx = 0  # 0-based 全局约束索引

    # ---- 7.1 Stage 0: 2N steering (k=1, involves v_0 and u_0) ----
    v_hat = nu_hat_anchor[:, 0]             # v_hat_0 = v0
    u_hat_k = u_anchor[:, 0]                # û_1
    for n in range(num_wheels):
        H_n = Hn[n]
        M = H_n.T @ H_n                     # (3,3) 对称
        for gg in range(2):
            Rg = R1 if gg == 0 else R2
            rot_name = 'R1' if gg == 0 else 'R2'
            T = (np.eye(2) + Rg) @ H_n      # (2,3)
            U = Rg @ H_n                    # (2,3)
            ell = T @ v_hat + U @ u_hat_k

            Qq_stack[:, :, idx] = 2.0 * M
            Sq_stack[:, :, idx] = M
            Rq_stack[:, :, idx] = M
            qq_stack[:, idx] = -(T.T @ ell)
            rq_stack[:, idx] = -(U.T @ ell)
            uq_stack[idx] = -0.5 * float(ell @ ell)

            _meta_append(metadata, 'steering', 1, n + 1, rot_name, 0, idx,
                         idx - int(nq_per_stage[:0].sum()))
            idx += 1

    # ---- 7.2 Stage 1..K-1: N wheel (k=s) + 2N steering (k=s+1) ----
    for s in range(1, K):
        # N wheel 约束 (k=s, on v_s)
        for n in range(num_wheels):
            H_n = Hn[n]
            M = H_n.T @ H_n
            Qq_stack[:, :, idx] = 2.0 * M
            uq_stack[idx] = vimax ** 2
            _meta_append(metadata, 'wheel', s, n + 1, '', s, idx,
                         idx - int(nq_per_stage[:s].sum()))
            idx += 1

        # 2N steering 约束 (k=s+1, involves v_s and u_s)
        v_hat = nu_hat_anchor[:, s]
        u_hat_k = u_anchor[:, s]
        for n in range(num_wheels):
            H_n = Hn[n]
            M = H_n.T @ H_n
            for gg in range(2):
                Rg = R1 if gg == 0 else R2
                rot_name = 'R1' if gg == 0 else 'R2'
                T = (np.eye(2) + Rg) @ H_n
                U = Rg @ H_n
                ell = T @ v_hat + U @ u_hat_k

                Qq_stack[:, :, idx] = 2.0 * M
                Sq_stack[:, :, idx] = M
                Rq_stack[:, :, idx] = M
                qq_stack[:, idx] = -(T.T @ ell)
                rq_stack[:, idx] = -(U.T @ ell)
                uq_stack[idx] = -0.5 * float(ell @ ell)

                _meta_append(metadata, 'steering', s + 1, n + 1, rot_name, s, idx,
                             idx - int(nq_per_stage[:s].sum()))
                idx += 1

    # ---- 7.3 Stage K: N wheel (k=K, terminal) ----
    for n in range(num_wheels):
        H_n = Hn[n]
        M = H_n.T @ H_n
        Qq_stack[:, :, idx] = 2.0 * M
        # 终端 stage nu=0: Sq/Rq/rq 不设置 (保持零, 求解器端跳过)
        uq_stack[idx] = vimax ** 2
        _meta_append(metadata, 'wheel', K, n + 1, '', K, idx,
                     idx - int(nq_per_stage[:K].sum()))
        idx += 1

    assert idx == total_nq, f'约束总数不匹配: idx={idx}, total_nq={total_nq}'

    # ===== 8. 返回 =====
    return {
        'A': A,
        'B': B,
        'b_stack': b_stack,
        'r_stack': r_stack,
        'Q_stack': Q_stack,
        'S_eff': S_eff,
        'R_eff': R_eff,
        'q_stack': q_stack,
        'const': const,
        'K': K,
        'N_stages': N_stages,
        'nx': np.full(N_stages, nx, dtype=np.int32),
        'nu': np.array([nu] * K + [0], dtype=np.int32),
        'nbx': np.array([nx] + [0] * K, dtype=np.int32),
        'nq': nq_per_stage.copy(),
        'nq_per_stage': nq_per_stage.copy(),
        'Qq_stack': Qq_stack,
        'Sq_stack': Sq_stack,
        'Rq_stack': Rq_stack,
        'qq_stack': qq_stack,
        'rq_stack': rq_stack,
        'uq_stack': uq_stack,
        'x0': x0,
        'idxbx': np.arange(3, dtype=np.int32),
        'v_ref': v_ref,
        'metadata': metadata,
    }


def _meta_append(metadata: dict, kind: str, k: int, wheel: int,
                 rotation: str, stage: int, global_idx: int, local_idx: int):
    metadata['kind'].append(kind)
    metadata['k'].append(k)
    metadata['wheel'].append(wheel)
    metadata['rotation'].append(rotation)
    metadata['stage'].append(stage)
    metadata['local_index'].append(local_idx)


def control_rss_vel(path: np.ndarray, step: int, state_dot: np.ndarray,
                    state: np.ndarray, params: AlgorithmParams,
                    verbose: bool = True):
    """单个 MPC 步的 velocity-only RSS 控制求解 (OCP QCQP + SCP).

    签名/返回与 control_rss_ocpqcqp 完全一致 (仿真器无感切换):
        返回 (u, new_state_dot, velocity, diagnostics)
    SCP 结构也一致: max_iter 次 outer 迭代 (默认 3), 每次恰好 1 次
    solve_ocp_qcqp, status==0 且解有限才更新 incumbent.
    """
    path = np.asarray(path, dtype=np.float64)
    state = np.asarray(state, dtype=np.float64).flatten()
    state_dot = np.asarray(state_dot, dtype=np.float64).flatten()

    # ================= 参数提取 =================
    K = int(params.K)
    psi0 = float(state[2])
    v0 = state_dot.copy()   # 闭环中实际传入车体系速度 nu_0

    # 逐步 rho 调度 (CLI --rho), 与 control_rss_ocpqcqp 相同的保存/恢复模式
    rho_orig = float(params.weights.rho)
    if getattr(params, 'rho_schedule', None):
        rho_idx = min(int(step) - 1, len(params.rho_schedule) - 1)
        params.weights.rho = float(params.rho_schedule[rho_idx])
    rho = float(params.weights.rho)

    # ================= 迭代 Setup =================
    max_iter = int(params.max_iter)   # 默认 3 (与原版一致)
    u_hat = np.zeros((3, K))          # u^(0) = 0

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

    # ================= Outer 循环 (严格 max_iter 次, 无 inner loop) =================
    for outer in range(max_iter):
        u_anchor = u_hat.copy()   # RSS 凸化锚点

        outer_status = 'Failed'
        outer_optval = float('nan')
        outer_solve_time = float('nan')
        outer_solver_name = 'HPIPM-OCPQCQP-VEL'
        outer_hpipm_iters = 0
        outer_nq = 0
        outer_max_viol = float('nan')
        outer_u_diff_val = float('nan')
        solver_calls_this = 0

        try:
            # ===== 构造 velocity-only 凸子问题 Q_K(u_anchor) =====
            ocp = construct_ocp_qcqp_vel(path, step, v0, state, u_anchor, params)
            outer_nq = int(np.sum(ocp['nq_per_stage']))

            # ===== 求解 u^(m+1) = S(u^(m)) =====
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

            # x_out = [u(:); ...] (u 部分提取与原版一致, 列主序)
            u_sol = x[0:3 * K].reshape((3, K), order='F')

        except Exception as exc:   # noqa: BLE001 - 与原版 catch 一致
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
            print(f'第{step}步第{outer + 1}次迭代 - HPIPM-OCPQCQP-VEL内部求解时间：'
                  f'{outer_solve_time:.6f}秒 (status={status_code})')
            print(f'最优代价: {outer_optval} | 求解状态: {outer_status}')

        # ===== 严格接受条件: status==0 且解有限 =====
        if status_code == 0 and np.all(np.isfinite(u_sol)):
            outer_u_diff = float(np.max(np.abs(u_sol - u_hat)))
            outer_u_diff_val = outer_u_diff

            # 精确约束检查 (固定 Q_K(u_anchor); 复用原版函数, 只吃 u/v)
            max_viol, wheel_viol, cone_viol = check_qk_violation(
                u_sol, u_anchor, v0, params)
            outer_max_viol = max_viol

            orig_wheel_viol, orig_cone_viol = compute_original_violation(
                u_sol, v0, params)
            diagnostics['iterations']['orig_wheel_viol'][outer] = orig_wheel_viol
            diagnostics['iterations']['orig_cone_viol'][outer] = orig_cone_viol
            diagnostics['iterations']['qk_wheel_viol'][outer] = wheel_viol
            diagnostics['iterations']['qk_cone_viol'][outer] = cone_viol

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

    # ================= 输出 =================
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

    # nu_1 = nu_0 + u_1 (与原版一致)
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

    diagnostics['total_solve_time'] = total_solve_time
    diagnostics['rho'] = rho

    # 恢复传入前的 weights.rho
    params.weights.rho = rho_orig

    return u, new_state_dot, velocity, diagnostics
