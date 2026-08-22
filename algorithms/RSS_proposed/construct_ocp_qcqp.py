"""误差状态 OCP QCQP 构造 (HPIPM ocp_qcqp, 精确凸二次约束).

OCP QCQP 矩阵构造的唯一实现 (原 MATLAB 版 construct_ocp_qcqp_from_rss.m 已删除)
(L1 fixture 测试保证 stacks 逐位一致). 本模块与 control_rss_ocpqcqp.py /
hpipm_qp_solver.py 同在 algorithms/RSS_proposed/ (算法实现归位算法包;
pipeline/ 只留轨迹/仿真/评估主干), 求解调用同目录
hpipm_qp_solver.solve_ocp_qcqp (不复制).

状态 x_n = [e_n; v_n] in R^6, n=0..K (K+1 stages)
    e_n = xi_n - xi_n^ref (跟踪误差, 3 dim)
    v_n = nu_n (车体系速度, 3 dim)
控制 u_n in R^3, n=0..K-1 (u_n = v_{n+1} - v_n)

二次约束 stage 分布 (num_wheels=N):
    stage 0:    2N steering (k=1)
    stage 1..K-1: N wheel (k=s) + 2N steering (k=s+1)
    stage K:    N wheel (terminal)
    总数 3NK (K=6, N=4 -> 72)

×2 约定 (与 MATLAB 版一致):
    Q, R, Qq, Rq 构造时×2 (HPIPM 0.5 前缀还原); S, q, r, Sq, qq, rq 不×2;
    uq, const 不×2.
"""

import numpy as np

from pipeline.params import AlgorithmParams
from pipeline.dynamics import wheel_matrices


def construct_ocp_qcqp_from_rss(path: np.ndarray, step: int, v0: np.ndarray,
                                state: np.ndarray, u_anchor: np.ndarray,
                                params: AlgorithmParams) -> dict:
    """构造 OCP QCQP 全部数据 (OCP QCQP + SCP).

    输入:
        path      : (3, N) 参考轨迹 [x; y; theta]
        step      : 当前 MPC 步号 (1-based, 与 MATLAB 一致)
        v0        : (3,) 当前车体系速度 nu_0
        state     : (3,) 世界系位姿 [x, y, psi]
        u_anchor  : (3, K) RSS 凸化锚点
        params    : AlgorithmParams

    返回 dict (键与 MATLAB ocp struct 对应, cell 数组改为堆叠矩阵):
        A (6,6), B (6,3), b_stack (6,K), r_stack (3,K),
        Q_stack (6, 6*(K+1)), S_eff (3,6), R_eff (3,3), q_stack (6, K+1),
        const, K, N_stages, nx, nu, nbx, nq, nq_per_stage,
        Qq_stack (6,6,3NK), Sq_stack (3,6,3NK), Rq_stack (3,3,3NK),
        qq_stack (6,3NK), rq_stack (3,3NK), uq_stack (3NK,),
        x0 (6,), idxbx (6,), metadata
    """
    path = np.asarray(path, dtype=np.float64)
    v0 = np.asarray(v0, dtype=np.float64).flatten()
    state = np.asarray(state, dtype=np.float64).flatten()
    u_anchor = np.asarray(u_anchor, dtype=np.float64)

    # ===== 参数提取 (论文 IV-A, 与 MATLAB 版一致) =====
    K = int(params.K)
    dt = float(params.dt)
    phidotmax = float(params.vehicle.phidotmax)
    vimax = float(params.vehicle.vimax)
    wheel_pos = np.asarray(params.vehicle.wheel_pos, dtype=np.float64)
    num_wheels = wheel_pos.shape[0]

    w_pos = float(params.weights.w_pos)
    w_psi = float(params.weights.w_psi)
    w_control = float(params.weights.w_control)
    rho = float(params.weights.rho)

    current_xy = state[0:2].copy()
    psi0 = float(state[2])
    # MATLAB: R_psi0 = 2x2 旋转矩阵 (论文 (1) 的平面部分)
    _c, _s = np.cos(psi0), np.sin(psi0)
    R_psi0 = np.array([[_c, -_s], [_s, _c]], dtype=np.float64)

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

    # ===== 1. 动力学矩阵 (LTI) =====
    # MATLAB: A = [eye(3), [R_psi0*dt, zeros(2,1); zeros(1,2), dt]; zeros(3,3), eye(3)]
    A = np.zeros((6, 6))
    A[0:3, 0:3] = np.eye(3)
    A[0:2, 3:5] = R_psi0 * dt
    A[2, 5] = dt
    A[3:6, 3:6] = np.eye(3)
    B = np.zeros((6, 3))
    B[3:6, 0:3] = np.eye(3)

    # ===== 2. dim 设置 =====
    N_stages = K + 1
    nx = 6
    nu = 3

    nq_per_stage = np.zeros(N_stages, dtype=np.int32)
    nq_per_stage[0] = 2 * num_wheels                 # stage 0: steering
    for s in range(1, K):                            # stage 1..K-1
        nq_per_stage[s] = num_wheels + 2 * num_wheels
    nq_per_stage[K] = num_wheels                     # stage K: terminal wheel
    total_nq = int(nq_per_stage.sum())

    # ===== 3. 时变量: b_stack, r_stack, const =====
    const = 0.0
    b_stack = np.zeros((6, K))
    r_stack = np.zeros((3, K))
    for n in range(K):
        ref_idx_k = min(num_path_pts, step + n)          # xi_k^ref
        ref_idx_kp1 = min(num_path_pts, step + n + 1)    # xi_{k+1}^ref
        ref_k = path[:, ref_idx_k - 1]
        ref_kp1 = path[:, ref_idx_kp1 - 1]
        b_xi = ref_k - ref_kp1
        b_stack[:, n] = np.concatenate([b_xi, np.zeros(3)])
        r_stack[:, n] = -2.0 * rho * u_anchor[:, n]
        const += rho * float(u_anchor[:, n] @ u_anchor[:, n])

    # ===== 4. 初始状态 (box bounds 固定) =====
    ref_idx_0 = min(num_path_pts, step)
    ref_0 = path[:, ref_idx_0 - 1]
    e0 = np.array([
        current_xy[0] - ref_0[0],
        current_xy[1] - ref_0[1],
        psi0 - ref_0[2],
    ], dtype=np.float64)
    x0 = np.concatenate([e0, v0])

    # ===== 5. 逐 stage 代价矩阵 =====
    Q_pos_psi = 2.0 * np.diag([w_pos, w_pos, w_psi, 0.0, 0.0, 0.0])  # stage 2..K
    Q_psi_only = 2.0 * np.diag([0.0, 0.0, w_psi, 0.0, 0.0, 0.0])     # stage 1
    Q_zero = np.zeros((6, 6))                                         # stage 0

    Q_stack = np.zeros((nx, nx * N_stages))
    Q_stack[:, 0:6] = Q_zero
    Q_stack[:, 6:12] = Q_psi_only
    for s in range(3, N_stages + 1):
        Q_stack[:, (s - 1) * nx:s * nx] = Q_pos_psi

    R_eff = 2.0 * (w_control * np.eye(nu) + rho * np.eye(nu))
    S_eff = np.zeros((nu, nx))
    q_stack = np.zeros((nx, N_stages))

    # ===== 6. 锚点速度序列 v_hat (转向锥凸化 B/L 项用) =====
    nu_hat_anchor = np.zeros((3, K + 1))
    nu_hat_anchor[:, 0] = v0
    for k in range(1, K + 1):
        nu_hat_anchor[:, k] = nu_hat_anchor[:, k - 1] + u_anchor[:, k - 1]

    # ===== 7. 构造全部精确二次约束 =====
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

    # ---- 7.1 Stage 0: 2N steering (k=1, involves x_0 and u_0) ----
    s = 0
    k = 1
    v_hat = nu_hat_anchor[:, k - 1]      # v_hat_0 = v0
    u_hat_k = u_anchor[:, k - 1]         # û_1
    for n in range(num_wheels):
        H_n = Hn[n]
        M = H_n.T @ H_n                  # (3,3) 对称
        for gg in range(2):
            Rg = R1 if gg == 0 else R2
            rot_name = 'R1' if gg == 0 else 'R2'
            T = (np.eye(2) + Rg) @ H_n   # (2,3)
            U = Rg @ H_n                 # (2,3)
            ell = T @ v_hat + U @ u_hat_k

            Qq = np.zeros((nx, nx))
            Qq[3:6, 3:6] = 2.0 * M
            Qq_stack[:, :, idx] = Qq
            Sq = np.zeros((nu, nx))
            Sq[:, 3:6] = M
            Sq_stack[:, :, idx] = Sq
            Rq_stack[:, :, idx] = M
            qq = np.zeros(nx)
            qq[3:6] = -(T.T @ ell)
            qq_stack[:, idx] = qq
            rq_stack[:, idx] = -(U.T @ ell)
            uq_stack[idx] = -0.5 * float(ell @ ell)

            _meta_append(metadata, 'steering', k, n + 1, rot_name, s, idx,
                         idx - int(nq_per_stage[:s].sum()))
            idx += 1

    # ---- 7.2 Stage 1..K-1: N wheel (k=s) + 2N steering (k=s+1) ----
    for s in range(1, K):
        # N wheel 约束 (k=s, on x_s)
        k_wheel = s
        for n in range(num_wheels):
            H_n = Hn[n]
            M = H_n.T @ H_n
            Qq = np.zeros((nx, nx))
            Qq[3:6, 3:6] = 2.0 * M
            Qq_stack[:, :, idx] = Qq
            uq_stack[idx] = vimax ** 2
            _meta_append(metadata, 'wheel', k_wheel, n + 1, '', s, idx,
                         idx - int(nq_per_stage[:s].sum()))
            idx += 1

        # 2N steering 约束 (k=s+1, involves x_s and u_s)
        k_steer = s + 1
        v_hat = nu_hat_anchor[:, k_steer - 1]    # v_hat_s
        u_hat_k = u_anchor[:, k_steer - 1]       # û_{s+1}
        for n in range(num_wheels):
            H_n = Hn[n]
            M = H_n.T @ H_n
            for gg in range(2):
                Rg = R1 if gg == 0 else R2
                rot_name = 'R1' if gg == 0 else 'R2'
                T = (np.eye(2) + Rg) @ H_n
                U = Rg @ H_n
                ell = T @ v_hat + U @ u_hat_k

                Qq = np.zeros((nx, nx))
                Qq[3:6, 3:6] = 2.0 * M
                Qq_stack[:, :, idx] = Qq
                Sq = np.zeros((nu, nx))
                Sq[:, 3:6] = M
                Sq_stack[:, :, idx] = Sq
                Rq_stack[:, :, idx] = M
                qq = np.zeros(nx)
                qq[3:6] = -(T.T @ ell)
                qq_stack[:, idx] = qq
                rq_stack[:, idx] = -(U.T @ ell)
                uq_stack[idx] = -0.5 * float(ell @ ell)

                _meta_append(metadata, 'steering', k_steer, n + 1, rot_name, s, idx,
                             idx - int(nq_per_stage[:s].sum()))
                idx += 1

    # ---- 7.3 Stage K: N wheel (k=K, terminal) ----
    s = K
    k_wheel = K
    for n in range(num_wheels):
        H_n = Hn[n]
        M = H_n.T @ H_n
        Qq = np.zeros((nx, nx))
        Qq[3:6, 3:6] = 2.0 * M
        Qq_stack[:, :, idx] = Qq
        # 终端 stage nu=0: Sq/Rq/rq 不设置 (保持零, 求解器端跳过)
        uq_stack[idx] = vimax ** 2
        _meta_append(metadata, 'wheel', k_wheel, n + 1, '', s, idx,
                     idx - int(nq_per_stage[:s].sum()))
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
        'idxbx': np.arange(6, dtype=np.int32),
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
