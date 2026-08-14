"""
hpipm_qp_solver.py
HPIPM dense QCQP / OCP QCQP / OCP QP 求解器封装 (Python 接口)

论文: "hpipm: a high-performance quadratic programming framework for
       model predictive control" (arXiv:2003.02547)

本文件对应 RSS26 论文 Algorithm 1 line 5 "Solve u^(m+1) = S(û) with a convex solver"
的实现: 将 MATLAB 端构造好的 dense QCQP / OCP QCQP 矩阵传入 HPIPM 求解。

接口:
    solve_qcqp(...)      -> dict  求解 dense QCQP (Golden oracle)
    solve_ocp_qcqp(...)  -> dict  求解 OCP QCQP (原生逐阶段, 与 dense 数学等价)
    solve_ocp_qp(...)    -> dict  求解 OCP QP (legacy, 线性化二次约束)
"""

import sys
import os
import time

# Windows + MATLAB: Python 环境可能包含 sjtu-agent venv 的 site-packages
# venv 的 delvewheel scipy/numpy 有 DLL 加载问题 (WinError 206)
# 在 import numpy 之前清理 venv 路径, 确保从系统 Python 加载
_venv_paths_removed = []
for _p in list(sys.path):
    if 'site-packages' in _p and 'sjtu-agent' in _p.lower():
        sys.path.remove(_p)
        _venv_paths_removed.append(_p)

# 显式添加系统 Python site-packages (scipy, cvxpy 装在这里)
_sys_site = os.path.join(sys.base_prefix, 'Lib', 'site-packages')
if os.path.isdir(_sys_site) and _sys_site not in sys.path:
    sys.path.insert(0, _sys_site)

import numpy as np

# HPIPM Python 接口
_hpipm_path = None

# 定位 HPIPM Python wrapper 路径 (third_party/hpipm/interfaces/python/)
_candidate = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'third_party', 'hpipm', 'interfaces', 'python', 'hpipm_python')
if not os.path.isdir(_candidate):
    _candidate = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'third_party', 'hpipm', 'interfaces', 'python')
if os.path.isdir(_candidate):
    if _candidate not in sys.path:
        sys.path.insert(0, _candidate)
    _hpipm_path = _candidate

# Windows: 把 libhpipm.dll 所在目录加入 DLL 搜索路径
_hpipm_lib_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'third_party', 'hpipm', 'lib')
_hpipm_lib_dir = os.path.abspath(_hpipm_lib_dir)
if os.path.isdir(_hpipm_lib_dir):
    if sys.platform.startswith('win'):
        try:
            os.add_dll_directory(_hpipm_lib_dir)
        except (OSError, FileNotFoundError):
            pass
        _path_env = os.environ.get('PATH', '')
        if _hpipm_lib_dir not in _path_env.split(os.pathsep):
            os.environ['PATH'] = _hpipm_lib_dir + os.pathsep + _path_env
    else:
        _ld = os.environ.get('LD_LIBRARY_PATH', '')
        if _hpipm_lib_dir not in _ld.split(os.pathsep):
            os.environ['LD_LIBRARY_PATH'] = _hpipm_lib_dir + os.pathsep + _ld

# 导入 HPIPM dense QCQP 接口类
try:
    from hpipm_python import (
        hpipm_dense_qcqp_dim,
        hpipm_dense_qcqp,
        hpipm_dense_qcqp_sol,
        hpipm_dense_qcqp_solver_arg,
        hpipm_dense_qcqp_solver,
    )
    _HPIPM_OK = True
except Exception as _e:
    _HPIPM_OK = False
    _HPIPM_ERR = str(_e)

# 导入 HPIPM OCP QCQP 接口类
try:
    from hpipm_python import (
        hpipm_ocp_qcqp_dim,
        hpipm_ocp_qcqp,
        hpipm_ocp_qcqp_sol,
        hpipm_ocp_qcqp_solver_arg,
        hpipm_ocp_qcqp_solver,
    )
    _HPIPM_OCP_OK = True
except Exception as _e:
    _HPIPM_OCP_OK = False
    _HPIPM_OCP_ERR = str(_e)

# 导入 HPIPM OCP QP 接口类 (legacy, 线性约束)
try:
    from hpipm_python import (
        hpipm_ocp_qp_dim,
        hpipm_ocp_qp,
        hpipm_ocp_qp_sol,
        hpipm_ocp_qp_solver_arg,
        hpipm_ocp_qp_solver,
    )
    _HPIPM_OCP_QP_OK = True
except Exception as _e:
    _HPIPM_OCP_QP_OK = False
    _HPIPM_OCP_QP_ERR = str(_e)


def solve_qcqp(H, g, A, b, Hq, gq, uq, verbose=False):
    """
    求解 dense QCQP (Golden oracle).

    HPIPM dense QCQP 标准形式 (硬约束, 无 slack):
        min  0.5 x^T H x + g^T x
        s.t. A x = b                                  (等式, 论文 (20c) 动力学)
             0.5 x^T Hq_i x + gq_i^T x <= uq_i        (二次不等式, 72 条)

    返回 dict:
        x, status, status_str, obj_value, solve_time, iters
    """
    if not _HPIPM_OK:
        raise RuntimeError(f"HPIPM 不可用: {_HPIPM_ERR}")

    H = np.asarray(H, dtype=np.float64)
    g = np.asarray(g, dtype=np.float64).flatten()
    n = H.shape[0]

    has_eq = (A is not None and b is not None
              and hasattr(A, '__len__') and len(b) > 0)
    if has_eq:
        A = np.asarray(A, dtype=np.float64)
        b = np.asarray(b, dtype=np.float64).flatten()
        ne = A.shape[0]
    else:
        ne = 0

    if Hq is not None and len(Hq) > 0:
        if isinstance(Hq, np.ndarray) and Hq.ndim == 3:
            Hq_list = [Hq[:, :, i] for i in range(Hq.shape[2])]
        elif isinstance(Hq, np.ndarray) and Hq.ndim == 2:
            nq = Hq.shape[1] // n
            Hq_list = [Hq[:, i*n:(i+1)*n] for i in range(nq)]
        else:
            Hq_list = [np.asarray(h, dtype=np.float64) for h in Hq]
        nq = len(Hq_list)

        if isinstance(gq, np.ndarray) and gq.ndim == 2:
            gq_list = [gq[:, i] for i in range(gq.shape[1])]
        else:
            gq_list = [np.asarray(gq_i, dtype=np.float64).flatten() for gq_i in gq]

        uq_arr = np.asarray(uq, dtype=np.float64).flatten()
    else:
        nq = 0
        Hq_list = []
        gq_list = []
        uq_arr = np.zeros(0)

    dim = hpipm_dense_qcqp_dim()
    dim.set('nv', n)
    dim.set('ne', ne)
    dim.set('nb', 0)
    dim.set('ng', 0)
    dim.set('nq', nq)

    qcqp = hpipm_dense_qcqp(dim)
    qcqp.set('H', H)
    qcqp.set('g', g)

    if ne > 0:
        qcqp.set('A', A)
        qcqp.set('b', b)

    if nq > 0:
        Hq_stacked = np.hstack([np.asarray(hq, dtype=np.float64) for hq in Hq_list])
        gq_stacked = np.column_stack(gq_list)
        qcqp.set('Hq', Hq_stacked)
        qcqp.set('gq', gq_stacked)
        qcqp.set('uq', uq_arr)

    qcqp_sol = hpipm_dense_qcqp_sol(dim)

    # 容差 1e-8 (对齐论文 benchmark: ECOS 默认容差 abstol/reltol/feastol=1e-8)
    arg = hpipm_dense_qcqp_solver_arg(dim, 'balance')
    arg.set('iter_max', 1000)
    arg.set('tol_stat', 1e-8)
    arg.set('tol_eq', 1e-8)
    arg.set('tol_ineq', 1e-8)
    arg.set('tol_comp', 1e-8)

    solver = hpipm_dense_qcqp_solver(dim, arg)
    t0 = time.perf_counter()
    solver.solve(qcqp, qcqp_sol)
    solve_time = time.perf_counter() - t0

    x = qcqp_sol.get('v').flatten()
    status = int(solver.get('status'))
    iters = int(solver.get('iter')) if hasattr(solver, 'get') else 0
    obj_value = float(0.5 * x @ H @ x + g @ x)

    if status != 0:
        if verbose:
            print(f"[hpipm] balance 失败 (status={status}), 重试 robust...", file=sys.stderr)
        arg2 = hpipm_dense_qcqp_solver_arg(dim, 'robust')
        arg2.set('iter_max', 2000)
        arg2.set('tol_stat', 1e-8)
        arg2.set('tol_eq', 1e-8)
        arg2.set('tol_ineq', 1e-8)
        arg2.set('tol_comp', 1e-8)
        solver2 = hpipm_dense_qcqp_solver(dim, arg2)
        t1 = time.perf_counter()
        solver2.solve(qcqp, qcqp_sol)
        solve_time += time.perf_counter() - t1
        x = qcqp_sol.get('v').flatten()
        status = int(solver2.get('status'))
        obj_value = float(0.5 * x @ H @ x + g @ x)

    status_str = 'Solved' if status == 0 else f'Failed({status})'

    return {
        'x': x,
        'status': status,
        'status_str': status_str,
        'obj_value': obj_value,
        'solve_time': float(solve_time),
        'iters': iters,
    }


def solve_ocp_qcqp(A, B, b_stack, Q_stack, S_stack, R_stack, q_stack, r_stack,
                   nx_arr, nu_arr, nq_arr, nbx_arr, nq_per_stage,
                   Qq_stack, Sq_stack, Rq_stack, qq_stack, rq_stack, uq_stack,
                   x0, idxbx, const=0.0, verbose=False):
    """
    求解 OCP QCQP (HPIPM ocp_qcqp 接口, 逐阶段结构).

    与 solve_qcqp (dense QCQP) 数学等价, 但利用 OCP 块三对角结构.

    HPIPM OCP QCQP 标准形式 (每 stage n=0..N):
        动力学:  x_{n+1} = A_n x_n + B_n u_n + b_n      (n=0..N-1)
        代价:    min Σ_{n=0}^{N-1} [0.5 x'Q_n x + x'S_n' u + 0.5 u'R_n u + q_n'x + r_n'u]
                     + 0.5 x_N'Q_N x_N + q_N'x_N
        二次约束(每 stage nq 条): 0.5 x'Qq x + x'Sq' u + 0.5 u'Rq u + qq'x + rq'u <= uq
        box 约束: lbx <= x[idxbx] <= ubx

    约定 (与 construct_ocp_qcqp_from_rss.m 一致):
        - 带 HPIPM 0.5 前缀的矩阵 (Q, R, Qq, Rq): 设为 2 × (数学系数)
          例: 若数学中 v'Mv (系数 1), 则 Qq(4:6,4:6) = 2*M; 若 0.5*u'Mu (系数 0.5), 则 Rq = M
        - 不带 0.5 前缀的矩阵/向量 (S, q, r, Sq, qq, rq): 直接设为数学系数
        - 标量 (uq, const): 直接设为数学值, 不×2

    参数:
        A, B        : 动力学矩阵 (常数, 所有 stage 相同)
                      A: (nx, nx), B: (nx, nu)
        b_stack     : (nx, K) 动力学 bias, 每列一个 stage (n=0..K-1)
        Q_stack     : (nx, nx*N_stages) 逐 stage 状态代价 Hessian (已×2)
        S_stack     : (nu, nx*K) 逐 stage 交叉项 (已×2, 终端无 S)
        R_stack     : (nu, nu*K) 逐 stage 控制代价 Hessian (已×2)
        q_stack     : (nx, N_stages) 逐 stage 状态线性项 (已×2)
        r_stack     : (nu, K) 逐 stage 控制线性项 (已×2)
        nx_arr      : (N_stages,) int, 每 stage 状态维数
        nu_arr      : (N_stages,) int, 每 stage 控制维数
        nq_arr      : (N_stages,) int, 每 stage 二次约束数 (dim 用, = nq_per_stage)
        nbx_arr     : (N_stages,) int, 每 stage box 约束数
        nq_per_stage: (N_stages,) int, 每 stage 实际非空二次约束数
        Qq_stack    : (nx, nx, total_qcqp) 3D, 按 stage 顺序排列所有约束的状态二次项 (已×2)
        Sq_stack    : (nu, nx, total_qcqp) 3D, 状态-控制交叉 (HPIPM 原生方向 nu×nx, 已×2)
                      注: 终端 stage (nu=0) 的约束不应出现在 Sq_stack 中 (其 slice 为空)
        Rq_stack    : (nu, nu, total_qcqp) 3D, 控制二次项 (已×2)
        qq_stack    : (nx, total_qcqp) 2D, 状态线性项 (已×2)
        rq_stack    : (nu, total_qcqp) 2D, 控制线性项 (已×2)
        uq_stack    : (total_qcqp,) 1D, 约束上界 (不×2)
        x0          : (nx,) 初始状态 (box bounds: lbx=ubx=x0)
        idxbx       : (nbx[0],) int, box 约束索引 (0-indexed, 索引 x 内位置)
        const       : float, 常数项 (不×2, 仅用于 obj_value_full 比较)

    返回 dict:
        x                 : (36,) [u(:); nu(:)]
        status            : int      0=成功
        status_str        : str
        obj_value_full    : float    完整目标 (含 const)
        obj_value_reduced : float    去除 const 的目标 (与 Dense solve_qcqp 日志一致)
        solve_time        : float
        iters             : int
        solver_call_count : int      Golden 路径恒为 1
        max_res_stat      : float
        max_res_eq        : float
        max_res_ineq      : float
        max_res_comp      : float
    """
    if not _HPIPM_OCP_OK:
        raise RuntimeError(
            f"HPIPM OCP QCQP 不可用: {_HPIPM_OCP_ERR}\n\n"
            "需要 HPIPM 编译时启用 OCP QCQP 支持 (d_ocp_qcqp_* 符号)."
        )

    # === 统一转为 numpy float64 / int32 ===
    # 注意: MATLAB 传 (1,N) 或 (N,1) 数组到 Python 时, numpy 可能压缩为 1D (N,)
    # 对于已知维度的矩阵, 先获取维度信息再 reshape
    nx_arr = np.asarray(nx_arr, dtype=np.int32).flatten()
    nu_arr = np.asarray(nu_arr, dtype=np.int32).flatten()
    nq_arr = np.asarray(nq_arr, dtype=np.int32).flatten()
    nbx_arr = np.asarray(nbx_arr, dtype=np.int32).flatten()
    nq_per_stage = np.asarray(nq_per_stage, dtype=np.int32).flatten()

    nx0 = int(nx_arr[0])
    nu0 = int(nu_arr[0])
    total_nq = int(nq_per_stage.sum())
    N_stages = len(nx_arr)
    K = N_stages - 1

    # A: (nx, nx), B: (nx, nu) — 根据 nx0/nu0 reshape
    A = np.asarray(A, dtype=np.float64)
    if A.size == nx0 * nx0:
        A = A.reshape((nx0, nx0), order='F')
    else:
        A = np.atleast_2d(A)
    B = np.asarray(B, dtype=np.float64)
    if B.size == nx0 * nu0:
        B = B.reshape((nx0, nu0), order='F')
    else:
        B = np.atleast_2d(B)

    # b_stack: (nx, K)
    b_stack = np.asarray(b_stack, dtype=np.float64)
    if b_stack.size == nx0 * K:
        b_stack = b_stack.reshape((nx0, K), order='F')
    else:
        b_stack = np.atleast_2d(b_stack)

    # Q_stack: (nx, nx*N_stages) 或 (nx, nx)
    Q_stack = np.asarray(Q_stack, dtype=np.float64)
    if Q_stack.size == nx0 * nx0 * N_stages:
        Q_stack = Q_stack.reshape((nx0, nx0 * N_stages), order='F')
    elif Q_stack.size == nx0 * nx0:
        Q_stack = Q_stack.reshape((nx0, nx0), order='F')
    else:
        Q_stack = np.atleast_2d(Q_stack)

    # S_stack: (nu, nx*K) 或 (nu, nx)
    S_stack = np.asarray(S_stack, dtype=np.float64)
    if S_stack.size == nu0 * nx0 * K:
        S_stack = S_stack.reshape((nu0, nx0 * K), order='F')
    elif S_stack.size == nu0 * nx0:
        S_stack = S_stack.reshape((nu0, nx0), order='F')
    else:
        S_stack = np.atleast_2d(S_stack)

    # R_stack: (nu, nu*K) 或 (nu, nu)
    R_stack = np.asarray(R_stack, dtype=np.float64)
    if R_stack.size == nu0 * nu0 * K:
        R_stack = R_stack.reshape((nu0, nu0 * K), order='F')
    elif R_stack.size == nu0 * nu0:
        R_stack = R_stack.reshape((nu0, nu0), order='F')
    else:
        R_stack = np.atleast_2d(R_stack)

    # q_stack: (nx, N_stages)
    q_stack = np.asarray(q_stack, dtype=np.float64)
    if q_stack.size == nx0 * N_stages:
        q_stack = q_stack.reshape((nx0, N_stages), order='F')
    else:
        q_stack = np.atleast_2d(q_stack)

    # r_stack: (nu, K)
    r_stack = np.asarray(r_stack, dtype=np.float64)
    if r_stack.size == nu0 * K:
        r_stack = r_stack.reshape((nu0, K), order='F')
    else:
        r_stack = np.atleast_2d(r_stack)

    # 3D stacks: 根据 total_nq reshape
    Qq_stack = np.asarray(Qq_stack, dtype=np.float64)
    if total_nq > 0 and Qq_stack.size == nx0 * nx0 * total_nq:
        Qq_stack = Qq_stack.reshape((nx0, nx0, total_nq), order='F')
    Sq_stack = np.asarray(Sq_stack, dtype=np.float64)
    if total_nq > 0 and Sq_stack.size == nu0 * nx0 * total_nq:
        Sq_stack = Sq_stack.reshape((nu0, nx0, total_nq), order='F')
    Rq_stack = np.asarray(Rq_stack, dtype=np.float64)
    if total_nq > 0 and Rq_stack.size == nu0 * nu0 * total_nq:
        Rq_stack = Rq_stack.reshape((nu0, nu0, total_nq), order='F')
    qq_stack = np.asarray(qq_stack, dtype=np.float64)
    if total_nq > 0 and qq_stack.size == nx0 * total_nq:
        qq_stack = qq_stack.reshape((nx0, total_nq), order='F')
    rq_stack = np.asarray(rq_stack, dtype=np.float64)
    if total_nq > 0 and rq_stack.size == nu0 * total_nq:
        rq_stack = rq_stack.reshape((nu0, total_nq), order='F')
    uq_stack = np.asarray(uq_stack, dtype=np.float64).flatten()
    x0 = np.asarray(x0, dtype=np.float64).flatten()
    idxbx = np.asarray(idxbx, dtype=np.int32).flatten()
    const = float(const)

    # === 维度推断 ===
    K = b_stack.shape[1]
    N_stages = K + 1
    nx0 = int(nx_arr[0])
    nu0 = int(nu_arr[0])
    total_nq = int(nq_per_stage.sum())

    # === 显式 shape 检查 + isfinite 检查 ===
    assert A.shape == (nx0, nx0), f"A shape {A.shape} != ({nx0},{nx0})"
    assert B.shape == (nx0, nu0), f"B shape {B.shape} != ({nx0},{nu0})"
    assert b_stack.shape == (nx0, K), f"b_stack shape {b_stack.shape} != ({nx0},{K})"
    # Q_stack: (nx, nx*N_stages) 逐 stage 或 (nx, nx) 向后兼容
    if Q_stack.shape == (nx0, nx0):
        Q_per_stage = False
    elif Q_stack.shape[0] == nx0 and Q_stack.shape[1] == nx0 * N_stages:
        Q_per_stage = True
    elif Q_stack.shape[1] == nx0 and Q_stack.shape[0] == nx0 * N_stages:
        Q_per_stage = True
        Q_stack = Q_stack.T.copy()
    else:
        raise ValueError(f"Q_stack shape {Q_stack.shape} 不合法")
    # S_stack: (nu, nx*K) 逐 stage 或 (nu, nx) 向后兼容
    if S_stack.shape == (nu0, nx0):
        S_per_stage = False
    elif S_stack.shape[0] == nu0 and S_stack.shape[1] == nx0 * K:
        S_per_stage = True
    elif S_stack.shape[1] == nu0 and S_stack.shape[0] == nx0 * K:
        S_per_stage = True
        S_stack = S_stack.T.copy()
    else:
        raise ValueError(f"S_stack shape {S_stack.shape} 不合法")
    # R_stack: (nu, nu*K) 逐 stage 或 (nu, nu) 向后兼容
    if R_stack.shape == (nu0, nu0):
        R_per_stage = False
    elif R_stack.shape[0] == nu0 and R_stack.shape[1] == nu0 * K:
        R_per_stage = True
    elif R_stack.shape[1] == nu0 and R_stack.shape[0] == nu0 * K:
        R_per_stage = True
        R_stack = R_stack.T.copy()
    else:
        raise ValueError(f"R_stack shape {R_stack.shape} 不合法")
    # q_stack: (nx, N_stages)
    if q_stack.ndim == 1 and q_stack.shape[0] == nx0:
        q_stack = np.zeros((nx0, N_stages))
    assert q_stack.shape == (nx0, N_stages), f"q_stack shape {q_stack.shape} != ({nx0},{N_stages})"
    # r_stack: (nu, K)
    assert r_stack.shape == (nu0, K), f"r_stack shape {r_stack.shape} != ({nu0},{K})"
    # 二次约束 stacks
    if total_nq > 0:
        assert Qq_stack.shape == (nx0, nx0, total_nq), \
            f"Qq_stack shape {Qq_stack.shape} != ({nx0},{nx0},{total_nq})"
        assert Sq_stack.shape == (nu0, nx0, total_nq), \
            f"Sq_stack shape {Sq_stack.shape} != ({nu0},{nx0},{total_nq})"
        assert Rq_stack.shape == (nu0, nu0, total_nq), \
            f"Rq_stack shape {Rq_stack.shape} != ({nu0},{nu0},{total_nq})"
        assert qq_stack.shape == (nx0, total_nq), \
            f"qq_stack shape {qq_stack.shape} != ({nx0},{total_nq})"
        assert rq_stack.shape == (nu0, total_nq), \
            f"rq_stack shape {rq_stack.shape} != ({nu0},{total_nq})"
        assert uq_stack.shape == (total_nq,), \
            f"uq_stack shape {uq_stack.shape} != ({total_nq},)"

    # isfinite 检查
    for name, arr in [('A', A), ('B', B), ('b_stack', b_stack),
                      ('Q_stack', Q_stack), ('S_stack', S_stack),
                      ('R_stack', R_stack), ('q_stack', q_stack),
                      ('r_stack', r_stack), ('x0', x0)]:
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"solve_ocp_qcqp: {name} contains NaN/Inf")
    if total_nq > 0:
        for name, arr in [('Qq_stack', Qq_stack), ('Sq_stack', Sq_stack),
                          ('Rq_stack', Rq_stack), ('qq_stack', qq_stack),
                          ('rq_stack', rq_stack), ('uq_stack', uq_stack)]:
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"solve_ocp_qcqp: {name} contains NaN/Inf")

    # === HPIPM 维度设置 ===
    # 必须显式设置所有 dim 字段 (nx, nu, nbx, nbu, ng, nq, ns, nbxe, nbue, nge, nqe)
    # 官方示例 ocp_qcqp_data.c 明确将 nbu/ng/ns/nbxe/nbue/nge/nqe 全部设为 0
    # 若不显式设置, d_ocp_qcqp_dim_create 可能不初始化这些字段, 导致 IPM 访问垃圾内存 → status=3 (NAN_SOL)
    dim = hpipm_ocp_qcqp_dim(K)
    for s in range(N_stages):
        dim.set('nx', int(nx_arr[s]), s)
        dim.set('nu', int(nu_arr[s]), s)
        dim.set('nbx', int(nbx_arr[s]), s)
        dim.set('nbu', 0, s)           # 无控制 box 约束
        dim.set('ng', 0, s)            # 无一般约束
        dim.set('nq', int(nq_arr[s]), s)
        dim.set('ns', 0, s)            # 无 soft 约束
        dim.set('nbxe', 0, s)          # 无 box equality (用 lbx=ubx=x0 代替)
        dim.set('nbue', 0, s)
        dim.set('nge', 0, s)
        dim.set('nqe', 0, s)

    # === HPIPM OCP QCQP 数据设置 (逐 stage) ===
    qp = hpipm_ocp_qcqp(dim)

    # 动力学 (n=0..K-1): A, B, b
    for n in range(K):
        qp.set('A', A, n)
        qp.set('B', B, n)
        qp.set('b', b_stack[:, n], n)

    # 代价 (逐 stage Q, q; n<K: R, r, S)
    for n in range(N_stages):
        Q_n = Q_stack[:, n*nx0:(n+1)*nx0] if Q_per_stage else Q_stack
        qp.set('Q', Q_n, n)
        qp.set('q', q_stack[:, n], n)
        if n < K:
            S_n = S_stack[:, n*nx0:(n+1)*nx0] if S_per_stage else S_stack
            R_n = R_stack[:, n*nu0:(n+1)*nu0] if R_per_stage else R_stack
            qp.set('S', S_n, n)
            qp.set('R', R_n, n)
            qp.set('r', r_stack[:, n], n)

    # 二次约束 (逐 stage)
    offset = 0
    for s in range(N_stages):
        nq_s = int(nq_per_stage[s])
        if nq_s > 0:
            nx_s = int(nx_arr[s])
            nu_s = int(nu_arr[s])
            # Qq: 水平堆叠 (nx, nx*nq_s)
            Qq_s = np.hstack([Qq_stack[:, :, offset + j] for j in range(nq_s)])
            # qq: 列堆叠 (nx, nq_s)
            qq_s = np.column_stack([qq_stack[:, offset + j] for j in range(nq_s)])
            # uq: (nq_s,)
            uq_s = uq_stack[offset:offset + nq_s]
            qp.set('Qq', Qq_s, s)
            qp.set('qq', qq_s, s)
            qp.set('uq', uq_s, s)
            # Sq/Rq/rq: 仅当 nu_s > 0 时设置 (终端 stage nu=0 时跳过这些 setter)
            if nu_s > 0:
                # Sq: HPIPM C 层 OCP_QCQP_SET_SQ 用 CVT_TRAN_MAT2STRMAT
                # 读取列优先 (nu, nx) 矩阵, 转置后存入 Hq 下左块 (nx, nu)
                # Python wrapper 用 np.ravel('F') 列优先展开
                # Sq_stack 已经是 (nu, nx, nq) 原生方向, 直接 reshape 为 (nu, nx*nq_s)
                Sq_s = Sq_stack[:, :, offset:offset + nq_s].reshape((nu_s, nx_s * nq_s), order='F')
                # Rq: (nu, nu*nq_s) 水平堆叠
                Rq_s = np.hstack([Rq_stack[:, :, offset + j] for j in range(nq_s)])
                # rq: (nu, nq_s) 列堆叠
                rq_s = np.column_stack([rq_stack[:, offset + j] for j in range(nq_s)])
                qp.set('Sq', Sq_s, s)
                qp.set('Rq', Rq_s, s)
                qp.set('rq', rq_s, s)
            # 终端 stage (nu=0): 不调用 Sq/Rq/rq setter (按 HPIPM 官方示例)
        offset += nq_s

    # 初始状态 box 约束 (stage 0: lbx=ubx=x0, idxbx=0..5)
    if int(nbx_arr[0]) > 0:
        qp.set('idxbx', idxbx, 0)
        qp.set('lbx', x0, 0)
        qp.set('ubx', x0, 0)

    # === 诊断 (verbose) ===
    if verbose:
        print(f"\n[hpipm_ocp_qcqp_diag] K={K}, N_stages={N_stages}", file=sys.stderr)
        print(f"[hpipm_ocp_qcqp_diag] nx_arr={nx_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qcqp_diag] nu_arr={nu_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qcqp_diag] nq_arr={nq_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qcqp_diag] nbx_arr={nbx_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qcqp_diag] nq_per_stage={nq_per_stage.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qcqp_diag] total_nq={total_nq}", file=sys.stderr)
        if total_nq > 0:
            # 检查每条约束 Hessian 的对称性和 PSD 性
            for j in range(min(total_nq, 8)):
                Qq_j = Qq_stack[:, :, j]
                sym_err = np.max(np.abs(Qq_j - Qq_j.T))
                eig_min = np.linalg.eigvalsh(0.5 * (Qq_j + Qq_j.T)).min()
                print(f"[hpipm_ocp_qcqp_diag] 约束 {j}: uq={uq_stack[j]:.6e}, sym_err={sym_err:.3e}, "
                      f"Qq eig_min={eig_min:.3e}", file=sys.stderr)
            # Hessian PSD 检查 (允许 -1e-12 浮点误差)
            for j in range(total_nq):
                Qq_j = Qq_stack[:, :, j]
                eig_min = np.linalg.eigvalsh(0.5 * (Qq_j + Qq_j.T)).min()
                if eig_min < -1e-10:
                    print(f"[hpipm_ocp_qcqp_diag] !!! 约束 {j} Qq 非 PSD: eig_min={eig_min:.3e} !!!",
                          file=sys.stderr)

    # === 求解 (Golden 路径: 1 次 solve, 不重试) ===
    qp_sol = hpipm_ocp_qcqp_sol(dim)

    # 求解器模式从环境变量读取 (默认 speed)
    # speed=1, balance=2, robust=3
    _ocp_qcqp_mode = os.environ.get('HPIPM_OCP_QCQP_MODE', 'speed')
    arg = hpipm_ocp_qcqp_solver_arg(dim, _ocp_qcqp_mode)
    arg.set('iter_max', 1000)
    arg.set('tol_stat', 1e-8)
    arg.set('tol_eq', 1e-8)
    arg.set('tol_ineq', 1e-8)
    arg.set('tol_comp', 1e-8)
    arg.set('mu0', 10.0)

    # warm_start=2: 选中"全量裁剪初始化"分支, 绕过 DLL C 层 bug
    # (d_ocp_qcqp_ipm_arg_set_t0_init 误写 t_lam_min 而非 t0_init,
    #  导致 nq>0 时零初始化除零 → status=3; warm_start=2 是实测唯一可行组合)
    arg.set('warm_start', 2)

    # verbose: 打印完整 C 结构体 (诊断 status=3)
    if verbose:
        print("\n[hpipm_ocp_qcqp_diag] === qp.print_C_struct() ===", file=sys.stderr)
        try:
            qp.print_C_struct()
        except Exception as e:
            print(f"[hpipm_ocp_qcqp_diag] print_C_struct failed: {e}", file=sys.stderr)

    solver = hpipm_ocp_qcqp_solver(dim, arg)
    t0 = time.perf_counter()
    solver.solve(qp, qp_sol)
    solve_time = time.perf_counter() - t0

    status = int(solver.get('status'))
    iters = int(solver.get('iter'))
    solver_call_count = 1  # Golden 路径: 恒为 1

    # 残差 (尽力提取, 不所有版本都支持)
    max_res_stat = float('nan')
    max_res_eq = float('nan')
    max_res_ineq = float('nan')
    max_res_comp = float('nan')
    try:
        max_res_stat = float(solver.get('max_res_stat'))
        max_res_eq = float(solver.get('max_res_eq'))
        max_res_ineq = float(solver.get('max_res_ineq'))
        max_res_comp = float(solver.get('max_res_comp'))
    except Exception:
        pass

    if verbose:
        print(f"[hpipm_ocp_qcqp_diag] status={status}, iters={iters}, "
              f"res_stat={max_res_stat:.3e}, res_eq={max_res_eq:.3e}, "
              f"res_ineq={max_res_ineq:.3e}, res_comp={max_res_comp:.3e}", file=sys.stderr)

    # 提取解: u_n (n=0..K-1), x_n (n=0..K)
    u_list = []
    x_list = []
    for n in range(N_stages):
        x_n = np.asarray(qp_sol.get('x', n)).flatten()
        x_list.append(x_n)
        if n < K:
            u_n = np.asarray(qp_sol.get('u', n)).flatten()
            u_list.append(u_n)

    # 拼成 [u(:); nu(:)] 兼容 control_RSS 的提取逻辑
    # nu(:) = 速度部分, 对应 x_n 的后 3 维 (RSS: x_n = [e_n; v_n], v_n = x_n[3:6])
    # 对 nx != 6 的情况, 取后 min(3, nx) 维; 若 nx <= 3, 取全部
    u_flat = np.concatenate(u_list) if u_list else np.zeros(0)
    v_dim = min(3, nx0) if nx0 > 3 else nx0
    v_start = nx0 - v_dim  # RSS: 6-3=3, 对应 0-indexed [3:6]
    nu_flat = np.concatenate([x_list[k][v_start:v_start+v_dim] for k in range(1, N_stages)])
    x_out = np.concatenate([u_flat, nu_flat])

    # 手动计算 obj_value_reduced (Q, S, R, q, r 已×2; 0.5 前缀还原; 不含 const)
    obj_value_reduced = 0.0
    for n in range(N_stages):
        xn = x_list[n]
        Q_n = Q_stack[:, n*nx0:(n+1)*nx0] if Q_per_stage else Q_stack
        obj_value_reduced += 0.5 * float(xn @ Q_n @ xn)
        obj_value_reduced += float(q_stack[:, n] @ xn)
        if n < K:
            un = u_list[n]
            S_n = S_stack[:, n*nx0:(n+1)*nx0] if S_per_stage else S_stack
            R_n = R_stack[:, n*nu0:(n+1)*nu0] if R_per_stage else R_stack
            obj_value_reduced += float(xn @ S_n.T @ un)  # x'S'u
            obj_value_reduced += 0.5 * float(un @ R_n @ un)
            obj_value_reduced += float(r_stack[:, n] @ un)

    obj_value_full = obj_value_reduced + const

    status_str = 'Solved' if status == 0 else f'Failed({status})'

    return {
        'x': x_out,
        'x_stages': x_list,  # list of x_n (n=0..K), 每个 (nx,)
        'u_stages': u_list,  # list of u_n (n=0..K-1), 每个 (nu,)
        'status': status,
        'status_str': status_str,
        'obj_value_full': float(obj_value_full),
        'obj_value_reduced': float(obj_value_reduced),
        'solve_time': float(solve_time),
        'iters': iters,
        'solver_call_count': solver_call_count,
        'max_res_stat': float(max_res_stat),
        'max_res_eq': float(max_res_eq),
        'max_res_ineq': float(max_res_ineq),
        'max_res_comp': float(max_res_comp),
    }


def solve_ocp_qp(A, B, Q_eff, R_eff, S_eff, b_stack, r_stack, q_stack,
                 nx_arr, nu_arr, ng_arr, nbx_arr, ng_per_stage,
                 Cmat_stack, Dmat_stack, lg_stack, ug_stack,
                 x0, idxbx, const=0.0, verbose=False,
                 x_init=None, u_init=None, warm_start=False,
                 allow_retry=True):
    """
    求解 OCP QP (HPIPM ocp_qp 接口, 一般线性约束).

    Legacy 求解路径: 二次约束在 MATLAB 端做一阶泰勒线性化, 转为一般线性约束.
    三次 OCP QP 是对原凸二次子问题的有限切平面近似;
    可在特定问题上接近 Dense QCQP, 但不保证一般性严格等价.

    HPIPM OCP QP 标准形式 (每 stage n=0..N):
        动力学:  x_{n+1} = A_n x_n + B_n u_n + b_n      (n=0..N-1)
        代价:    min Σ [0.5 x'Qx + 0.5 u'Ru + r'u] + 0.5 x_N'Q_N x_N
        一般线性约束 (每 stage ng 条): lg <= C*x + D*u <= ug
        box 约束: lbx <= x[idxbx] <= ubx

    ×2 约定 (与 construct_ocp_qp_from_rss.m 一致):
        Q_eff, R_eff, r 构造时已×2; HPIPM 的 0.5 前缀使其还原为原始系数.
        C, D, lg, ug, const 不×2 (线性约束和常数项无 1/2 前缀).

    返回 dict:
        x, status, status_str, obj_value, solve_time, iters, solver_call_count
    """
    if not _HPIPM_OCP_QP_OK:
        raise RuntimeError(
            f"HPIPM OCP QP 不可用: {_HPIPM_OCP_QP_ERR}\n\n"
            "需要 HPIPM 编译时启用 OCP QP 支持 (d_ocp_qp_* 符号)."
        )

    # === 统一转为 numpy float64 / int32 ===
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    Q_eff = np.asarray(Q_eff, dtype=np.float64)
    R_eff = np.asarray(R_eff, dtype=np.float64)
    S_eff = np.asarray(S_eff, dtype=np.float64)
    q_stack = np.asarray(q_stack, dtype=np.float64)
    b_stack = np.asarray(b_stack, dtype=np.float64)
    r_stack = np.asarray(r_stack, dtype=np.float64)
    nx_arr = np.asarray(nx_arr, dtype=np.int32).flatten()
    nu_arr = np.asarray(nu_arr, dtype=np.int32).flatten()
    ng_arr = np.asarray(ng_arr, dtype=np.int32).flatten()
    nbx_arr = np.asarray(nbx_arr, dtype=np.int32).flatten()
    ng_per_stage = np.asarray(ng_per_stage, dtype=np.int32).flatten()
    Cmat_stack = np.asarray(Cmat_stack, dtype=np.float64)
    Dmat_stack = np.asarray(Dmat_stack, dtype=np.float64)
    lg_stack = np.asarray(lg_stack, dtype=np.float64).flatten()
    ug_stack = np.asarray(ug_stack, dtype=np.float64).flatten()
    x0 = np.asarray(x0, dtype=np.float64).flatten()
    idxbx = np.asarray(idxbx, dtype=np.int32).flatten()

    K = b_stack.shape[1]                  # horizon (=6)
    N_stages = K + 1                      # n=0..K
    nx0 = int(nx_arr[0])                  # 6
    nu0 = int(nu_arr[0])                  # 3

    # === HPIPM 维度设置 ===
    dim = hpipm_ocp_qp_dim(K)             # N=K (horizon)
    for s in range(N_stages):
        dim.set('nx', int(nx_arr[s]), s)
        dim.set('nu', int(nu_arr[s]), s)
        dim.set('nbx', int(nbx_arr[s]), s)
        dim.set('ng', int(ng_arr[s]), s)

    # === HPIPM OCP QP 数据设置 (逐 stage) ===
    qp = hpipm_ocp_qp(dim)

    # 动力学 (n=0..K-1): A, B, b
    for n in range(K):
        qp.set('A', A, n)
        qp.set('B', B, n)
        qp.set('b', b_stack[:, n], n)

    # 代价 (n=0..K): Q (逐stage可不同), q=0; (n=0..K-1): R, r, S=0; 终端 n=K: 仅 Q, q
    if Q_eff.shape == (nx0, nx0):
        Q_eff_per_stage = False
    elif Q_eff.shape[0] == nx0 and Q_eff.shape[1] == 6 * N_stages:
        Q_eff_per_stage = True   # (6, 42) 行堆叠
    elif Q_eff.shape[1] == nx0 and Q_eff.shape[0] == 6 * N_stages:
        Q_eff_per_stage = True   # (42, 6) 列堆叠 (MATLAB 转置)
        Q_eff = Q_eff.T          # 转回 (6, 42)
    else:
        Q_eff_per_stage = False
    for n in range(N_stages):
        Q_n = Q_eff[:, n*nx0:(n+1)*nx0] if Q_eff_per_stage else Q_eff
        qp.set('Q', Q_n, n)
        qp.set('q', q_stack[:, n], n)
        if n < K:
            qp.set('R', R_eff, n)
            qp.set('r', r_stack[:, n], n)
            qp.set('S', S_eff, n)

    # 一般线性约束 (逐 stage, 按 ng_per_stage 切片 Cmat_stack/Dmat_stack)
    offset = 0
    for s in range(N_stages):
        ng_s = int(ng_per_stage[s])
        if ng_s > 0:
            nx_s = int(nx_arr[s])
            nu_s = int(nu_arr[s])
            C_s = Cmat_stack[offset:offset + ng_s, :nx_s]
            qp.set('C', C_s, s)
            if nu_s > 0:
                D_s = Dmat_stack[offset:offset + ng_s, :nu_s]
                qp.set('D', D_s, s)
            lg_s = lg_stack[offset:offset + ng_s]
            ug_s = ug_stack[offset:offset + ng_s]
            qp.set('lg', lg_s, s)
            qp.set('ug', ug_s, s)
        offset += ng_s

    # 初始状态 box 约束 (stage 0: lbx=ubx=x0, idxbx=0..5)
    if int(nbx_arr[0]) > 0:
        qp.set('idxbx', idxbx, 0)
        qp.set('lbx', x0, 0)
        qp.set('ubx', x0, 0)

    # === 诊断 ===
    if verbose:
        print(f"\n[hpipm_ocp_qp_diag] K={K}, N_stages={N_stages}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] nx_arr={nx_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] nu_arr={nu_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] ng_arr={ng_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] nbx_arr={nbx_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] ng_per_stage={ng_per_stage.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] total_ng={int(ng_per_stage.sum())}", file=sys.stderr)
        for name, arr in [('A', A), ('B', B), ('Q_eff', Q_eff), ('R_eff', R_eff),
                          ('b_stack', b_stack), ('r_stack', r_stack),
                          ('x0', x0), ('Cmat_stack', Cmat_stack),
                          ('Dmat_stack', Dmat_stack), ('lg_stack', lg_stack),
                          ('ug_stack', ug_stack)]:
            if not np.all(np.isfinite(arr)):
                print(f"[hpipm_ocp_qp_diag] !!! {name} contains NaN/Inf !!!", file=sys.stderr)

    # === 求解 ===
    qp_sol = hpipm_ocp_qp_sol(dim)

    arg = hpipm_ocp_qp_solver_arg(dim, 'balance')
    arg.set('iter_max', 2000)
    arg.set('tol_stat', 1e-6)
    arg.set('tol_eq', 1e-6)
    arg.set('tol_ineq', 1e-6)
    arg.set('tol_comp', 1e-6)

    # Warm start
    if warm_start and x_init is not None and u_init is not None:
        try:
            x_init_arr = np.asarray(x_init, dtype=np.float64)
            u_init_arr = np.asarray(u_init, dtype=np.float64)
            for n in range(N_stages):
                qp_sol.set('x', x_init_arr[:, n], n)
                if n < K:
                    qp_sol.set('u', u_init_arr[:, n], n)
            arg.set('warm_start', True)
            if verbose:
                print(f"[hpipm_ocp_qp] warm start enabled", file=sys.stderr)
        except Exception as e:
            if verbose:
                print(f"[hpipm_ocp_qp] warm start failed: {e}", file=sys.stderr)

    solver = hpipm_ocp_qp_solver(dim, arg)
    t0 = time.perf_counter()
    solver.solve(qp, qp_sol)
    solve_time = time.perf_counter() - t0
    try:
        st = float(arg.get('solve_time'))
        if st > 0:
            solve_time = st
    except Exception:
        pass

    status = int(solver.get('status'))
    iters = int(solver.get('iter'))
    if verbose:
        print(f"[hpipm_ocp_qp_diag] status={status}, iters={iters}", file=sys.stderr)

    # 提取解
    u_list = []
    x_list = []
    for n in range(N_stages):
        x_n = np.asarray(qp_sol.get('x', n)).flatten()
        x_list.append(x_n)
        if n < K:
            u_n = np.asarray(qp_sol.get('u', n)).flatten()
            u_list.append(u_n)

    u_flat = np.concatenate(u_list) if u_list else np.zeros(0)
    nu_flat = np.concatenate([x_list[k][3:6] for k in range(1, N_stages)])
    x_out = np.concatenate([u_flat, nu_flat])

    # 手动计算 obj_value
    obj_value = const
    for n in range(N_stages):
        xn = x_list[n]
        Q_n = Q_eff[:, n*nx0:(n+1)*nx0] if Q_eff_per_stage else Q_eff
        obj_value += 0.5 * float(xn @ Q_n @ xn)
        if n < K:
            un = u_list[n]
            obj_value += 0.5 * float(un @ R_eff @ un)
            obj_value += float(r_stack[:, n] @ un)

    # balance 失败回退 robust (仅对严重错误 status>=2, 不对 MAX_ITER status=1 重试)
    solver_call_count = 1
    if status >= 2 and allow_retry:
        if verbose:
            print(f"[hpipm_ocp_qp] balance 严重失败 (status={status}), 重试 robust...", file=sys.stderr)
        arg2 = hpipm_ocp_qp_solver_arg(dim, 'robust')
        arg2.set('iter_max', 5000)
        arg2.set('tol_stat', 1e-8)
        arg2.set('tol_eq', 1e-8)
        arg2.set('tol_ineq', 1e-8)
        arg2.set('tol_comp', 1e-8)
        solver2 = hpipm_ocp_qp_solver(dim, arg2)
        t1 = time.perf_counter()
        solver2.solve(qp, qp_sol)
        solve_time += time.perf_counter() - t1
        try:
            st2 = float(arg2.get('solve_time'))
            if st2 > 0:
                solve_time = solve_time - (time.perf_counter() - t1) + st2
        except Exception:
            pass
        status = int(solver2.get('status'))
        iters = int(solver2.get('iter'))
        solver_call_count = 2
        u_list = []
        x_list = []
        for n in range(N_stages):
            x_n = np.asarray(qp_sol.get('x', n)).flatten()
            x_list.append(x_n)
            if n < K:
                u_n = np.asarray(qp_sol.get('u', n)).flatten()
                u_list.append(u_n)
        u_flat = np.concatenate(u_list) if u_list else np.zeros(0)
        nu_flat = np.concatenate([x_list[k][3:6] for k in range(1, N_stages)])
        x_out = np.concatenate([u_flat, nu_flat])
        obj_value = const
        for n in range(N_stages):
            xn = x_list[n]
            Q_n = Q_eff[:, n*nx0:(n+1)*nx0] if Q_eff_per_stage else Q_eff
            obj_value += 0.5 * float(xn @ Q_n @ xn)
            if n < K:
                un = u_list[n]
                obj_value += 0.5 * float(un @ R_eff @ un)
                obj_value += float(r_stack[:, n] @ un)

    status_str = 'Solved' if status == 0 else f'Failed({status})'

    return {
        'x': x_out,
        'status': status,
        'status_str': status_str,
        'obj_value': float(obj_value),
        'solve_time': float(solve_time),
        'iters': iters,
        'solver_call_count': solver_call_count,
    }
