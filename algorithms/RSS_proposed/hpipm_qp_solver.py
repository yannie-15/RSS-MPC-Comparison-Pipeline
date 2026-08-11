"""
hpipm_qp_solver.py
HPIPM dense QCQP 求解器封装 (Python 接口)

论文: "hpipm: a high-performance quadratic programming framework for
       model predictive control" (arXiv:2003.02547)

本文件对应 RSS26 论文 Algorithm 1 line 5 "Solve u^(m+1) = S(û) with a convex solver"
的实现: 将 MATLAB 端构造好的 dense QCQP 矩阵传入 HPIPM 求解。

接口:
    solve_qcqp(H, g, A, b, Hq, gq, uq) -> dict(x, status, obj_value, solve_time, iters)

HPIPM dense QCQP 标准形式:

[HPIPM 论文 Section 2.1 公式 (1) — 完整 dense QP (线性约束, 含 slack)]
    min_{v,s}  1/2 [v;1]^T [H g; g^T 0] [v;1]
              + 1/2 [s^l;s^u;1]^T [Z^l 0 z^l; 0 Z^u z^u; (z^l)^T (z^u)^T 0] [s^l;s^u;1]
    s.t. A v = b                                                              (等式)
         [v_; d_] <= [J^{b,v}; C] v + [J^{s,v}; J^{s,g}] s^l                  (下界+slack)
         [J^{b,v}; C] v - [J^{s,v}; J^{s,g}] s^u <= [v^; d^]                 (上界+slack)
         s^l >= s^l_lb,  s^u >= s^u_lb                                         (slack 非负)

[dense QCQP 扩展 — 在 dense QP 基础上增加二次约束]
    0.5 v^T Hq_i v + gq_i^T v <= uq_i    (二次不等式, 亦可带 slack)

[本代码使用硬约束子集 (nb=0, ng=0, ns=0, 无 slack)]
    min  0.5 x^T H x + g^T x
    s.t. A x = b                                  (等式, 论文 (20c) 动力学)
         0.5 x^T Hq_i x + gq_i^T x <= uq_i        (二次不等式, 论文 (20a)+(20b))
即 HPIPM 的 slack/box/一般线性约束均未启用, 退化为纯 QCQP (硬约束)

HPIPM 求解器模式 (对应 HPIPM 论文 Section III):
    - 'balance': 平衡模式 (默认, 论文 IV-B 中 ECOS 默认参数的等价)
    - 'robust':  鲁棒模式 (balance 失败时回退)
"""

import sys
import time
import numpy as np

# HPIPM Python 接口
_hpipm_path = None
import os
import sys

# 定位 HPIPM Python wrapper 路径 (third_party/hpipm/interfaces/python/)
_candidate = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'third_party', 'hpipm', 'interfaces', 'python', 'hpipm_python')
if not os.path.isdir(_candidate):
    _candidate = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'third_party', 'hpipm', 'interfaces', 'python')
if os.path.isdir(_candidate):
    if _candidate not in sys.path:
        sys.path.insert(0, _candidate)
    _hpipm_path = _candidate

# Windows: 把 libhpipm.dll 所在目录加入 DLL 搜索路径
# HPIPM Python wrapper 内部用 ctypes.CDLL('libhpipm.dll') 加载共享库
_hpipm_lib_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'third_party', 'hpipm', 'lib')
_hpipm_lib_dir = os.path.abspath(_hpipm_lib_dir)
if os.path.isdir(_hpipm_lib_dir):
    if sys.platform.startswith('win'):
        # Python 3.8+: 使用 add_dll_directory 添加 DLL 搜索路径
        try:
            os.add_dll_directory(_hpipm_lib_dir)
        except (OSError, FileNotFoundError):
            pass
        # 同时加入 PATH, 作为部分 ctypes 版本的 fallback
        _path_env = os.environ.get('PATH', '')
        if _hpipm_lib_dir not in _path_env.split(os.pathsep):
            os.environ['PATH'] = _hpipm_lib_dir + os.pathsep + _path_env
    else:
        _ld = os.environ.get('LD_LIBRARY_PATH', '')
        if _hpipm_lib_dir not in _ld.split(os.pathsep):
            os.environ['LD_LIBRARY_PATH'] = _hpipm_lib_dir + os.pathsep + _ld

# 导入 HPIPM dense QCQP 接口类 (对应 HPIPM 论文中的 dense QP 数据结构)
try:
    from hpipm_python import (
        hpipm_dense_qcqp_dim,       # 维度对象 (nv, ne, nb, ng, nq)
        hpipm_dense_qcqp,           # QCQP 问题数据 (H, g, A, b, Hq, gq, uq)
        hpipm_dense_qcqp_sol,       # 解对象 (存储 v, lam 等对偶变量)
        hpipm_dense_qcqp_solver_arg,  # 求解器参数 (mode, tol, iter_max)
        hpipm_dense_qcqp_solver,    # 求解器对象 (solve 方法)
    )
    _HPIPM_OK = True
except Exception as _e:
    _HPIPM_OK = False
    _HPIPM_ERR = str(_e)

# 导入 HPIPM OCP QCQP 接口类 (对应 HPIPM ocp_qcqp 逐阶段数据结构)
try:
    from hpipm_python import (
        hpipm_ocp_qcqp_dim,         # OCP 维度对象 (N, nx, nu, nbx, ng, nq)
        hpipm_ocp_qcqp,             # OCP QCQP 问题数据 (A, B, b, Q, S, R, q, r, Qq, Sq, Rq, qq, rq, uq)
        hpipm_ocp_qcqp_sol,         # OCP 解对象 (存储 x, u 逐阶段)
        hpipm_ocp_qcqp_solver_arg,  # OCP 求解器参数 (mode, tol, iter_max)
        hpipm_ocp_qcqp_solver,      # OCP 求解器对象 (solve 方法)
    )
    _HPIPM_OCP_OK = True
except Exception as _e:
    _HPIPM_OCP_OK = False
    _HPIPM_OCP_ERR = str(_e)

# 导入 HPIPM OCP QP 接口类 (线性约束, 对应 ocp_qp 逐阶段数据结构)
# 用于替代 ocp_qcqp (OCP QCQP solver 存在 bug 导致 status=3 NAN_SOL)
try:
    from hpipm_python import (
        hpipm_ocp_qp_dim,           # OCP QP 维度对象 (N, nx, nu, nbx, ng, nbu, ns)
        hpipm_ocp_qp,               # OCP QP 问题数据 (A, B, b, Q, S, R, q, r, C, D, lg, ug)
        hpipm_ocp_qp_sol,           # OCP QP 解对象 (存储 x, u 逐阶段)
        hpipm_ocp_qp_solver_arg,    # OCP QP 求解器参数 (mode, tol, iter_max)
        hpipm_ocp_qp_solver,        # OCP QP 求解器对象 (solve 方法)
    )
    _HPIPM_OCP_QP_OK = True
except Exception as _e:
    _HPIPM_OCP_QP_OK = False
    _HPIPM_OCP_QP_ERR = str(_e)


def solve_qcqp(H, g, A, b, Hq, gq, uq, verbose=False):
    """
    求解 dense QCQP (对应 RSS26 论文 Algorithm 1 line 5 的求解步骤).

    参数 (均可为 numpy 数组或 MATLAB py.numpy.array 传入):
        H  : (n, n)   目标函数 Hessian (PSD), 论文 (18)+(19) 展开后的二次型
        g  : (n,)     目标函数线性项
        A  : (ne, n)  等式约束矩阵 (论文 (20c) 动力学递推)
        b  : (ne,)    等式约束右端
        Hq : list of (n, n)  二次约束 Hessian 列表 (论文 (20a) 转向锥 + (20b) 轮速)
        gq : list of (n,)    二次约束线性项列表
        uq : (nq,)           二次约束右端

    返回 dict:
        x          : (n,)     最优解 (x = [u(:); nu(:)], 36维)
        status     : int      0=成功 (HPIPM 论文: status=0 表示最优解找到)
        status_str : str      状态字符串
        obj_value  : float    目标函数值 (0.5*x'Hx + g'x)
        solve_time : float    求解耗时 (秒)
        iters      : int      迭代次数
    """
    if not _HPIPM_OK:
        raise RuntimeError(
            f"HPIPM 不可用: {_HPIPM_ERR}\n\n"
            "可能原因: libhpipm.dll (Windows) / libhpipm.so (Linux/Mac) 未编译.\n"
            "解决方法 (Windows):\n"
            "  1. 安装 MSYS2: https://www.msys2.org/\n"
            "  2. 在 MSYS2 UCRT64 终端安装工具链:\n"
            "       pacboy sync:mman-git ucrt64/toolchain msys/make msys/bc\n"
            "  3. 运行编译脚本 (PowerShell 中):\n"
            "       C:\\msys64\\usr\\bin\\env.exe MSYSTEM=UCRT64 /usr/bin/bash -lc \"/d/PROJECT/RSS_V2/scripts/build_hpipm_windows.sh\"\n"
            "  4. 验证: python -c \"import sys; sys.path.insert(0,'python'); import hpipm_qp_solver; print(hpipm_qp_solver._HPIPM_OK)\""
        )

    # === 统一转为 numpy float64 数组 ===
    H = np.asarray(H, dtype=np.float64)
    g = np.asarray(g, dtype=np.float64).flatten()
    n = H.shape[0]  # 决策变量维度 (36: u(18)+nu(18))

    # 等式约束 (论文 (20c): ν_{k+1}=ν_k+u_{k+1})
    has_eq = (A is not None and b is not None
              and hasattr(A, '__len__') and len(b) > 0)
    if has_eq:
        A = np.asarray(A, dtype=np.float64)
        b = np.asarray(b, dtype=np.float64).flatten()
        ne = A.shape[0]  # 等式约束数 (18)
    else:
        ne = 0

    # 二次约束 (论文 (20a) 转向锥 48 + (20b) 轮速 24 = 72)
    if Hq is not None and len(Hq) > 0:
        # 支持 list of 2D 或 3D 数组 (n, n, nq)
        if isinstance(Hq, np.ndarray) and Hq.ndim == 3:
            Hq_list = [Hq[:, :, i] for i in range(Hq.shape[2])]
        elif isinstance(Hq, np.ndarray) and Hq.ndim == 2:
            # 水平堆叠 (n, n*nq) -> 拆分
            nq = Hq.shape[1] // n
            Hq_list = [Hq[:, i*n:(i+1)*n] for i in range(nq)]
        else:
            Hq_list = [np.asarray(h, dtype=np.float64) for h in Hq]
        nq = len(Hq_list)  # 二次约束数 (72)

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

    # === HPIPM 维度设置 (对应 HPIPM 论文 Section II-B) ===
    # nv: 变量数, ne: 等式数, nb: box bounds (不用), ng: 一般线性 (不用), nq: 二次约束数
    dim = hpipm_dense_qcqp_dim()
    dim.set('nv', n)    # 36
    dim.set('ne', ne)   # 18
    dim.set('nb', 0)    # 无 box 约束
    dim.set('ng', 0)    # 无一般线性约束
    dim.set('nq', nq)   # 72

    # === HPIPM QCQP 数据设置 ===
    # 对应 HPIPM 论文中的 dense QP 数据结构 (Section II-B)
    qcqp = hpipm_dense_qcqp(dim)
    qcqp.set('H', H)    # 目标 Hessian
    qcqp.set('g', g)    # 目标线性项

    if ne > 0:
        qcqp.set('A', A)  # 等式约束矩阵
        qcqp.set('b', b)  # 等式约束右端

    if nq > 0:
        # HPIPM 要求 Hq 水平堆叠为 (n, n*nq), gq 列堆叠为 (n, nq)
        # 对应 HPIPM 论文中的多约束存储格式
        Hq_stacked = np.hstack([np.asarray(hq, dtype=np.float64) for hq in Hq_list])
        gq_stacked = np.column_stack(gq_list)
        qcqp.set('Hq', Hq_stacked)  # 二次约束 Hessian (n × n*nq)
        qcqp.set('gq', gq_stacked)  # 二次约束线性项 (n × nq)
        qcqp.set('uq', uq_arr)      # 二次约束右端 (nq,)

    # === HPIPM 求解 (对应 RSS26 论文 Alg.1 line 5) ===
    qcqp_sol = hpipm_dense_qcqp_sol(dim)

    # 求解器参数: balance 模式 + tol=1e-6 + iter_max=300
    # HPIPM 论文 Section III: 'balance' 模式平衡速度与鲁棒性
    arg = hpipm_dense_qcqp_solver_arg(dim, 'balance')
    arg.set('iter_max', 300)
    arg.set('tol_stat', 1e-6)   # 状态可行性容差
    arg.set('tol_eq', 1e-6)     # 等式约束容差
    arg.set('tol_ineq', 1e-6)   # 不等式约束容差
    arg.set('tol_comp', 1e-6)   # 互补性容差

    # 创建求解器并求解
    solver = hpipm_dense_qcqp_solver(dim, arg)
    t0 = time.time()
    solver.solve(qcqp, qcqp_sol)
    solve_time = time.time() - t0

    # 提取解 (HPIPM 中决策变量名为 'v', 对应论文的 u ∈ R^{3×K})
    x = qcqp_sol.get('v').flatten()
    status = int(solver.get('status'))  # 0=SUCCESS (HPIPM 论文 Table I)
    iters = int(solver.get('iter')) if hasattr(solver, 'get') else 0
    obj_value = float(0.5 * x @ H @ x + g @ x)

    # balance 失败时回退到 robust 模式
    # HPIPM 论文 Section III: 'robust' 模式更稳定但更慢
    if status != 0:
        if verbose:
            print(f"[hpipm] balance 失败 (status={status}), 重试 robust...", file=sys.stderr)
        arg2 = hpipm_dense_qcqp_solver_arg(dim, 'robust')
        arg2.set('iter_max', 500)  # 增加最大迭代数
        arg2.set('tol_stat', 1e-6)
        arg2.set('tol_eq', 1e-6)
        arg2.set('tol_ineq', 1e-6)
        arg2.set('tol_comp', 1e-6)
        solver2 = hpipm_dense_qcqp_solver(dim, arg2)
        t1 = time.time()
        solver2.solve(qcqp, qcqp_sol)
        solve_time += time.time() - t1
        x = qcqp_sol.get('v').flatten()
        status = int(solver2.get('status'))
        obj_value = float(0.5 * x @ H @ x + g @ x)

    status_str = 'Solved' if status == 0 else f'Failed({status})'

    return {
        'x': x,                # 最优解 (36维: [u(:); nu(:)])
        'status': status,      # 0=成功
        'status_str': status_str,
        'obj_value': obj_value,
        'solve_time': float(solve_time),
        'iters': iters,
    }


def solve_ocp_qcqp(A, B, Q_eff, R_eff, b_stack, r_stack,
                   nx_arr, nu_arr, nq_arr, nbx_arr, nq_per_stage,
                   Qq_stack, Sq_stack, Rq_stack, qq_stack, rq_stack, uq_stack,
                   x0, idxbx, const=0.0, verbose=False):
    """
    求解 OCP QCQP (HPIPM ocp_qcqp 接口, 逐阶段结构).

    与 solve_qcqp (dense QCQP) 数学等价, 但利用 OCP 块三对角结构, 更高效.

    HPIPM OCP QCQP 标准形式 (每 stage n=0..N):
        动力学:  x_{n+1} = A_n x_n + B_n u_n + b_n      (n=0..N-1)
        代价:    min Σ_{n=0}^{N-1} [0.5 x'Qx + x'S'u + 0.5 u'Ru + q'x + r'u] + 0.5 x_N'Q_N x_N + q_N'x_N
        二次约束(每 stage nq 条): 0.5 v'Qq v + v'Sq'u + 0.5 u'Rq u + qq'v + rq'u <= uq
        box 约束: lbx <= x[idxbx] <= ubx

    ×2 约定 (与 construct_ocp_qp_from_rss.m 一致):
        Q_eff, R_eff, r, Qq, Sq, Rq, qq, rq 构造时已×2; HPIPM 的 0.5 前缀使其还原为原始系数.
        uq, const 不×2.

    参数 (均可为 numpy 数组或 MATLAB py.numpy.array 传入):
        A, B        : 动力学矩阵 (常数, 所有 stage 相同)
                      A: (nx, nx), B: (nx, nu)
        Q_eff       : (nx, nx) 状态代价 Hessian (已×2, 所有 stage 相同, 含终端)
        R_eff       : (nu, nu) 控制代价 Hessian (已×2, 含 R+rho*I)
        b_stack     : (nx, K) 动力学 bias, 每列一个 stage (n=0..K-1)
        r_stack     : (nu, K) 控制线性项, 每列一个 stage (已×2)
        nx_arr      : (N_stages,) int, 每 stage 状态维数
        nu_arr      : (N_stages,) int, 每 stage 控制维数
        nq_arr      : (N_stages,) int, 每 stage 二次约束数 (dim 用, = nq_per_stage)
        nbx_arr     : (N_stages,) int, 每 stage box 约束数
        nq_per_stage: (N_stages,) int, 每 stage 实际非空二次约束数
        Qq_stack    : (nx, nx, total_qcqp) 3D, 按 stage 顺序排列所有约束的状态二次项
        Sq_stack    : (nx, nu, total_qcqp) 3D, 状态-控制交叉 (MATLAB 约定 nx×nu, 内部转置为 nu×nx)
        Rq_stack    : (nu, nu, total_qcqp) 3D, 控制二次项
        qq_stack    : (nx, total_qcqp) 2D, 状态线性项
        rq_stack    : (nu, total_qcqp) 2D, 控制线性项
        uq_stack    : (total_qcqp,) 1D, 约束上界 (不×2)
        x0          : (nx,) 初始状态 (box bounds: lbx=ubx=x0)
        idxbx       : (nbx[0],) int, box 约束索引 (0-indexed, 索引 x 内位置)
        const       : float, 常数项 (不×2, 仅用于 obj 比较)

    返回 dict (与 solve_qcqp 兼容, 便于 control_RSS.m 复用):
        x          : (36,) [u(:); nu(:)], u 为控制序列, nu 为速度序列
        status     : int      0=成功
        status_str : str
        obj_value  : float    目标函数值
        solve_time : float    求解耗时 (秒)
        iters      : int      迭代次数
    """
    if not _HPIPM_OCP_OK:
        raise RuntimeError(
            f"HPIPM OCP QCQP 不可用: {_HPIPM_OCP_ERR}\n\n"
            "需要 HPIPM 编译时启用 OCP QCQP 支持 (d_ocp_qcqp_* 符号)."
        )

    # === 统一转为 numpy float64 / int32 ===
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    Q_eff = np.asarray(Q_eff, dtype=np.float64)
    R_eff = np.asarray(R_eff, dtype=np.float64)
    b_stack = np.asarray(b_stack, dtype=np.float64)
    r_stack = np.asarray(r_stack, dtype=np.float64)
    nx_arr = np.asarray(nx_arr, dtype=np.int32).flatten()
    nu_arr = np.asarray(nu_arr, dtype=np.int32).flatten()
    nq_arr = np.asarray(nq_arr, dtype=np.int32).flatten()
    nbx_arr = np.asarray(nbx_arr, dtype=np.int32).flatten()
    nq_per_stage = np.asarray(nq_per_stage, dtype=np.int32).flatten()
    Qq_stack = np.asarray(Qq_stack, dtype=np.float64)
    Sq_stack = np.asarray(Sq_stack, dtype=np.float64)
    Rq_stack = np.asarray(Rq_stack, dtype=np.float64)
    qq_stack = np.asarray(qq_stack, dtype=np.float64)
    rq_stack = np.asarray(rq_stack, dtype=np.float64)
    uq_stack = np.asarray(uq_stack, dtype=np.float64).flatten()
    x0 = np.asarray(x0, dtype=np.float64).flatten()
    idxbx = np.asarray(idxbx, dtype=np.int32).flatten()

    K = b_stack.shape[1]                  # horizon (=6)
    N_stages = K + 1                      # n=0..K
    nx0 = int(nx_arr[0])                  # 6
    nu0 = int(nu_arr[0])                  # 3

    # === HPIPM 维度设置 ===
    dim = hpipm_ocp_qcqp_dim(K)           # N=K (horizon)
    for s in range(N_stages):
        dim.set('nx', int(nx_arr[s]), s)
        dim.set('nu', int(nu_arr[s]), s)
        dim.set('nbx', int(nbx_arr[s]), s)
        dim.set('nq', int(nq_arr[s]), s)
        # ng=0, ns=0 (默认, 无一般线性/软约束)

    # === HPIPM OCP QCQP 数据设置 (逐 stage) ===
    qp = hpipm_ocp_qcqp(dim)

    # 动力学 (n=0..K-1): A, B, b
    for n in range(K):
        qp.set('A', A, n)
        qp.set('B', B, n)
        qp.set('b', b_stack[:, n], n)

    # 代价 (n=0..K): Q, q=0; (n=0..K-1): R, r, S=0; 终端 n=K: 仅 Q, q
    q_zero = np.zeros(nx0)
    S_zero = np.zeros((nu0, nx0))         # HPIPM S 是 (nu, nx)
    for n in range(N_stages):
        qp.set('Q', Q_eff, n)
        qp.set('q', q_zero, n)
        if n < K:
            qp.set('R', R_eff, n)
            qp.set('r', r_stack[:, n], n)
            qp.set('S', S_zero, n)

    # 二次约束 (逐 stage 堆叠后 set)
    offset = 0
    for s in range(N_stages):
        nq_s = int(nq_per_stage[s])
        if nq_s > 0:
            nx_s = int(nx_arr[s])
            nu_s = int(nu_arr[s])
            # Qq: 水平堆叠 (nx, nx*nq_s) — 维度仅依赖 nx_s, 不受 nu_s 影响
            Qq_s = np.hstack([Qq_stack[:, :, offset + j] for j in range(nq_s)])
            # qq: 列堆叠 (nx, nq_s) — 同上
            qq_s = np.column_stack([qq_stack[:, offset + j] for j in range(nq_s)])
            # uq: (nq_s,) 向量
            uq_s = uq_stack[offset:offset + nq_s]
            qp.set('Qq', Qq_s, s)
            qp.set('qq', qq_s, s)
            qp.set('uq', uq_s, s)
            # Sq/Rq/rq 维度依赖 nu_s: 终端 stage nu_s=0 时必须传 0 行数组,
            # 否则 (3, ...) 与 HPIPM 期望的 (0, ...) 不匹配 → 堆缓冲区溢出 → NaN/Inf
            if nu_s > 0:
                # Sq: HPIPM C 层 OCP_QCQP_SET_SQ 用 CVT_TRAN_MAT2STRMAT(=blasfeo_pack_tran_dmat)
                #     读取列优先 (nu, nx) 矩阵, 转置后存入 Hq 下左块 (nx, nu)
                # Python wrapper 用 np.ravel('F') 列优先展开
                # 需传入 (nu, nx*nq_s) 数组, 使 ravel('F') = 列优先 (nu, nx*nq_s)
                #     = 每约束 nu*nx 元素的列优先 (nu, nx) 布局 = C 期望格式
                # Sq_stack 是 (nx, nu, nq) [MATLAB (nx,nu) 约定], 需先转置为 (nu, nx, nq) 再 reshape
                Sq_s = np.transpose(Sq_stack[:, :, offset:offset + nq_s], (1, 0, 2)).reshape((nu_s, nx_s * nq_s), order='F')
                # Rq: (nu, nu) 对称, 'F' ravel 的 nu 元素组 = C 行优先 nu 元素组 (对称阵等价)
                Rq_s = np.hstack([Rq_stack[:, :, offset + j] for j in range(nq_s)])
                # rq: (nu, nq_s) 列堆叠, 'F' ravel = [c0(nu), c1(nu), ...] = C 期望
                rq_s = np.column_stack([rq_stack[:, offset + j] for j in range(nq_s)])
                qp.set('Sq', Sq_s, s)
                qp.set('Rq', Rq_s, s)
                qp.set('rq', rq_s, s)
            else:
                # 终端 stage (nu=0): 传 0 行数组, 维度与 HPIPM 期望 (nu=0) 一致
                qp.set('Sq', np.zeros((0, nx_s * nq_s)), s)
                qp.set('Rq', np.zeros((0, 0)), s)
                qp.set('rq', np.zeros((0, nq_s)), s)
        offset += nq_s

    # 初始状态 box 约束 (stage 0: lbx=ubx=x0, idxbx=0..5)
    if int(nbx_arr[0]) > 0:
        qp.set('idxbx', idxbx, 0)
        qp.set('lbx', x0, 0)
        qp.set('ubx', x0, 0)

    # === 诊断: 打印维度信息 (临时, 用于排查 status=3) ===
    if verbose:
        print(f"\n[hpipm_diag] K={K}, N_stages={N_stages}", file=sys.stderr)
        print(f"[hpipm_diag] nx_arr={nx_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_diag] nu_arr={nu_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_diag] nq_arr={nq_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_diag] nbx_arr={nbx_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_diag] nq_per_stage={nq_per_stage.tolist()}", file=sys.stderr)
        print(f"[hpipm_diag] total_qcqp={int(nq_per_stage.sum())}", file=sys.stderr)
        # 检查数据有限性
        for name, arr in [('A', A), ('B', B), ('Q_eff', Q_eff), ('R_eff', R_eff),
                          ('b_stack', b_stack), ('r_stack', r_stack),
                          ('x0', x0), ('Qq_stack', Qq_stack), ('Sq_stack', Sq_stack),
                          ('Rq_stack', Rq_stack), ('qq_stack', qq_stack),
                          ('rq_stack', rq_stack), ('uq_stack', uq_stack)]:
            if not np.all(np.isfinite(arr)):
                print(f"[hpipm_diag] !!! {name} contains NaN/Inf !!!", file=sys.stderr)
        print(f"[hpipm_diag] uq_stack={uq_stack}", file=sys.stderr)
        # 检查每个约束在 x=0, u=0 处的值: 应为 0 <= uq (uq>0 时可行)
        # 约束形式: 0.5*x'Qq*x + ... + qq'x + rq'u <= uq
        # 在 x=0, u=0 处: 0 <= uq, 所以 uq 必须 >= 0
        neg_uq = np.where(uq_stack < 0)[0]
        if len(neg_uq) > 0:
            print(f"[hpipm_diag] !!! uq_stack 有 {len(neg_uq)} 个负值: 索引 {neg_uq.tolist()}", file=sys.stderr)
        else:
            print(f"[hpipm_diag] uq_stack 全部 >= 0 (约束在原点可行)", file=sys.stderr)
        # 手动验证 stage 0 约束 0 在 x=x0, u=0 处的约束值
        # 约束形式: 0.5*x'Qq*x + x'Sq'u + 0.5*u'Rq*u + qq'x + rq'u <= uq
        # stage 0: x=x0 (固定), u=0
        s = 0
        nq_s0 = int(nq_per_stage[s])
        nx_s0 = int(nx_arr[s])
        nu_s0 = int(nu_arr[s])
        if nq_s0 > 0:
            print(f"[hpipm_diag] --- stage 0 约束验证 (x=x0, u=0) ---", file=sys.stderr)
            for j in range(min(nq_s0, 4)):  # 只打印前 4 条
                Qq_j = Qq_stack[:, :, j]
                qq_j = qq_stack[:, j]
                uq_j = uq_stack[j]
                # x=x0, u=0
                val_x = 0.5 * x0 @ Qq_j @ x0 + qq_j @ x0
                print(f"[hpipm_diag]   约束 {j}: 0.5*x0'Qq*x0 + qq'x0 = {val_x:.6e}, uq = {uq_j:.6e}, 违反量 = {val_x - uq_j:.6e}", file=sys.stderr)
        # 验证 stage 1 (有 u) 约束 0 在 x=0, u=0
        s = 1
        nq_s1 = int(nq_per_stage[s])
        if nq_s1 > 0:
            print(f"[hpipm_diag] --- stage 1 约束验证 (x=0, u=0) ---", file=sys.stderr)
            offset1 = int(nq_per_stage[0])
            for j in range(min(nq_s1, 4)):
                Qq_j = Qq_stack[:, :, offset1 + j]
                Rq_j = Rq_stack[:, :, offset1 + j]
                uq_j = uq_stack[offset1 + j]
                # x=0, u=0
                val = 0.0  # 所有项为 0
                print(f"[hpipm_diag]   约束 {j}: val(0,0)=0, uq = {uq_j:.6e}, Qq PSD eigvals = {np.linalg.eigvalsh(Qq_j).min():.3e}/{np.linalg.eigvalsh(Qq_j).max():.3e}, Rq PSD eigvals = {np.linalg.eigvalsh(Rq_j).min():.3e}/{np.linalg.eigvalsh(Rq_j).max():.3e}", file=sys.stderr)
        # 打印 HPIPM 内部结构到文件 (避免被 MATLAB 吞掉)
        try:
            with open('hpipm_struct_dump.txt', 'w') as f:
                import io
                from contextlib import redirect_stdout
                buf = io.StringIO()
                with redirect_stdout(buf):
                    qp.print_C_struct()
                f.write(buf.getvalue())
            print(f"[hpipm_diag] HPIPM C struct 已写入 hpipm_struct_dump.txt", file=sys.stderr)
        except Exception as e:
            print(f"[hpipm_diag] print_C_struct 失败: {e}", file=sys.stderr)

    # === 求解 ===
    qp_sol = hpipm_ocp_qcqp_sol(dim)

    # balance 模式 + tol=1e-6 + iter_max=300
    arg = hpipm_ocp_qcqp_solver_arg(dim, 'balance')
    arg.set('iter_max', 300)
    arg.set('tol_stat', 1e-6)
    arg.set('tol_eq', 1e-6)
    arg.set('tol_ineq', 1e-6)
    arg.set('tol_comp', 1e-6)

    solver = hpipm_ocp_qcqp_solver(dim, arg)
    t0 = time.time()
    solver.solve(qp, qp_sol)
    solve_time = time.time() - t0

    status = int(solver.get('status'))
    iters = int(solver.get('iter'))
    if verbose:
        print(f"[hpipm_diag] status={status}, iters={iters}", file=sys.stderr)

    # 提取解: u_n (n=0..K-1), x_n (n=0..K)
    u_list = []
    x_list = []
    for n in range(N_stages):
        x_n = np.asarray(qp_sol.get('x', n)).flatten()
        x_list.append(x_n)
        if n < K:
            u_n = np.asarray(qp_sol.get('u', n)).flatten()
            u_list.append(u_n)

    # 拼成 36 维 [u(:); nu(:)] 兼容 control_RSS.m 的提取逻辑
    # u_n (n=0..K-1) = 论文 u_{n+1} (k=1..K), 列顺序对齐
    u_flat = np.concatenate(u_list) if u_list else np.zeros(0)   # 18 维
    # nu (速度序列): x_k 后 3 维 (k=1..K), 即 v_k
    nu_flat = np.concatenate([x_list[k][3:6] for k in range(1, N_stages)])  # 18 维
    x_out = np.concatenate([u_flat, nu_flat])

    # 手动计算 obj_value (Q_eff, R_eff, r 已×2; 0.5 前缀还原)
    obj_value = const
    for n in range(N_stages):
        xn = x_list[n]
        obj_value += 0.5 * float(xn @ Q_eff @ xn)     # 状态代价 (已×2 → 0.5*2=1)
        if n < K:
            un = u_list[n]
            obj_value += 0.5 * float(un @ R_eff @ un)  # 控制代价
            obj_value += float(r_stack[:, n] @ un)      # 控制线性项 (已×2)

    # balance 失败回退 robust
    if status != 0:
        if verbose:
            print(f"[hpipm_ocp] balance 失败 (status={status}), 重试 robust...", file=sys.stderr)
        arg2 = hpipm_ocp_qcqp_solver_arg(dim, 'robust')
        arg2.set('iter_max', 500)
        arg2.set('tol_stat', 1e-6)
        arg2.set('tol_eq', 1e-6)
        arg2.set('tol_ineq', 1e-6)
        arg2.set('tol_comp', 1e-6)
        solver2 = hpipm_ocp_qcqp_solver(dim, arg2)
        t1 = time.time()
        solver2.solve(qp, qp_sol)
        solve_time += time.time() - t1
        status = int(solver2.get('status'))
        iters = int(solver2.get('iter'))
        # 重新提取解
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
            obj_value += 0.5 * float(xn @ Q_eff @ xn)
            if n < K:
                un = u_list[n]
                obj_value += 0.5 * float(un @ R_eff @ un)
                obj_value += float(r_stack[:, n] @ un)

    status_str = 'Solved' if status == 0 else f'Failed({status})'

    return {
        'x': x_out,                 # (36,) [u(:); nu(:)]
        'status': status,           # 0=成功
        'status_str': status_str,
        'obj_value': float(obj_value),
        'solve_time': float(solve_time),
        'iters': iters,
    }


def solve_ocp_qp(A, B, Q_eff, R_eff, b_stack, r_stack,
                 nx_arr, nu_arr, ng_arr, nbx_arr, ng_per_stage,
                 Cmat_stack, Dmat_stack, lg_stack, ug_stack,
                 x0, idxbx, const=0.0, verbose=False):
    """
    求解 OCP QP (HPIPM ocp_qp 接口, 一般线性约束).

    用于替代 solve_ocp_qcqp — HPIPM OCP QCQP IPM solver (d_ocp_qcqp_ipm_solve)
    存在 bug 导致 status=3 (NAN_SOL), 而 OCP QP 接口工作正常.
    二次约束在 MATLAB 端 (construct_ocp_qp_from_rss.m) 做一阶泰勒线性化,
    转为一般线性约束 lg <= C*x + D*u <= ug.
    外层 SCP 迭代 (control_RSS.m 中 m=1..3) 保证收敛到原二次约束的解.

    HPIPM OCP QP 标准形式 (每 stage n=0..N):
        动力学:  x_{n+1} = A_n x_n + B_n u_n + b_n      (n=0..N-1)
        代价:    min Σ [0.5 x'Qx + 0.5 u'Ru + r'u] + 0.5 x_N'Q_N x_N
        一般线性约束 (每 stage ng 条): lg <= C*x + D*u <= ug
        box 约束: lbx <= x[idxbx] <= ubx

    ×2 约定 (与 construct_ocp_qp_from_rss.m 一致):
        Q_eff, R_eff, r 构造时已×2; HPIPM 的 0.5 前缀使其还原为原始系数.
        C, D, lg, ug, const 不×2 (线性约束和常数项无 1/2 前缀).

    参数 (均可为 numpy 数组或 MATLAB py.numpy.array 传入):
        A, B        : 动力学矩阵 (常数, 所有 stage 相同)
                      A: (nx, nx), B: (nx, nu)
        Q_eff       : (nx, nx) 状态代价 Hessian (已×2, 所有 stage 相同, 含终端)
        R_eff       : (nu, nu) 控制代价 Hessian (已×2, 含 R+rho*I)
        b_stack     : (nx, K) 动力学 bias, 每列一个 stage (n=0..K-1)
        r_stack     : (nu, K) 控制线性项, 每列一个 stage (已×2)
        nx_arr      : (N_stages,) int, 每 stage 状态维数
        nu_arr      : (N_stages,) int, 每 stage 控制维数
        ng_arr      : (N_stages,) int, 每 stage 一般线性约束数 (dim 用)
        nbx_arr     : (N_stages,) int, 每 stage box 约束数
        ng_per_stage: (N_stages,) int, 每 stage 实际非空线性约束数
        Cmat_stack  : (total_ng, nx) 所有线性约束的状态系数, 按 stage 顺序行堆叠
                      每行是一条约束的 C 系数 (1×nx), HPIPM set('C') 期望 (ng, nx)
        Dmat_stack  : (total_ng, nu) 所有线性约束的控制系数, 按 stage 顺序行堆叠
                      每行是一条约束的 D 系数 (1×nu), HPIPM set('D') 期望 (ng, nu)
        lg_stack    : (total_ng,) 下界 (按 stage 顺序)
        ug_stack    : (total_ng,) 上界 (按 stage 顺序)
        x0          : (nx,) 初始状态 (box bounds: lbx=ubx=x0)
        idxbx       : (nbx[0],) int, box 约束索引 (0-indexed, 索引 x 内位置)
        const       : float, 常数项 (不×2, 仅用于 obj 比较)
        verbose     : bool, 是否打印诊断信息

    返回 dict (与 solve_ocp_qcqp 兼容, 便于 control_RSS.m 复用):
        x          : (36,) [u(:); nu(:)], u 为控制序列, nu 为速度序列
        status     : int      0=成功
        status_str : str
        obj_value  : float    目标函数值
        solve_time : float    求解耗时 (秒)
        iters      : int      迭代次数
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
        # nbu=0, ns=0 (默认, 无控制 box/软约束)

    # === HPIPM OCP QP 数据设置 (逐 stage) ===
    qp = hpipm_ocp_qp(dim)

    # 动力学 (n=0..K-1): A, B, b
    for n in range(K):
        qp.set('A', A, n)
        qp.set('B', B, n)
        qp.set('b', b_stack[:, n], n)

    # 代价 (n=0..K): Q, q=0; (n=0..K-1): R, r, S=0; 终端 n=K: 仅 Q, q
    q_zero = np.zeros(nx0)
    S_zero = np.zeros((nu0, nx0))         # HPIPM S 是 (nu, nx)
    for n in range(N_stages):
        qp.set('Q', Q_eff, n)
        qp.set('q', q_zero, n)
        if n < K:
            qp.set('R', R_eff, n)
            qp.set('r', r_stack[:, n], n)
            qp.set('S', S_zero, n)

    # 一般线性约束 (逐 stage, 按 ng_per_stage 切片 Cmat_stack/Dmat_stack)
    # HPIPM OCP QP: set('C', (ng, nx), stage), set('D', (ng, nu), stage)
    #               set('lg', (ng,), stage), set('ug', (ng,), stage)
    # 约束形式: lg <= C*x + D*u <= ug
    offset = 0
    for s in range(N_stages):
        ng_s = int(ng_per_stage[s])
        if ng_s > 0:
            nx_s = int(nx_arr[s])
            nu_s = int(nu_arr[s])
            # C: (ng_s, nx_s) — 从 Cmat_stack 切片行
            C_s = Cmat_stack[offset:offset + ng_s, :nx_s]
            qp.set('C', C_s, s)
            # D: (ng_s, nu_s) — 终端 stage nu_s=0 时跳过 (HPIPM 自动处理)
            if nu_s > 0:
                D_s = Dmat_stack[offset:offset + ng_s, :nu_s]
                qp.set('D', D_s, s)
            # lg/ug: (ng_s,) 向量
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

    # === 诊断: 打印维度信息 (仅首步首迭代) ===
    if verbose:
        print(f"\n[hpipm_ocp_qp_diag] K={K}, N_stages={N_stages}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] nx_arr={nx_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] nu_arr={nu_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] ng_arr={ng_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] nbx_arr={nbx_arr.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] ng_per_stage={ng_per_stage.tolist()}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] total_ng={int(ng_per_stage.sum())}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] Cmat_stack.shape={Cmat_stack.shape}", file=sys.stderr)
        print(f"[hpipm_ocp_qp_diag] Dmat_stack.shape={Dmat_stack.shape}", file=sys.stderr)
        # 检查数据有限性
        for name, arr in [('A', A), ('B', B), ('Q_eff', Q_eff), ('R_eff', R_eff),
                          ('b_stack', b_stack), ('r_stack', r_stack),
                          ('x0', x0), ('Cmat_stack', Cmat_stack),
                          ('Dmat_stack', Dmat_stack), ('lg_stack', lg_stack),
                          ('ug_stack', ug_stack)]:
            if not np.all(np.isfinite(arr)):
                print(f"[hpipm_ocp_qp_diag] !!! {name} contains NaN/Inf !!!", file=sys.stderr)
        # 打印 HPIPM 内部结构到文件
        try:
            with open('hpipm_ocp_qp_struct_dump.txt', 'w') as f:
                import io
                from contextlib import redirect_stdout
                buf = io.StringIO()
                with redirect_stdout(buf):
                    qp.print_C_struct()
                f.write(buf.getvalue())
            print(f"[hpipm_ocp_qp_diag] HPIPM C struct 已写入 hpipm_ocp_qp_struct_dump.txt", file=sys.stderr)
        except Exception as e:
            print(f"[hpipm_ocp_qp_diag] print_C_struct 失败: {e}", file=sys.stderr)

    # === 求解 ===
    qp_sol = hpipm_ocp_qp_sol(dim)

    # balance 模式 + tol=1e-6 + iter_max=300
    arg = hpipm_ocp_qp_solver_arg(dim, 'balance')
    arg.set('iter_max', 300)
    arg.set('tol_stat', 1e-6)
    arg.set('tol_eq', 1e-6)
    arg.set('tol_ineq', 1e-6)
    arg.set('tol_comp', 1e-6)

    solver = hpipm_ocp_qp_solver(dim, arg)
    t0 = time.time()
    solver.solve(qp, qp_sol)
    solve_time = time.time() - t0

    status = int(solver.get('status'))
    iters = int(solver.get('iter'))
    if verbose:
        print(f"[hpipm_ocp_qp_diag] status={status}, iters={iters}", file=sys.stderr)

    # 提取解: u_n (n=0..K-1), x_n (n=0..K)
    u_list = []
    x_list = []
    for n in range(N_stages):
        x_n = np.asarray(qp_sol.get('x', n)).flatten()
        x_list.append(x_n)
        if n < K:
            u_n = np.asarray(qp_sol.get('u', n)).flatten()
            u_list.append(u_n)

    # 拼成 36 维 [u(:); nu(:)] 兼容 control_RSS.m 的提取逻辑
    # u_n (n=0..K-1) = 论文 u_{n+1} (k=1..K), 列顺序对齐
    u_flat = np.concatenate(u_list) if u_list else np.zeros(0)   # 18 维
    # nu (速度序列): x_k 后 3 维 (k=1..K), 即 v_k
    nu_flat = np.concatenate([x_list[k][3:6] for k in range(1, N_stages)])  # 18 维
    x_out = np.concatenate([u_flat, nu_flat])

    # 手动计算 obj_value (Q_eff, R_eff, r 已×2; 0.5 前缀还原)
    obj_value = const
    for n in range(N_stages):
        xn = x_list[n]
        obj_value += 0.5 * float(xn @ Q_eff @ xn)     # 状态代价 (已×2 → 0.5*2=1)
        if n < K:
            un = u_list[n]
            obj_value += 0.5 * float(un @ R_eff @ un)  # 控制代价
            obj_value += float(r_stack[:, n] @ un)      # 控制线性项 (已×2)

    # balance 失败回退 robust
    if status != 0:
        if verbose:
            print(f"[hpipm_ocp_qp] balance 失败 (status={status}), 重试 robust...", file=sys.stderr)
        arg2 = hpipm_ocp_qp_solver_arg(dim, 'robust')
        arg2.set('iter_max', 500)
        arg2.set('tol_stat', 1e-6)
        arg2.set('tol_eq', 1e-6)
        arg2.set('tol_ineq', 1e-6)
        arg2.set('tol_comp', 1e-6)
        solver2 = hpipm_ocp_qp_solver(dim, arg2)
        t1 = time.time()
        solver2.solve(qp, qp_sol)
        solve_time += time.time() - t1
        status = int(solver2.get('status'))
        iters = int(solver2.get('iter'))
        # 重新提取解
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
            obj_value += 0.5 * float(xn @ Q_eff @ xn)
            if n < K:
                un = u_list[n]
                obj_value += 0.5 * float(un @ R_eff @ un)
                obj_value += float(r_stack[:, n] @ un)

    status_str = 'Solved' if status == 0 else f'Failed({status})'

    return {
        'x': x_out,                 # (36,) [u(:); nu(:)]
        'status': status,           # 0=成功
        'status_str': status_str,
        'obj_value': float(obj_value),
        'solve_time': float(solve_time),
        'iters': iters,
    }
