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
