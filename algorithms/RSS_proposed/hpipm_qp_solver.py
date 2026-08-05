"""
hpipm_qp_solver.py
通用 QCQP 求解器 (HPIPM dense QCQP)

仅负责求解, 不构造 QP. QP 矩阵由 MATLAB 端 (construct_complete_qp_from_rss.m) 构造后传入.

接口:
    solve_qcqp(H, g, A, b, Hq, gq, uq) -> dict(x, status, obj_value, solve_time, iters)

QP 形式:
    min  0.5*x'Hx + g'x
    s.t. A*x = b
         0.5*x'Hq_i*x + gq_i'*x <= uq_i  (i = 1..nq)
"""

import sys
import time
import numpy as np

# HPIPM Python 接口
_hpipm_path = None
import os
import sys

_candidate = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'third_party', 'hpipm', 'interfaces', 'python', 'hpipm_python')
if not os.path.isdir(_candidate):
    _candidate = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'third_party', 'hpipm', 'interfaces', 'python')
if os.path.isdir(_candidate):
    if _candidate not in sys.path:
        sys.path.insert(0, _candidate)
    _hpipm_path = _candidate

# Windows: 把 libhpipm.dll 所在目录加入 DLL 搜索路径
# (hpipm Python wrapper 内部用 CDLL('libhpipm.dll') 加载, 需保证目录在搜索路径中)
_hpipm_lib_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'third_party', 'hpipm', 'lib')
_hpipm_lib_dir = os.path.abspath(_hpipm_lib_dir)
if os.path.isdir(_hpipm_lib_dir):
    if sys.platform.startswith('win'):
        # Python 3.8+: 使用 add_dll_directory
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


def solve_qcqp(H, g, A, b, Hq, gq, uq, verbose=False):
    """
    求解 dense QCQP.

    参数 (均可为 numpy 数组或 MATLAB py.numpy.array 传入):
        H  : (n, n)   目标函数 Hessian (PSD)
        g  : (n,)     目标函数线性项
        A  : (ne, n)  等式约束矩阵 (无则 None/空)
        b  : (ne,)    等式约束右端
        Hq : list of (n, n)  二次约束 Hessian 列表
        gq : list of (n,)    二次约束线性项列表
        uq : (nq,)           二次约束右端

    返回 dict:
        x          : (n,)     最优解
        status     : int      0=成功
        status_str : str
        obj_value  : float
        solve_time : float
        iters      : int
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

    # === 统一转为 numpy 数组 ===
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

    # Hq / gq / uq
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

    # === 维度 ===
    dim = hpipm_dense_qcqp_dim()
    dim.set('nv', n)
    dim.set('ne', ne)
    dim.set('nb', 0)
    dim.set('ng', 0)
    dim.set('nq', nq)

    # === QCQP 数据 ===
    qcqp = hpipm_dense_qcqp(dim)
    qcqp.set('H', H)
    qcqp.set('g', g)

    if ne > 0:
        qcqp.set('A', A)
        qcqp.set('b', b)

    if nq > 0:
        # 关键: Hq 必须水平堆叠为 (n, n*nq), gq 必须列堆叠为 (n, nq)
        Hq_stacked = np.hstack([np.asarray(hq, dtype=np.float64) for hq in Hq_list])
        gq_stacked = np.column_stack(gq_list)
        qcqp.set('Hq', Hq_stacked)
        qcqp.set('gq', gq_stacked)
        qcqp.set('uq', uq_arr)

    # === 求解 ===
    qcqp_sol = hpipm_dense_qcqp_sol(dim)

    # 与参考仓库一致: balance 模式 + tol=1e-6 + iter_max=300
    arg = hpipm_dense_qcqp_solver_arg(dim, 'balance')
    arg.set('iter_max', 300)
    arg.set('tol_stat', 1e-6)
    arg.set('tol_eq', 1e-6)
    arg.set('tol_ineq', 1e-6)
    arg.set('tol_comp', 1e-6)

    solver = hpipm_dense_qcqp_solver(dim, arg)
    t0 = time.time()
    solver.solve(qcqp, qcqp_sol)
    solve_time = time.time() - t0

    x = qcqp_sol.get('v').flatten()
    status = int(solver.get('status'))
    iters = int(solver.get('iter')) if hasattr(solver, 'get') else 0
    obj_value = float(0.5 * x @ H @ x + g @ x)

    # balance 失败时回退到 robust 模式 (与参考仓库一致)
    if status != 0:
        if verbose:
            print(f"[hpipm] balance 失败 (status={status}), 重试 robust...", file=sys.stderr)
        arg2 = hpipm_dense_qcqp_solver_arg(dim, 'robust')
        arg2.set('iter_max', 500)
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
        'x': x,
        'status': status,
        'status_str': status_str,
        'obj_value': obj_value,
        'solve_time': float(solve_time),
        'iters': iters,
    }
