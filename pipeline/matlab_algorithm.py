"""MATLAB Engine 算法桥: e-lmpc / interior-point / active-set 每步求解.

架构 (pipeline 统一框架, 4 算法同一条主干):
    trajectory_generator.py + params.py   (Python) 参考轨迹 + 场景/车辆参数
    simulator.py                          (Python) 原始动力学闭环, 每步调用算法
        |- proposed-3iter        -> algorithms/RSS_proposed/control_rss_ocpqcqp.py
        |- e-lmpc 等 3 个 MATLAB 算法 -> 本桥 (MATLAB 只做算法求解)
              control() -> MATLAB Engine 常驻会话
                        -> pipeline/matlab_control_bridge.m
                        -> algorithms/{RSS_sqp|RSS_fmincon|RSS_active_set}/control_RSS.m
    结果直接进 Python SimResult, metrics/plotting/三件套与纯 Python 路径完全统一.

MATLAB Engine 会话生命周期: __init__ 启动一次 (约 30s), close() 退出;
每步 control() 仅函数级调用 (毫秒级序列化开销), 轨迹等大数组经 base
workspace 一次性传入, 不逐步序列化.

依赖: matlabengine 包 (R2026a 支持 Python 3.9-3.13), MATLAB 在 PATH。
安装: cd <matlabroot>/extern/engines/python && python setup.py install
     (或 pip install matlabengine, 版本须与 MATLAB 发行版匹配)
"""

from pathlib import Path
import tempfile
import os

import numpy as np
from scipy.io import savemat

MATLAB_BRIDGE_ALGORITHMS = ('e-lmpc', 'interior-point', 'active-set')


def _sanitize_for_matlab(obj):
    """递归清洗 run_config 供 scipy savemat 序列化.

    savemat 不支持 None (如未传 --rho 时 to_dict() 的 rho_schedule=None),
    遇到会抛异常导致 PIPELINE_CONFIG 注入失败、MATLAB 侧退回旧架构兜底;
    None 值的键直接移除 (MATLAB 侧 isfield 检查自然走默认值),
    list/tuple 统一转 numpy 数组 (含嵌套 dict 内的 vehicle.wheel_pos).
    """
    if isinstance(obj, dict):
        return {k: _sanitize_for_matlab(v) for k, v in obj.items()
                if v is not None}
    if isinstance(obj, (list, tuple)):
        arr = np.asarray(obj)
        # numpy unicode 数组同样不被 savemat 支持: 保持原生 list (转 cell)
        return obj if arr.dtype.kind == 'U' else arr
    return obj


class MatlabAlgorithmBridge:
    """MATLAB 算法的每步求解桥 (供 simulator 循环调用)."""

    ALGO_SUBMODULE = {
        'e-lmpc': 'RSS_sqp',
        'interior-point': 'RSS_fmincon',
        'active-set': 'RSS_active_set',
    }

    def __init__(self, algorithm: str, trajectory, verbose: bool = True,
                 repo_root=None, K: float = None, run_config: dict | None = None):
        """启动 MATLAB Engine 并完成一次性初始化.

        trajectory: pipeline.trajectory_generator.Trajectory (提供 path/seed_id,
                    场景数据与 MATLAB scenario_bank/defaultConfig 同源).
        """
        import matlab.engine  # 延迟导入: 纯 Python 路径无需安装 matlabengine

        algorithm = algorithm.lower()
        if algorithm not in MATLAB_BRIDGE_ALGORITHMS:
            raise ValueError(
                f'MatlabAlgorithmBridge 仅服务 {MATLAB_BRIDGE_ALGORITHMS}, '
                f'收到: {algorithm}')

        self.algorithm = algorithm
        self.verbose = verbose
        self.eng = None
        self.run_config = run_config

        repo_root = Path(repo_root).resolve() if repo_root else \
            Path(__file__).resolve().parent.parent
        path = np.asarray(trajectory.path, dtype=np.float64)

        if verbose:
            print(f'[matlab_engine] 启动 MATLAB Engine (约 30s): '
                  f'{self.algorithm} ({self.ALGO_SUBMODULE[algorithm]}) ...')
        self.eng = matlab.engine.start_matlab()
        self.eng.addpath(str(repo_root / 'pipeline').replace('\\', '/'))
        # 初始化数据写入 base workspace (matlab_control_bridge.m 首次调用读取)
        self.eng.workspace['PIPELINE_ALGORITHM'] = self.algorithm
        self.eng.workspace['PIPELINE_SEED'] = float(trajectory.seed_id)
        self.eng.workspace['PIPELINE_PATH'] = matlab.double(path.tolist())
        # 将请求的预测时域 K 写入 base workspace 供 matlab_control_bridge 读取
        if K is not None:
            try:
                self.eng.workspace['PIPELINE_K'] = float(K)
            except Exception:
                # 保底：如果赋值失败，仍可由 MATLAB 端使用默认 K=6
                pass
        # 如果外部传入了完整 run_config，生成 .mat 并加载到 base workspace
        run_config = None
        try:
            # if caller passed run_config in kwargs (simulate does), retrieve it
            # Note: trajectory may carry much of the info; but prefer explicit run_config
            caller_locals = None
        except Exception:
            pass
        # attempt to find attribute 'run_config' on self if set externally
        if hasattr(self, 'run_config') and self.run_config is not None:
            run_config = self.run_config
        # If initializer received run_config via keyword args, Python won't set it automatically,
        # so check trajectory for algorithm_params to assemble a run_config fallback
        if run_config is None:
            try:
                alg_params = getattr(trajectory, 'algorithm_params', None)
                run_config = {
                    'seed_id': int(trajectory.seed_id),
                    'algorithm': self.algorithm,
                    'K': int(K) if K is not None else None,
                    'dt': float(getattr(trajectory, 'dt', 0.01)),
                    'num_steps': int(getattr(trajectory, 'num_steps', trajectory.num_points if hasattr(trajectory, 'num_points') else 0)),
                    'trajectory_source': getattr(trajectory, 'source', ''),
                    'scenario_name': getattr(trajectory, 'scenario_name', ''),
                    'algorithm_params': alg_params.to_dict() if alg_params is not None else {},
                }
            except Exception:
                run_config = None

        if run_config is not None:
            try:
                fd, tmp = tempfile.mkstemp(prefix='pipeline_run_config_', suffix='.mat')
                os.close(fd)
                # savemat 不支持 None: 递归剔除 (rho_schedule=None 等),
                # 否则注入失败后 MATLAB 侧会退回旧架构兜底 (core/ 已归档, 必崩)
                matdict = {'PIPELINE_CONFIG': _sanitize_for_matlab(run_config)}
                savemat(tmp, matdict, do_compression=False)
                try:
                    self.eng.load(tmp, nargout=0)
                except Exception:
                    # fallback to eval load if needed
                    tmp_fwd = tmp.replace('\\', '/')
                    self.eng.eval(f"load('{tmp_fwd}')", nargout=0)
            except Exception as e:
                if verbose:
                    print(f'[matlab_engine] 警告: 未能写入/加载 PIPELINE_CONFIG '
                          f'到 MATLAB ({e})')
        if verbose:
            print('[matlab_engine] 就绪 (每步函数级调用, '
                  '轨迹/闭环/评估均在 Python 侧)')

    def control(self, path, k, last_body_velocity, state, params=None,
                verbose: bool = True):
        """单步求解, 签名与 control_rss_ocpqcqp 统一.

        返回 (u_full(3,K), world_velocity(3,), body_velocity(3,), diagnostics).
        """
        if self.eng is None:
            raise RuntimeError('MATLAB Engine 会话已关闭')

        vel = np.asarray(last_body_velocity, dtype=np.float64).reshape(3, 1)
        st = np.asarray(state, dtype=np.float64).reshape(1, 3)
        # 逐步 rho 注入 (CLI --rho 序列): 第 k 步取 rho_k 写入 base workspace,
        # matlab_control_bridge 每步读取覆盖 config.rho; 未传 --rho 时不写,
        # MATLAB 侧保持各算法默认 rho (不改变基线行为)
        if params is not None and getattr(params, 'rho_schedule', None):
            rho_idx = min(int(k) - 1, len(params.rho_schedule) - 1)
            self.eng.workspace['PIPELINE_RHO'] = float(
                params.rho_schedule[rho_idx])
        world_velocity, body_velocity, solve_time, iter_num, u_col, log_text = \
            self.eng.matlab_control_bridge(
                float(k),
                matlab_double(vel),
                matlab_double(st),
                nargout=6)

        if verbose and log_text:
            print(str(log_text).rstrip())

        world_velocity = np.asarray(world_velocity, dtype=np.float64).flatten()
        body_velocity = np.asarray(body_velocity, dtype=np.float64).flatten()
        u = np.asarray(u_col, dtype=np.float64).flatten()

        K = int(params.K) if params is not None else 6
        u_full = np.zeros((3, K))
        u_full[:, 0] = u

        diagnostics = {
            'total_solve_time': float(solve_time),
            'solver_call_count': int(iter_num),   # fmincon 迭代数
            'matlab_algorithm': True,
        }
        return u_full, world_velocity, body_velocity, diagnostics

    def close(self):
        """退出 MATLAB Engine 会话 (幂等)."""
        if getattr(self, 'eng', None) is not None:
            try:
                self.eng.quit()
            except Exception:
                pass
            self.eng = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


def matlab_double(arr):
    """numpy (r,c) -> matlab.double (避免顶层 import matlab, 供 control 热路径用)."""
    import matlab
    return matlab.double(np.asarray(arr).tolist())
