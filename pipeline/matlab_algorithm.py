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

import numpy as np

MATLAB_BRIDGE_ALGORITHMS = ('e-lmpc', 'interior-point', 'active-set')


class MatlabAlgorithmBridge:
    """MATLAB 算法的每步求解桥 (供 simulator 循环调用)."""

    ALGO_SUBMODULE = {
        'e-lmpc': 'RSS_sqp',
        'interior-point': 'RSS_fmincon',
        'active-set': 'RSS_active_set',
    }

    def __init__(self, algorithm: str, trajectory, verbose: bool = True,
                 repo_root=None):
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
