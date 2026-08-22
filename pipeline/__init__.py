"""RSS-MPC 纯 Python 仿真 pipeline.

架构 (main.py 入口):
    main.py --seed N --algorithm proposed-3iter --K 6     [dt (tau) 写死 0.01 s]
      -> trajectory_generator.py   (seed, K -> 参考轨迹 + 场景参数)
      -> simulator.py              (原始动力学闭环, 循环调用算法)
      -> 结果三件套: run_config.json (参数记录)
                     + metrics.json (评估, 含求解时长)
                     + figures/ (画图)
                     + simulation_data.npz (原始数组)

算法 (proposed-3iter) 的实现归位算法包 algorithms/RSS_proposed/ (原 MATLAB 版已删除):
    construct_ocp_qcqp.py   (OCP QCQP 矩阵构造)
    control_rss_ocpqcqp.py  (SCP 外层控制器)
    hpipm_qp_solver.py      (HPIPM ctypes 接口)
pipeline/ 只留轨迹/仿真/评估主干 (simulator.py 加载上述控制器).
"""

from pipeline.params import DT, VehicleParams, CostWeights, AlgorithmParams, RunConfig
from pipeline.dynamics import (
    rotation_matrix,
    propagate_state,
    compute_wheel_outputs,
    wheel_matrices,
)

__all__ = [
    'DT', 'VehicleParams', 'CostWeights', 'AlgorithmParams', 'RunConfig',
    'rotation_matrix', 'propagate_state',
    'compute_wheel_outputs', 'wheel_matrices',
]
