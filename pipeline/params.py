"""运行参数定义。

dt (论文符号 tau) 按要求写死为 0.01 s, 不作为命令行可配置项。
其余默认值与 algorithms/RSS_proposed/config.m (golden 基准来源) 一致。
"""

from dataclasses import dataclass, field, asdict
import numpy as np

# tau: 离散化步长 (论文 IV-A: 0.01 s), 固定常量
DT = 0.01


@dataclass
class VehicleParams:
    """车辆几何与约束参数 (轮子参数), 对应 config.m 的 Lx/Ly/wheel_pos/vimax/phidotmax."""
    Lx: float = 0.655
    Ly: float = 0.335
    vimax: float = 5.0          # z_max: 最大轮速 (m/s), 论文 (20b)
    phidotmax: float = 5.0 * np.pi  # omega_max: 最大转向角速率 (rad/s), 论文 (5)
    wheel_pos: np.ndarray = field(default=None, repr=False)  # (N,2) 各轮位置 [dx, dy]

    def __post_init__(self):
        if self.wheel_pos is None:
            a = self.Lx / 2.0
            b = self.Ly / 2.0
            # 与 config.m / defaultConfig.m 相同的轮位顺序
            self.wheel_pos = np.array([
                [a, b],
                [-a, b],
                [-a, -b],
                [a, -b],
            ], dtype=np.float64)

    @property
    def num_wheels(self) -> int:
        return int(self.wheel_pos.shape[0])

    def to_dict(self) -> dict:
        return {
            'Lx': float(self.Lx),
            'Ly': float(self.Ly),
            'vimax': float(self.vimax),
            'phidotmax': float(self.phidotmax),
            'wheel_pos': self.wheel_pos.tolist(),
        }


@dataclass
class CostWeights:
    """代价函数权重 (论文公式 18) 与 RSS 正则化 (论文 17)."""
    w_pos: float = 30.0     # Q 位置分量
    w_psi: float = 1.0      # Q 姿态分量
    w_control: float = 0.3  # R 对角元
    rho: float = 0.01       # RSS 强凸正则化

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class SolverSettings:
    """HPIPM OCP QCQP 求解器设置 (记录用; 实际值写死于 hpipm_qp_solver.solve_ocp_qcqp).

    speed 模式 + tol 1e-8 对齐 ECOS 默认 (用户确认); warm_start=2 绕 DLL C 层
    d_ocp_qcqp_ipm_arg_set_t0_init 误写 t_lam_min 的 bug.

    integrator: 真实动力学 (plant) 状态推进方式, 4 算法统一生效, 不影响
    算法内部预测模型 (HPIPM A/B 保持论文欧拉离散):
        'euler' 定步长显式欧拉 (默认, golden 基准口径)
        'ode45' scipy RK45 (Dormand-Prince 4(5)) 单步积分, 纯 Python
                无 MATLAB 依赖 (dynamics.propagate_state_ode45)
    """
    mode: str = 'speed'
    tol: float = 1e-8
    iter_max: int = 1000
    mu0: float = 10.0
    warm_start: int = 2
    integrator: str = 'euler'

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class AlgorithmParams:
    """算法 (proposed-3iter) 每步求解所需的全部参数.

    由 main.py 组装, 经仿真器传给 algorithms/RSS_proposed/ 下的
    control_rss_ocpqcqp / construct_ocp_qcqp, 替代 MATLAB 版内部 config() 的角色.
    """
    K: int = 6                                  # 预测时域 (论文 IV-A: K=6; 现可参数化)
    dt: float = DT                              # 离散化步长 (固定 0.01)
    vehicle: VehicleParams = field(default_factory=VehicleParams)
    weights: CostWeights = field(default_factory=CostWeights)
    solver: SolverSettings = field(default_factory=SolverSettings)
    max_iter: int = 3                           # SCP 外层迭代数 (论文 IV-B: 固定 3 次)
    # 逐步 rho 序列 (CLI --rho 显式传入时非空):
    #   None        -> 用 weights.rho (各算法默认, 不改变基线行为)
    #   [v]         -> 常数 rho = v
    #   [v1, v2,..] -> 第 k 步取 v_k, 序列短于总步数时保持末值
    rho_schedule: list = None

    def to_dict(self) -> dict:
        return {
            'K': int(self.K),
            'dt': float(self.dt),
            'max_iter': int(self.max_iter),
            'vehicle': self.vehicle.to_dict(),
            'weights': self.weights.to_dict(),
            'solver': self.solver.to_dict(),
            'rho_schedule': (list(self.rho_schedule)
                             if self.rho_schedule else None),
        }


@dataclass
class RunConfig:
    """单次 pipeline 运行的完整配置 (写入 run_config.json)."""
    seed_id: int                                # 场景 seed (seedid)
    algorithm: str                              # 算法名 (目前仅 proposed-3iter)
    K: int                                      # 预测时域
    dt: float = DT                              # 写死 0.01
    num_steps: int = 0                          # 闭环仿真总步数 (场景决定)
    trajectory_source: str = ''                 # 场景来源 ('scenario_bank' / 'paper_fixed')
    scenario_name: str = ''
    algorithm_params: AlgorithmParams = field(default_factory=AlgorithmParams)

    def to_dict(self) -> dict:
        return {
            'seed_id': int(self.seed_id),
            'algorithm': self.algorithm,
            'K': int(self.K),
            'dt': float(self.dt),
            'num_steps': int(self.num_steps),
            'trajectory_source': self.trajectory_source,
            'scenario_name': self.scenario_name,
            'algorithm_params': self.algorithm_params.to_dict(),
        }
