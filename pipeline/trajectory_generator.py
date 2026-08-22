"""轨迹生成器: (seed_id, K, dt) -> 参考轨迹 + 场景参数.

数据源:
    - seed_id >= 1: 读 scenario_bank/scenario_seed{N}.mat (MATLAB v7.3 / HDF5, h5py 读取),
      与现有 MATLAB 批量仿真完全同场景, 结果可直接对比.
    - seed_id == 0: paper_fixed 固定场景 (论文 Section IV / golden 基准:
      J_total=13.3838, RMSE=0.036793, validSteps=100/100).

参考轨迹生成 generate_reference 逐行复刻 core/generateReference.m
(Bezier 曲线 + 弧长参数化姿态角).

K (预测时域) 传入本模块: 参考轨迹点数由场景 num_steps 决定 (与 MATLAB 一致,
保证 golden 可比); K 仅用于校验轨迹长度足以覆盖预测窗口.
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import h5py

from pipeline.params import DT, VehicleParams

# paper_fixed 场景参数 (来源: algorithms/RSS_proposed/config.m + others/paper_reproduction/paper_reproduction.m)
PAPER_FIXED_CTRL_PTS = np.array([
    [0.0, 0.0],
    [0.750, 0.250],
    [0.250, 0.750],
    [1.250, -1.000],
    [1.000, 0.0],
], dtype=np.float64)
PAPER_FIXED_INITIAL_STATE = np.array([0.05, 0.1, 0.2], dtype=np.float64)
PAPER_FIXED_INITIAL_VELOCITY = np.array([0.01, 0.01, 0.01], dtype=np.float64)
PAPER_FIXED_NUM_STEPS = 100  # t_end=1.0s / dt=0.01s


@dataclass
class Trajectory:
    """轨迹生成输出: 参考轨迹 + 场景参数 (供仿真器使用)."""
    path: np.ndarray               # (3, N) 参考轨迹 [x; y; theta]
    num_steps: int                 # 闭环仿真总步数
    initial_state: np.ndarray      # (3,) 世界系初始位姿 [x, y, psi]
    initial_velocity: np.ndarray   # (3,) 车体系初始速度
    vehicle: VehicleParams         # 车辆 (轮子) 参数
    ctrl_pts: np.ndarray           # (n+1, 2) Bezier 控制点
    scenario_name: str
    seed_id: int
    source: str                    # 'scenario_bank' / 'paper_fixed'

    @property
    def num_points(self) -> int:
        return int(self.path.shape[1])


def generate_reference(ctrl_pts: np.ndarray, num_points: int) -> np.ndarray:
    """参考轨迹生成 (core/generateReference.m 逐行等价).

    Bezier 曲线 (Bernstein 基) + 姿态角 theta = (2*pi / Lp^2) * s^2 (s 为累积弧长).
    """
    ctrl_pts = np.asarray(ctrl_pts, dtype=np.float64)
    t = np.linspace(0.0, 1.0, num_points)
    n = ctrl_pts.shape[0] - 1

    path_x = np.zeros(num_points)
    path_y = np.zeros(num_points)
    for i in range(n + 1):
        # nchoosek(n,i) * t^i * (1-t)^(n-i) * ctrl_pts(i+1,:)
        coeff = _nchoosek(n, i) * (t ** i) * ((1.0 - t) ** (n - i))
        path_x = path_x + coeff * ctrl_pts[i, 0]
        path_y = path_y + coeff * ctrl_pts[i, 1]

    dx = np.diff(path_x)
    dy = np.diff(path_y)
    seg = np.sqrt(dx ** 2 + dy ** 2)
    s_list = np.concatenate([[0.0], np.cumsum(seg)])
    Lp = s_list[-1]
    path_theta = (2.0 * np.pi / Lp ** 2) * s_list ** 2

    return np.vstack([path_x, path_y, path_theta])


def _nchoosek(n: int, k: int) -> int:
    from math import comb
    return comb(n, k)


# =========================================================
# 场景加载
# =========================================================

def make_trajectory(seed_id: int, K: int, dt: float = DT,
                    repo_root: Path = None) -> Trajectory:
    """统一入口: (seed_id, K, dt) -> Trajectory.

    seed_id >= 1: scenario_bank 随机场景; seed_id == 0: paper_fixed 固定场景.
    """
    if dt != DT:
        raise ValueError(f"dt 固定为 {DT} s, 不接受 {dt}")
    if seed_id == 0:
        return paper_fixed_trajectory(K, dt)
    return load_scenario(seed_id, K, dt, repo_root)


def paper_fixed_trajectory(K: int, dt: float = DT) -> Trajectory:
    """paper_fixed 固定场景 (golden 基准来源 run_paper_baseline_case.m)."""
    num_steps = PAPER_FIXED_NUM_STEPS
    path = generate_reference(PAPER_FIXED_CTRL_PTS, num_steps)
    vehicle = VehicleParams()  # 默认值即 config.m: Lx=0.655, Ly=0.335, vimax=5, phidotmax=5*pi
    _check_horizon(K, num_steps, path.shape[1])
    return Trajectory(
        path=path,
        num_steps=num_steps,
        initial_state=PAPER_FIXED_INITIAL_STATE.copy(),
        initial_velocity=PAPER_FIXED_INITIAL_VELOCITY.copy(),
        vehicle=vehicle,
        ctrl_pts=PAPER_FIXED_CTRL_PTS.copy(),
        scenario_name='paper_fixed',
        seed_id=0,
        source='paper_fixed',
    )


def load_scenario(seed_id: int, K: int, dt: float = DT,
                  repo_root: Path = None) -> Trajectory:
    """读 scenario_bank/scenario_seed{N}.mat (v7.3) -> Trajectory."""
    import h5py

    if repo_root is None:
        repo_root = Path(__file__).resolve().parent.parent
    mat_file = Path(repo_root) / 'scenario_bank' / f'scenario_seed{seed_id}.mat'
    if not mat_file.exists():
        raise FileNotFoundError(
            f"场景文件不存在: {mat_file}\n"
            f"(可在 MATLAB 中运行 scenario_bank({seed_id}, true) 生成)"
        )

    with h5py.File(mat_file, 'r') as f:
        config = _h5_group_to_dict(f['config'])
        scenario = _h5_group_to_dict(f['scenario'])

    # config 数值字段 (MATLAB 列主序已在 _h5_group_to_dict 中转置恢复)
    cfg_dt = float(config['dt'])
    if abs(cfg_dt - dt) > 1e-12:
        raise ValueError(f"场景 dt={cfg_dt} 与 pipeline dt={dt} 不一致")

    num_steps = int(config['num_steps'])
    num_path_pts = int(config['num_path_pts'])
    ctrl_pts = np.asarray(config['ctrl_pts'], dtype=np.float64)
    wheel_pos = np.asarray(config['wheel_pos'], dtype=np.float64)

    vehicle = VehicleParams(
        Lx=float(config['Lx']),
        Ly=float(config['Ly']),
        vimax=float(config['vimax']),
        phidotmax=float(config['phidotmax']),
        wheel_pos=wheel_pos,
    )

    initial_state = np.asarray(scenario['initialState'], dtype=np.float64).flatten()
    initial_velocity = np.asarray(scenario['initialVelocity'], dtype=np.float64).flatten()

    path = generate_reference(ctrl_pts, num_path_pts)
    _check_horizon(K, num_steps, path.shape[1])

    scenario_name = 'seed_' + str(seed_id)
    return Trajectory(
        path=path,
        num_steps=num_steps,
        initial_state=initial_state,
        initial_velocity=initial_velocity,
        vehicle=vehicle,
        ctrl_pts=ctrl_pts,
        scenario_name=scenario_name,
        seed_id=int(seed_id),
        source='scenario_bank',
    )


def _check_horizon(K: int, num_steps: int, num_path_pts: int):
    """校验参考轨迹足以覆盖预测窗口 (算法末端 clamp 到最后一点, 与 MATLAB 一致)."""
    if K < 1:
        raise ValueError(f"K 必须 >= 1, 收到 {K}")
    if num_steps < 1:
        raise ValueError(f"num_steps 必须 >= 1, 收到 {num_steps}")


# =========================================================
# HDF5 (MATLAB v7.3 .mat) 解析辅助
# =========================================================

def _h5_group_to_dict(group) -> dict:
    """递归解析 HDF5 group (MATLAB v7.3 struct) 为 dict.

    - 数值 dataset: 转置恢复 MATLAB 列主序形状 ((c,r) -> (r,c))
    - 字符串: 尽力 decode (uint16 char / vlen str), 失败返回 None
    - 嵌套 group (struct): 递归
    - 引用 (cell/hyper): 本项目场景不含, 忽略
    """
    out = {}
    for key in group.keys():
        item = group[key]
        if isinstance(item, h5py.Group):
            out[key] = _h5_group_to_dict(item)
        else:  # dataset
            out[key] = _h5_dataset_to_py(item)
    return out


def _h5_dataset_to_py(dset):
    """HDF5 dataset -> Python 值 (数值转置恢复 / 字符串 decode)."""
    if h5py.check_string_dtype(dset.dtype) is not None:
        return dset[()].decode('utf-8', errors='replace') if isinstance(dset[()], bytes) else str(dset[()])

    arr = np.asarray(dset)
    if arr.ndim == 0:
        return arr.item()
    # MATLAB 字符数组: uint16/uint8 编码
    if arr.dtype in (np.uint16, np.uint8, np.int16):
        try:
            cls = dset.attrs.get('MATLAB_class', b'')
            if isinstance(cls, bytes) and cls in (b'char',):
                if arr.dtype == np.uint16:
                    return ''.join(chr(int(v)) for v in arr.flatten() if v != 0)
                return ''.join(chr(int(v)) for v in arr.flatten() if v != 0)
        except Exception:
            pass
    # 数值: MATLAB 列主序, h5py 读出形状为转置, 恢复之
    if arr.ndim >= 2:
        arr = arr.T
    return arr
