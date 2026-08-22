"""原始动力学 (仿真器使用, 非算法内部的展开/线性化误差动力学).

逐行对应 MATLAB 版:
    propagate_state       <- core/propagateState.m
    compute_wheel_outputs <- core/computeWheelOutputs.m (轮子正运动学)

状态: xi = [x; y; psi] (世界系位姿)
车体速度: nu = [vx; vy; omega] (车体系)
世界速度: xdot = R(psi) * nu
"""

import numpy as np


def rotation_matrix(psi: float) -> np.ndarray:
    """论文 (1): R(psi) = [cos -sin 0; sin cos 0; 0 0 1]."""
    c, s = np.cos(psi), np.sin(psi)
    return np.array([
        [c, -s, 0.0],
        [s, c, 0.0],
        [0.0, 0.0, 1.0],
    ], dtype=np.float64)


def propagate_state(state: np.ndarray, world_velocity: np.ndarray, dt: float) -> np.ndarray:
    """原始动力学积分 (core/propagateState.m 等价).

    x_{k+1} = x_k + v_world_x * dt
    y_{k+1} = y_k + v_world_y * dt
    psi_{k+1} = psi_k + omega_world * dt
    """
    state = np.asarray(state, dtype=np.float64).flatten()
    world_velocity = np.asarray(world_velocity, dtype=np.float64).flatten()
    next_state = state.copy()
    next_state[0:2] = state[0:2] + world_velocity[0:2] * dt
    next_state[2] = state[2] + world_velocity[2] * dt
    return next_state


def wheel_matrices(wheel_pos: np.ndarray) -> list:
    """论文 (3): H_n = [1, 0, -dy_n; 0, 1, dx_n], 每轮一个 (2,3) 矩阵."""
    wheel_pos = np.asarray(wheel_pos, dtype=np.float64)
    return [
        np.array([
            [1.0, 0.0, -wheel_pos[n, 1]],
            [0.0, 1.0, wheel_pos[n, 0]],
        ], dtype=np.float64)
        for n in range(wheel_pos.shape[0])
    ]


def compute_wheel_outputs(body_velocity: np.ndarray, wheel_pos: np.ndarray):
    """轮子正运动学 (core/computeWheelOutputs.m 等价).

    输入: body_velocity (3,) 车体系速度; wheel_pos (N,2)
    输出: wheel_speed (N,) 各轮线速度; wheel_angle (N,) 各轮转向角 phi = atan2(vy, vx)
    """
    body_velocity = np.asarray(body_velocity, dtype=np.float64).flatten()
    num_wheels = np.asarray(wheel_pos).shape[0]
    wheel_speed = np.zeros(num_wheels)
    wheel_angle = np.zeros(num_wheels)
    for i, Hj in enumerate(wheel_matrices(wheel_pos)):
        wv = Hj @ body_velocity
        wheel_speed[i] = np.sqrt(wv[0] ** 2 + wv[1] ** 2)
        wheel_angle[i] = np.arctan2(wv[1], wv[0])
    return wheel_speed, wheel_angle
