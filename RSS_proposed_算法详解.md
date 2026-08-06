# RSS_proposed 算法文件详解

> 本文严格依据 `algorithms/RSS_proposed/` 下的实际代码编写，所有物理量、参数、矩阵构造均与代码一一对应。

## 论文公式索引

| 公式 | 出处 | 在代码中的对应 |
|---|---|---|
| (1) | RSS26 — 车体系->世界系变换 | `control_RSS.m` 输出 `new_state_dot = R(psi)*(state_dot + k1*u(:,1))`，其中 `state_dot+u(:,1)` 是车体系速度，`new_state_dot` 是世界系状态导数 |
| (3) | RSS26 — 特征矩阵 H_n | `construct_complete_qp_from_rss.m` `Hn{n}` |
| (4) | RSS26 — 对齐条件 z_n=H_n*nu | 用于轮速 SOC 约束和转向锥约束 |
| (9) | RSS26 — 离散化动力学 nu(:,k+1)=nu(:,k)+u(:,k+1) | 等式约束 A*x=b (18 条) |
| (11) | RSS26 — delta_theta=phidotmax*dt | `delta_theta = dt*phidotmax = pi/20` |
| (12) | RSS26 — 原始双线性转向锥约束 (非凸) | 经 Prop.1 凸化后实现 |
| (13) | RSS26 — 代入 z_n=H_n*nu(:,k) 后的双线性形式 | 转向锥凸化约束 |
| (15) | RSS26 — 凸化约束 C=A-B-L <= 0 | 转向锥二次约束 (48 条) |
| (16) | RSS26 — A/B/L 定义 | 转向锥 Hq/gq/uq 填充 |
| (17) | RSS26 — 凸子问题 Q_K(u_hat) 目标 f=g+rho*||u-u_hat||^2 | 代价函数 H/g 构造 |
| (18) | RSS26 — 代价 g(u)=sum e_k'*Q*e_k + u_k'*R*u_k | 位置+姿态+控制正则 |
| (19) | RSS26 — 跟踪误差 e_k 的 SO(2) 一阶展开 | 位置/姿态误差展开 |
| (20a) | RSS26 — 转向锥约束 (凸化后) | 转向锥二次约束 (48 条) |
| (20b) | RSS26 — 轮速 SOC 约束 ||H_n*nu(:,k)||<=vimax | 轮速二次约束 (24 条) |
| (20c) | RSS26 — 等式约束集 (动力学) | 等式约束 A*x=b |
| Alg.1 | RSS26 — Trajectory optimizer 算法 | `control_RSS.m` SQP 外层循环 |
| (1) | HPIPM — 完整 dense QP (含 slack) | 本代码使用硬约束子集 (无 slack) |

---

## 1. 文件总览

| 文件 | 作用 |
|---|---|
| [control_RSS.m](control_RSS.m) | 算法主入口，SQP 外层循环 + 调用 Python 求解器 |
| [construct_complete_qp_from_rss.m](construct_complete_qp_from_rss.m) | QP/QCQP 矩阵显式构造（H, g, A, b, Hq, gq, uq） |
| [hpipm_qp_solver.py](hpipm_qp_solver.py) | HPIPM dense QCQP 求解器 Python 封装 |
| [config.m](config.m) | 算法参数（机器人几何、仿真参数、路径参数） |
| [build_hpipm_windows.sh](build_hpipm_windows.sh) | Windows MSYS2 编译 blasfeo + hpipm 共享库脚本 |
| [.gitignore](.gitignore) | 忽略 MATLAB 自动备份与验证产物 |

---

## 2. 物理量与参数（与代码严格一致）

### 2.1 几何与限值（来自 `config.m`）

| 代码变量 | 代码值 | 含义 |
|---|---|---|
| `params.Lx` | `0.655` | 车身长度 (m) |
| `params.Ly` | `0.335` | 车身宽度 (m) |
| `params.a` | `Lx/2 = 0.3275` | 车轮 x 向偏移 |
| `params.b` | `Ly/2 = 0.1675` | 车轮 y 向偏移 |
| `params.wheel_pos` | `[a,b; -a,b; -a,-b; a,-b]` | 4 个车轮位置 |
| `params.vimax` | `5` | 最大轮速 (m/s) |
| `params.phidotmax` | `5*pi` | 最大转向角速率 (rad/s) |
| `params.dt` | `0.01` | 离散化时间步长 (s) |
| `params.t_end` | `1` | 仿真总时长 (s) |
| `params.num_steps` | `100` | 仿真步数 |

### 2.2 算法参数（来自 `control_RSS.m` 硬编码，覆盖 `config.m`）

| 代码变量 | 代码值 | 含义 |
|---|---|---|
| `K` | `6` | 预测时域长度 |
| `rho` | `0.01` | 强凸正则化参数 rho |
| `k1` | `1` | 输出增益（**硬编码为 1**，覆盖 `config.m` 中的 `0.15`） |
| `epsilon` | `0` | 收敛阈值（未启用，固定跑满迭代） |
| `max_iter` | `3` | SQP 最大迭代次数 |

### 2.3 代价权重（来自 `construct_complete_qp_from_rss.m`）

| 代码变量 | 代码值 | 含义 |
|---|---|---|
| `w_pos` | `30` | 位置跟踪权重（Q 的 xy 分量） |
| `w_psi` | `1` | 姿态跟踪权重（Q 的 ψ 分量） |
| `w_control` | `0.3` | 控制正则化权重（R 的对角元） |
| `rho` | `0.01` | RSS 强凸正则化参数 |

### 2.4 派生量（代码中计算）

| 代码变量 | 计算式 | 值 |
|---|---|---|
| `delta_theta` | `dt * phidotmax` | `0.01 * 5*pi = pi/20 ~= 0.1571` rad |
| `num_wheels` | `size(wheel_pos, 1)` | `4` |
| `n_var` | `6*K` | `36`（u 18 维 + nu 18 维） |
| `nu_start` | `3*K + 1` | `19`（nu 在 x 中的起始索引） |
| `n_eq` | `3*K` | `18`（3 初始 + 3*(K-1) 递推） |
| 二次约束总数 | `4*K + 2*4*K` | `24 + 48 = 72` |

### 2.5 决策变量排列

```
x = [u(:); nu(:)]  (36维)
```

- `u(i,k)` 全局索引 = `(k-1)*3 + i`，i=1,2,3，k=1..K -> 占索引 1..18
- `nu(i,k)` 全局索引 = `nu_start + (k-1)*3 + (i-1)` = `18 + (k-1)*3 + i` -> 占索引 19..36

---

## 3. `control_RSS.m` — 算法主入口

### 3.1 接口

```matlab
function [u, new_state_dot, velocity, diagnostics] = control_RSS(path, step, state_dot, state)
```

- 输入：参考轨迹 `path`、当前步 `step`、当前车体系速度 `state_dot`（即 `current_nu`）、当前位姿 `state`
- 输出：控制增量序列 `u`、世界系状态导数 `new_state_dot`、车体系速度 `velocity`、诊断 `diagnostics`

### 3.2 SQP 外层循环（对应论文 Algorithm 1）

**论文 Algorithm 1 (Trajectory optimizer for SWMRs)**：
```
1: Initialize u^(0) in ri(D(P_K))    % static init: u_hat = zeros(3, K)
2: m <- 0
3: while m < max_iter
4:   Construct convex subproblem Q_K(u^(m))    % 公式 (17)
5:   Solve u^(m+1) = S(u^(m))                   % 凸求解器
6:   if ||u^(m+1) - u^(m)|| < epsilon: break
7:   m <- m + 1
8: end while
9:   u_hat <- u^(m+1)                            % 无条件更新
10: end
11: return velocity = state_dot + k1 * u(:,1)   % 输出 (k1=1 硬编码)
```

```matlab
max_iter = 3;
u_hat = zeros(3, K);  % static init: u^(0) = 0

for m = 1 : max_iter
    % 1. 构造凸子问题 Q_K(u^(m))
    qp = construct_complete_qp_from_rss(path, step, current_nu, state, u_hat, params);

    % 2. 数据格式转换: cell -> 3D numpy
    Hq_3d = zeros(n_var, n_var, n_qcqp);
    gq_2d = zeros(n_var, n_qcqp);
    for i = 1:n_qcqp
        Hq_3d(:, :, i) = qp.Hq{i};
        gq_2d(:, i)    = qp.gq{i};
    end

    % 3. 求解 u^(m+1) = S(u^(m))
    result = py.hpipm_qp_solver.solve_qcqp(...
        py.numpy.array(qp.H), py.numpy.array(qp.g), ...
        py.numpy.array(qp.A), py.numpy.array(qp.b), ...
        py.numpy.array(Hq_3d), py.numpy.array(gq_2d), ...
        py.numpy.array(qp.uq'), py.bool(false));

    x = double(result{'x'});           % 36维解
    u_sol = reshape(x(1:3*K), 3, K);   % 提取 u

    % 4. 无条件更新 u_hat (论文 Alg.1 无 break-on-failure)
    u = u_sol;
    u_hat = u;  % 下次迭代凸化使用
end
```

### 3.3 输出（对应论文 Alg.1 line 11: velocity = state_dot + k1*u(:,1)）

**论文公式 (1)** — 车体系 -> 世界系变换：
```
new_state_dot = [R(psi), 0; 0, 1] * (state_dot + k1 * u(:,1))
```
其中：
- `state_dot` = 输入的当前车体系速度（即 `current_nu`）
- `state_dot + k1*u(:,1)` = 输出车体系速度 `velocity`（论文 Alg.1 line 11）
- `R(psi)` = 2D 旋转矩阵（由 `state(3)` 构造）
- `new_state_dot` = 输出的世界系状态导数

**论文 Alg.1 line 11**：`velocity = state_dot + k1*u(:,1)`（代码中 `k1=1` 硬编码）

```matlab
new_state_dot = [cos(psi), -sin(psi), 0;
                 sin(psi),  cos(psi), 0;
                    0,        0,    1] * (state_dot + 1.00 * u(:, 1));
velocity = current_nu + u(:, 1);
```

注意 `k1=1`（硬编码）。`state_dot + 1*u(:,1)` 即车体系速度 `velocity`，再经旋转矩阵 `R(psi)` 旋到世界系得到 `new_state_dot`。

### 3.4 Python 环境管理

为避免 `libhpipm.dll` 内存泄漏（WinError 8），使用 `persistent` 变量保证 Python 模块只加载一次：

```matlab
persistent py_path_added py_reloaded;
if isempty(py_path_added)
    sys_mod = py.importlib.import_module('sys');
    py.getattr(sys_mod, 'path').append(script_path);
    py_path_added = true;
end
if isempty(py_reloaded)
    solver_mod = py.importlib.import_module('hpipm_qp_solver');
    py.importlib.reload(solver_mod);
    py_reloaded = true;
end
```

---

## 4. `construct_complete_qp_from_rss.m` — QP 矩阵构造详解

### 4.1 接口与总览

```matlab
function qp_problem = construct_complete_qp_from_rss(path, step, current_nu, state, u_hat, params)
```

- 输入：参考轨迹、当前步、当前车体速度 `current_nu`、当前位姿、上次迭代解 `u_hat`、参数
- 输出：`qp_problem` 结构体，含 `H, g, A, b, C(空), d(空), lb(空), ub(空), Hq(cell), gq(cell), uq, n_var, K, n_eq, n_qcqp, objective_constant`

#### 论文公式 → H 填充位置总览

下表把 RSS 论文公式逐项翻译为 H 矩阵中的填充位置与系数（代码细节见 4.2 节）：

| 论文公式 | 论文原始形式 | H 中的位置 | H 中的系数 | 代码行 |
|---|---|---|---|---|
| (18) 位置 `e_k'Q e_k` | `w_pos·‖R·Σ ν·dt‖^2` | nu_x, nu_y 跨阶段交叉 | `2·w_pos·dt^2` | [L166-167](algorithms/RSS_proposed/construct_complete_qp_from_rss.m#L166-L167) |
| (18) 姿态 `e_k'Q e_k` | `w_psi·(dt·Σ ν_ψ)^2` | nu_psi 跨阶段交叉 | `2·w_psi·dt^2` | [L207](algorithms/RSS_proposed/construct_complete_qp_from_rss.m#L207) |
| (18) 控制 `u'Ru` | `w_control·‖u‖^2` | u 块对角 | `2·w_control` | [L227](algorithms/RSS_proposed/construct_complete_qp_from_rss.m#L227) |
| (17) 强凸项 | `rho·‖u‖^2` | u 块对角 | `2·rho` | [L241](algorithms/RSS_proposed/construct_complete_qp_from_rss.m#L241) |

### 4.2 目标函数 H, g 构造

#### 4.2.1 论文原始代价（公式 17-19）

**论文公式 (17)** — 凸子问题 Q_K(u_hat) 的目标函数：
```
f(u, u_hat) = g(u) + rho * ||u - u_hat||^2
```
其中 `g(u)` 为跟踪代价，`rho*||u-u_hat||^2` 为 RSS 强凸正则化（保证 S(u_hat) 唯一解，论文 Theorem 1）。

**论文公式 (18)** — 代价函数 g(u)：
```
g(u) = sum_{k=1}^K [ e_k' * Q * e_k + u_k' * R * u_k ]
```
其中 Q = diag(30, 30, 1)（位置 30，姿态 1），R = diag(0.3, 0.3, 0.3)（控制正则），e_k 为跟踪误差。

**论文公式 (19)** — 跟踪误差 e_k 的 SO(2) 一阶展开：
```
e_k = state + [R(psi0), 0; 0, 1] * sum_{l=0}^{k-1} nu(:,l) * dt - path(:, step+k) + O(dt^2)
```
其中 O(dt^2) 为高阶余项，dt=0.01s 时为 1e-4，代码中略去。

由 (17)+(18) 可知，目标函数由 **4 部分代价累加**：
1. 位置跟踪 `e_k'Q e_k`（位置部分）
2. 姿态跟踪 `e_k'Q e_k`（姿态部分）
3. 控制正则化 `u'Ru`
4. RSS 强凸正则化 `rho*||u-u_hat||^2`

#### 4.2.2 标准形式转换（论文 → HPIPM）

HPIPM dense QCQP 求解器要求的目标函数标准形式（[hpipm_qp_solver.py:28](algorithms/RSS_proposed/hpipm_qp_solver.py#L28)）：
```
min  0.5 x' H x + g' x
```

但论文公式 (18) 写成 `||·||^2 = v'v` 形式（无 0.5 系数）。要把论文代价套进 HPIPM 模板，必须做系数换算：

| 论文原始项 | 标准形式贡献 | 换算依据 |
|---|---|---|
| `w · ‖S‖^2` (二次) | `H += 2·w` | `0.5·x'(2w·I)x = w·‖x‖^2` |
| `2·w·c'·S` (一次) | `g += 2·w·c` | `g'x = 2·w·c'·x` |
| `w·‖c‖^2` (常数) | 累加到 `objective_constant` | 不含决策变量 |

**关键点**：代码里 `H_mat(...) += 2*w_pos*dt^2`（[L166-167](algorithms/RSS_proposed/construct_complete_qp_from_rss.m#L166-L167)）、`g_vec(...) += 2*w_pos*dt*grad_dir`（[L176-177](algorithms/RSS_proposed/construct_complete_qp_from_rss.m#L176-L177)）等系数都带 `2`，就是因为这个换算——把 `‖v‖^2` 形式适配到 HPIPM 的 `0.5·x'Hx` 形式。

##### 表中 `S` 和 `c` 的具体含义

表中 `S` 和 `c` **不是论文原文符号**，而是把论文 `e_k'Q e_k` 展开后的**中间变量**，对应代码中的实际变量：

| 简化记号 | 位置部分代码 | 姿态部分代码 |
|---|---|---|
| `S` | `R_psi0·dt·S_k`（其中 `S_k = Σ_{j=1}^{k-1} nu(1:2,j)`） | `dt·sum(nu(3,1:k-1))` |
| `c` | `c_k = current_xy - ref_xy + R_psi0·current_nu(1:2)·dt` | `psi_c_k = psi0 - ref_psi + current_nu(3)·dt` |
| `w` | `w_pos = 30` | `w_psi = 1` |

**对齐说明**：两边的 `c` 形式统一为 `当前状态 - 参考 + 第0步位移`：

```
位置: c_k      = current_xy - ref_xy    + R_psi0·current_nu(1:2)·dt
                └─当前位置(2维)─┘   └─参考(2维)─┘   └──第0步位移(2维)──┘
姿态: psi_c_k  = psi0       - ref_psi   + current_nu(3)·dt
                └─当前姿态(1维)─┘  └─参考(1维)─┘  └─第0步位移(1维)┘
```

- 第1项：当前状态（位置 `current_xy` / 姿态 `psi0`）
- 第2项：参考轨迹（`ref_xy` / `ref_psi`）
- 第3项：第 0 步位移（l=0 的速度积分，用已知量 `current_nu`）
  - 位置部分需 `R_psi0` 旋转（车体系→世界系）
  - 姿态部分是标量，无需旋转

**展开过程（以位置为例）**：

跟踪误差展开（基于论文公式 19）：
```
e_k_xy = c_k + R_psi0·dt·S_k
         │      └────── S ──────┘
         └── c ──┘
```

代价 `w_pos·‖e_k_xy‖^2` 展开为三部分：
```
w_pos·‖c_k + R_psi0·dt·S_k‖^2
= w_pos·‖c_k‖^2                          ← w·‖c‖^2  (常数 → objective_constant)
+ 2·w_pos·(c_k)'·R_psi0·dt·S_k           ← 2·w·c'·S (一次 → g_vec)
+ w_pos·‖R_psi0·dt·S_k‖^2                ← w·‖S‖^2  (二次 → H_mat)
```

利用 `R_psi0'·R_psi0 = I`（旋转矩阵正交性），二次项中 R 自动消去：
```
‖R_psi0·dt·S_k‖^2 = dt^2·‖S_k‖^2 = dt^2·Σ_{i,j} nu(:,i)'·nu(:,j)
```

这就是 H 中 **nu 块跨阶段交叉项**的来源。

**代码里的实际变量名**：

位置部分（[L153-179](algorithms/RSS_proposed/construct_complete_qp_from_rss.m#L153-L179)）：
```matlab
c_k = current_xy - ref_xy + R_psi0 * current_nu(1:2) * dt;  % ← 这就是 c
grad_dir = R_psi0' * c_k;                                    % c 经 R' 转换
% ...
H_mat(...) += 2*w_pos*dt^2;        % 二次项 (S)
g_vec(...) += 2*w_pos*dt*grad_dir; % 一次项 (c)
% ...
objective_constant += w_pos*(c_k'*c_k);  % 常数项 (c) [L262]
```

姿态部分（[L199-216](algorithms/RSS_proposed/construct_complete_qp_from_rss.m#L199-L216)）：
```matlab
psi_c_k = psi0 + current_nu(3)*dt - ref_psi;  % ← 这就是 c
% ...
H_mat(...) += 2*w_psi*dt^2;        % 二次项 (S)
g_vec(...) += 2*w_psi*dt*psi_c_k;  % 一次项 (c)
% ...
objective_constant += w_psi*psi_c_k^2;  % 常数项 (c) [L270]
```

**总结**：
- `c` = `c_k`（位置）/ `psi_c_k`（姿态）= 跟踪误差的**常数部分**（不含决策变量）
- `S` = `R_psi0·dt·S_k`（位置）/ `dt·Σ nu(3,j)`（姿态）= 跟踪误差的**决策变量部分**
- `w` = `w_pos`（位置）/ `w_psi`（姿态）= Q 矩阵的对角权重

表里用 `S` 和 `c` 是为了抽象出共同的二次型展开模式 `‖c+S‖^2 = ‖c‖^2 + 2c'S + ‖S‖^2`，让位置和姿态两部分的换算逻辑统一。

#### 4.2.3 H 矩阵最终结构

```
                    u 块 (1..18)                nu 块 (19..36)
              ┌────────────────────────┬──────────────────────────┐
              │                        │                          │
   u 块       │ 2*(w_control + rho)*I  │      0                   │
   (1..18)    │  (纯对角)              │   (u 与 nu 无交叉)        │
              │                        │                          │
              ├────────────────────────┼──────────────────────────┤
              │                        │  位置 Hessian:           │
              │                        │   nu_x, nu_y 跨阶段交叉  │
   nu 块      │      0                 │  姿态 Hessian:           │
   (19..36)   │                        │   nu_psi 跨阶段交叉      │
              │                        │  (nu_x/y/psi 之间无交叉) │
              └────────────────────────┴──────────────────────────┘
```

- **u 块**：纯对角（控制正则 + RSS 正则），每个 u(i,k) 独立
- **nu 块**：有跨阶段交叉项（因为 `Σ ν_l·dt` 的平方展开后产生 (i,j) 对的交叉）
- **u 与 nu 无交叉项**

#### 4.2.4 逐项填充代码

##### (1) 位置跟踪代价（k=2..K，对应论文 (18) 位置 + (19) 展开）

跟踪误差的一阶展开（略去 O(dt^2)）：
```
position_error_k = current_xy - ref_xy_k + R(psi0)*(current_nu(1:2)*dt + sum_{j=1}^{k-1} nu(1:2,j)*dt)
                 = c_k + R_psi0*dt*S_k
```
其中 `c_k = current_xy - ref_xy + R_psi0*current_nu(1:2)*dt`（常数），`S_k = sum_{j=1}^{k-1} nu(1:2,j)`（决策变量）。

代价 `w_pos*||c_k + R_psi0*dt*S_k||^2` 展开（利用 `R_psi0'*R_psi0 = I`）：
```
= w_pos*||c_k||^2                          ← 常数 → objective_constant
+ 2*w_pos*dt*(R_psi0'*c_k)'*S_k            ← 一次项 → g_vec
+ w_pos*dt^2*||S_k||^2                     ← 二次项 → H_mat
```

```matlab
for k = 2:K
    ref_xy = path(1:2, min(size(path,2), step+k));
    c_k = current_xy - ref_xy + R_psi0 * current_nu(1:2) * dt;
    grad_dir = R_psi0' * c_k;  % 2×1

    % 二次项: H(nu_x(i), nu_x(j)) += 2*w_pos*dt^2 (对 i,j in [1, k-1])
    for i = 1:k-1
        for j = 1:k-1
            nu_x_i = nu_start + (i-1)*3;
            nu_x_j = nu_start + (j-1)*3;
            nu_y_i = nu_start + (i-1)*3 + 1;
            nu_y_j = nu_start + (j-1)*3 + 1;
            H_mat(nu_x_i, nu_x_j) = H_mat(nu_x_i, nu_x_j) + 2 * w_pos * dt^2;
            H_mat(nu_y_i, nu_y_j) = H_mat(nu_y_i, nu_y_j) + 2 * w_pos * dt^2;
        end
    end

    % 一次项: g(nu_x(j)) += 2*w_pos*dt*grad_dir(1) (对 j in [1, k-1])
    for j = 1:k-1
        nu_x_j = nu_start + (j-1)*3;
        nu_y_j = nu_start + (j-1)*3 + 1;
        g_vec(nu_x_j) = g_vec(nu_x_j) + 2 * w_pos * dt * grad_dir(1);
        g_vec(nu_y_j) = g_vec(nu_y_j) + 2 * w_pos * dt * grad_dir(2);
    end
end
```

填充位置：**nu 块的 x 行、y 行的跨阶段交叉**（i,j ∈ [1,k-1]）。

##### (2) 姿态跟踪代价（k=1..K，对应论文 (18) 姿态 + (19) 展开）

姿态误差：`psi(k) - ref_psi_k`，其中
```
psi(k) = psi0 + current_nu(3)*dt + dt*sum_{j=1}^{k-1} nu(3,j)
```
常数部分：`psi_c_k = psi0 + current_nu(3)*dt - ref_psi_k`

代价 `w_psi*(psi_c_k + dt*Σ nu_psi)^2` 展开：
```
= w_psi*psi_c_k^2                          ← 常数
+ 2*w_psi*dt*psi_c_k*Σ nu_psi              ← 一次项 → g_vec
+ w_psi*dt^2*(Σ nu_psi)^2                  ← 二次项 → H_mat
```

```matlab
for k = 1:K
    ref_psi = path(3, min(size(path,2), step+k));
    psi_c_k = psi0 + current_nu(3)*dt - ref_psi;

    % 二次项: H(nu_psi(i), nu_psi(j)) += 2*w_psi*dt^2
    for i = 1:k-1
        for j = 1:k-1
            nu_psi_i = nu_start + (i-1)*3 + 2;
            nu_psi_j = nu_start + (j-1)*3 + 2;
            H_mat(nu_psi_i, nu_psi_j) = H_mat(nu_psi_i, nu_psi_j) + 2 * w_psi * dt^2;
        end
    end

    % 一次项: g(nu_psi(j)) += 2*w_psi*dt*psi_c_k
    for j = 1:k-1
        nu_psi_j = nu_start + (j-1)*3 + 2;
        g_vec(nu_psi_j) = g_vec(nu_psi_j) + 2 * w_psi * dt * psi_c_k;
    end
end
```

**k=1 时**：循环 `1:0` 为空，不贡献二次/一次项，只贡献常数项 `w_psi*psi_c_1^2`。

填充位置：**nu 块的 ψ 行（nu(3,:)）的跨阶段交叉**。

##### (3) 控制正则化（k=1..K，对应论文 (18) 第二项）

```
u_k' * R * u_k,   R = diag(0.3, 0.3, 0.3)
```

代价 `w_control*||u||^2 = w_control*sum u(i,k)^2`，由于每个分量独立，H 只在对角有值：

```matlab
for k = 1:K
    for i = 1:3
        u_idx = (k-1)*3 + i;
        H_mat(u_idx, u_idx) = H_mat(u_idx, u_idx) + 2 * w_control;
    end
end
```

填充位置：**u 块对角**（u(1:18) 的对角元）。

##### (4) RSS 强凸正则化（k=1..K，对应论文 (17) 第二项）

```
rho * ||u - u_hat||^2    (提供强凸性, 保证 S(u_hat) 唯一解, 论文 Theorem 1)
```

展开 `rho*||u - u_hat||^2 = rho*(||u||^2 - 2*u'*u_hat + ||u_hat||^2)`：
- `rho*||u||^2` → 二次项进 H（u 对角）
- `-2*rho*u'*u_hat` → 一次项进 g
- `rho*||u_hat||^2` → 常数

```matlab
for k = 1:K
    for i = 1:3
        u_idx = (k-1)*3 + i;
        H_mat(u_idx, u_idx) = H_mat(u_idx, u_idx) + 2 * rho;
        g_vec(u_idx) = g_vec(u_idx) - 2 * rho * u_hat(i, k);
    end
end
```

填充位置：**u 块对角**（与项 (3) 叠加）。

#### 4.2.5 H 对称化与常数项收尾

```matlab
H_mat = 0.5 * (H_mat + H_mat');  % 保证对称 (数值稳定性)
```

**常数项 `objective_constant`** 收集所有不含决策变量的常数项（[L255-274](algorithms/RSS_proposed/construct_complete_qp_from_rss.m#L255-L274)）：
```
objective_constant = sum_{k=2}^K w_pos*||c_k||^2
                   + sum_{k=1}^K w_psi*psi_c_k^2
                   + rho*||u_hat||^2
```

它**不影响最优解 x\***（HPIPM 求解时不需要），只用于事后计算真实的 `obj_value` 做指标对比（论文 Table II 的 J_total）。

#### 4.2.6 拆解总结

```
0.5 x' H x  ← 二次项: 位置跟踪、姿态跟踪、控制正则化、RSS 强凸正则
   g' x     ← 一次项: 跟踪误差的常数偏移、RSS 正则的 -2*rho*u_hat 项
   const    ← 常数项: ||c_k||^2、psi_c_k^2、rho*||u_hat||^2 (不影响最优解, 仅用于指标)
```

这种分离是为了把论文公式 (18) 的原始代价 `e'Qe + u'Ru + rho*||u-u_hat||^2` 套进 HPIPM 求解器规定的 `0.5*x'Hx + g'x` 模板。

---

### 4.3 预备量

#### 4.3.1 特征矩阵 H_n（2×3）

**论文公式 (3)** — 第 n 个轮子的特征矩阵：
```
H_n = [1, 0, -dy_n; 0, 1, dx_n]
```

**论文公式 (4)** — 对齐条件（z_n 为轮子平面投影的车体速度）：
```
z_n = H_n * nu
```

```matlab
for n = 1:num_wheels
    Hn{n} = [1, 0, -wheel_pos(n, 2);   % -dy_n
                 0, 1,  wheel_pos(n, 1)];   %  dx_n
end
```

M_n = H_n' * H_n（3×3），用于轮速二次约束。

#### 4.3.2 旋转矩阵 R(psi0)

```matlab
R_psi0 = [cos(psi0), -sin(psi0); sin(psi0), cos(psi0)];
```

#### 4.3.3 delta_theta

**论文公式 (11)** — 单步最大转向角变化：
```
delta_theta = phidotmax * dt
```

```matlab
delta_theta = dt * phidotmax;  % = 0.01 * 5*pi = pi/20
```

---

### 4.4 等式约束 A*x = b（动力学递推，共 18 条）

**论文公式 (9)** — 离散化动力学（车体系速度递推）：
```
nu(:,k+1) = nu(:,k) + u(:,k+1),   k = 0, 1, ..., K-1
```

**论文公式 (20c)** — 等式约束集（包含在 Q_K 中）：
```
u in {u : nu(:,k+1) = nu(:,k) + u(:,k+1)}
```

动力学：`nu(:,k+1) = nu(:,k) + u(:,k+1)`

#### 4.4.1 初始条件（3 条）

```matlab
% nu(:,1) - u(:,1) = current_nu  即 nu(:,1) = current_nu + u(:,1)
for i = 1:3
    u_idx = i;                       % u(i,1)
    nu_idx = nu_start + i - 1;       % nu(i,1)
    A_eq(eq_row, u_idx)  = -1;
    A_eq(eq_row, nu_idx) = 1;
    b_eq(eq_row) = current_nu(i);
    eq_row = eq_row + 1;
end
```

#### 4.4.2 递推约束（3*(K-1) = 15 条）

```matlab
% nu(:,k+1) - nu(:,k) - u(:,k+1) = 0
for k = 1:K-1
    for i = 1:3
        u_idx      = k*3 + i;                      % u(i,k+1)
        nu_k_idx   = nu_start + (k-1)*3 + i - 1;   % nu(i,k)
        nu_kp1_idx = nu_start + k*3 + i - 1;       % nu(i,k+1)
        A_eq(eq_row, u_idx)      = -1;
        A_eq(eq_row, nu_k_idx)   = -1;
        A_eq(eq_row, nu_kp1_idx) = 1;
        b_eq(eq_row) = 0;
        eq_row = eq_row + 1;
    end
end
```

---

### 4.5 二次不等式约束（共 72 条）

QCQP 形式：`0.5*x'*Hq_i*x + gq_i'*x <= uq_i`

**论文公式 (20a)** — 转向锥约束（凸化后）+ **论文公式 (20b)** — 轮速 SOC 约束：
```
(20a): C^k_{i,n}(u, u_hat) <= 0     (转向锥, 凸化后)
(20b): ||H_n * nu(:,k)|| <= vimax   (轮速, 原本就是凸的)
```

#### 4.5.1 轮速 SOC 约束（4×K = 24 条）

**论文公式 (20b)** — 轮速约束：
```
||H_n * nu(:,k)|| <= vimax,   for all n in {1,...,N}, k in {1,...,K}
```

`||H_n*nu(:,k)|| <= vimax` -> `nu(:,k)'*M_n*nu(:,k) <= vimax^2`

```matlab
for k = 1:K
    for n = 1:num_wheels
        Hq_k = zeros(n_var, n_var);
        gq_k = zeros(n_var, 1);
        nu_k_start = nu_start + (k-1)*3;
        nu_k_end   = nu_start + (k-1)*3 + 2;
        % Hq 在 nu(:,k) 块 = 2*(H_n'*H_n) (因 0.5*x'*Hq*x 形式)
        Hq_k(nu_k_start:nu_k_end, nu_k_start:nu_k_end) = 2 * (Hn{n}' * Hn{n});
        Hq_list{end+1} = Hq_k;
        gq_list{end+1} = gq_k;
        uq_list(end+1) = vimax^2;  % = 25
    end
end
```

#### 4.5.2 转向锥凸化约束（2×4×K = 48 条）

**论文公式 (12)** — 原始双线性转向锥约束（非凸）：
```
(R_i * z_n(t_k))' * z_n(t_{k+1}) >= 0,   i in {1, 2}
R_1 = R(pi/2 - delta_theta),   R_2 = R_1'
```

**论文公式 (13)** — 代入 z_n = H_n*nu(:,k)（论文 (4)）后：
```
nu(:,k-1)' * H_n' * R_i' * H_n * (nu(:,k-1) + u(:,k)) >= 0   (非凸双线性)
```

**论文公式 (15)** — Proposition 1 凸化后的约束：
```
C^k_{i,n}(u, u_hat) = A^k_{i,n}(u) - B^k_{i,n}(u_hat) - L^k_{i,n}(u, u_hat) <= 0
```

**论文公式 (16)** — A, B, L 定义：
```
A^k_{i,n}(u) = 1/2 * ||H_n*nu(:,k-1)||^2 + 1/2 * ||H_n*(nu(:,k-1)+u(:,k))||^2   (凸, 关于 u)
B^k_{i,n}(u_hat) = 1/2 * ||(I+R_i)*H_n*nu_hat(:,k-1) + R_i*H_n*u_hat(:,k)||^2     (常数, 在 u_hat 处)
L^k_{i,n}(u, u_hat) = ∇_u B|_u_hat * (u - u_hat)                                 (B 的一阶展开, 线性)
```

两组旋转矩阵：
```matlab
R1 = [sin(delta_theta), -cos(delta_theta); cos(delta_theta),  sin(delta_theta)];  % R(pi/2-delta_theta)
R2 = [sin(delta_theta),  cos(delta_theta); -cos(delta_theta), sin(delta_theta)]; % R1'
```

构造 `nu_hat`（u_hat 对应的 nu_hat 序列，用于凸化的 B 项）：
```matlab
nu_hat = zeros(3, K);
nu_hat(:, 1) = current_nu + u_hat(:, 1);
for k = 1:K-1
    nu_hat(:, k+1) = nu_hat(:, k) + u_hat(:, k+1);
end
```

子函数 `add_steering_cone_constraints` 对每个 (k, n) 构造一条约束。约束形式（**论文公式 (15)**）：
```
C^k_{i,n}(u, u_hat) = A(u) - B(u_hat) - L(u, u_hat) <= 0
```

其中 A/B/L 由 **论文公式 (16)** 定义：
```
A(u) = 1/2*||H_n*nu(:,k-1)||^2 + 1/2*||H_n*(nu(:,k-1)+u(:,k))||^2    (二次, 凸)
B(u_hat) = 1/2*||(I+R_i)*H_n*nu_hat(:,k-1) + R_i*H_n*u_hat(:,k)||^2     (常数, 在 u_hat 处)
L(u, u_hat) = ∇_u B|_u_hat * (u - u_hat)                                (线性, B 的一阶展开)
```

##### k > 1 情况（nu(:,k-1) 是决策变量）

```matlab
% b = (I+R)*H_n*nu_hat(:,k-1) + R*H_n*u_hat(:,k)  (B 项的核)
lv = (eye(2) + R) * Hn * nu_hat(:, k-1) + R * Hn * u_hat(:, k);

% A 项 (二次): ||H_n*nu(:,k-1)||^2 + ||H_n*(nu(:,k-1)+u(:,k))||^2
%   = nu(:,k-1)'*M_n*nu(:,k-1) + (nu(:,k-1)+u(:,k))'*M_n*(nu(:,k-1)+u(:,k))
% 在 0.5*x'*Hq*x 形式:
%   Hq(nu(:,k-1), nu(:,k-1)) = 4*M_n  (两项各贡献 2*M_n)
%   Hq(nu(:,k-1), u(:,k))    = 2*M_n
%   Hq(u(:,k), u(:,k))        = 2*M_n
nu_k1_start = nu_start + (k-2)*3;
nu_k1_end   = nu_start + (k-2)*3 + 2;
u_k_start   = (k-1)*3 + 1;
u_k_end     = k*3;
Hq_k(nu_k1_start:nu_k1_end, nu_k1_start:nu_k1_end) = 4 * Mn;
Hq_k(nu_k1_start:nu_k1_end, u_k_start:u_k_end)     = 2 * Mn;
Hq_k(u_k_start:u_k_end, nu_k1_start:nu_k1_end)     = 2 * Mn;
Hq_k(u_k_start:u_k_end, u_k_start:u_k_end)         = 2 * Mn;

% L 项 (线性): -∇B*(u-u_hat)
%   对 u(:,l) (l<k): gq += -2*((I+R)'*lv)'*H_n
%   对 u(:,k):       gq += -2*(R'*lv)'*H_n
for l = 1:k-1
    coeff = -2 * ((eye(2) + R)' * lv)' * Hn;
    u_l_start = (l-1)*3 + 1;
    gq_k(u_l_start:u_l_start+2) = gq_k(u_l_start:u_l_start+2) + coeff';
end
coeff_k = -2 * (R' * lv)' * Hn;
gq_k(u_k_start:u_k_start+2) = gq_k(u_k_start:u_k_start+2) + coeff_k';

% B + L 中的 u_hat 常数 (移到 uq)
u_hat_const = 0;
for l = 1:k-1
    u_hat_const = u_hat_const + 2 * lv' * (eye(2) + R) * Hn * u_hat(:, l);
end
u_hat_const = u_hat_const + 2 * lv' * R * Hn * u_hat(:, k);
uq_val = lv' * lv - u_hat_const;  % ||b||^2 - L 中的 u_hat 常数部分
```

##### k = 1 情况（nu(:,0) = current_nu 是已知常数）

```matlab
% b = (I+R)*H_n*current_nu + R*H_n*u_hat(:,1)
lv = (eye(2) + R) * Hn * current_nu + R * Hn * u_hat(:, 1);

% A 项: 只有 u(:,1)'*M_n*u(:,1) 是决策变量部分 (current_nu 常数移到 uq)
Hq_k(1:3, 1:3) = 2 * Mn;

% 线性项: 2*current_nu'*M_n*u(:,1) (来自 A) - 2*b'*R*H_n*u(:,1) (来自 -L)
gq_k(1:3) = 2 * Mn' * current_nu - 2 * Hn' * R' * lv;

% 常数项: -2*current_nu'*M_n*current_nu (A 常数) + ||b||^2 (B) - 2*b'*R*H_n*u_hat(:,1) (L 中 u_hat 常数)
uq_val = -2 * current_nu' * Mn * current_nu + lv' * lv - 2 * lv' * R * Hn * u_hat(:, 1);
```

最后 Hq_k 对称化：
```matlab
Hq_k = 0.5 * (Hq_k + Hq_k');
```

---

### 4.6 返回结构体

```matlab
qp_problem.H  = H_mat;       % 代价 Hessian
qp_problem.g  = g_vec;       % 代价线性项
qp_problem.A  = A_eq;        % 等式约束矩阵 (18×36)
qp_problem.b  = b_eq;        % 等式约束右端 (18×1)
qp_problem.C  = C_ineq;      % 一般线性不等式矩阵 (空, ng=0)
qp_problem.d  = d_ineq;      % 一般线性不等式右端 (空)
qp_problem.lb = [];          % box 下界 (空, nb=0)
qp_problem.ub = [];          % box 上界 (空, nb=0)
qp_problem.Hq = Hq_list;     % 二次约束 Hessian cell (72 个 36×36)
qp_problem.gq = gq_list;     % 二次约束线性项 cell (72 个 36×1)
qp_problem.uq = uq_list;     % 二次约束右端 (72×1)
qp_problem.n_var = 36;
qp_problem.K = 6;
qp_problem.n_eq = 18;
qp_problem.n_qcqp = 72;
qp_problem.objective_constant = objective_constant;
```

**关于 `C`/`d`/`lb`/`ub` 为空**：对应 HPIPM dense QP 公式 (1) 中的：
- `C`, `d` — 一般线性不等式 `d_lb <= C*v <= d_ub`（本代码未启用，设 `ng=0`）
- `lb`, `ub` — box 约束 `v_lb <= v <= v_ub`（本代码未启用，设 `nb=0`）

本代码使用 HPIPM 硬约束子集（`nb=0, ng=0, ns=0, 无 slack`），退化为纯 QCQP。

---

## 5. `hpipm_qp_solver.py` — Python 求解器封装

**HPIPM 论文**（arXiv:2003.02547）Section 2.1 公式 (1) — 完整 dense QP 形式：
```
min_{v,s}  1/2 [v;1]^T [H  g; g^T  0] [v;1]
          + 1/2 [s^l;s^u;1]^T [Z^l  0  z^l; 0  Z^u  z^u; (z^l)^T (z^u)^T  0] [s^l;s^u;1]
s.t. A v = b                                                              (等式)
     [v_lb; d_lb] <= [J^{b,v}; C] v + [J^{s,v}; J^{s,g}] s^l              (下界+slack)
     [J^{b,v}; C] v - [J^{s,v}; J^{s,g}] s^u <= [v_ub; d_ub]              (上界+slack)
     s^l >= s^l_lb,  s^u >= s^u_lb                                         (slack 非负)
```

**dense QCQP 扩展** — 在 dense QP 基础上增加二次约束：
```
0.5 v^T Hq_i v + gq_i^T v <= uq_i   (二次不等式, 亦可带 slack)
```

**本代码使用硬约束子集**（`nb=0, ng=0, ns=0, 无 slack`）：
```
min  0.5 x^T H x + g^T x
s.t. A x = b                                  (等式, 论文 (20c) 动力学)
     0.5 x^T Hq_i x + gq_i^T x <= uq_i         (二次不等式, 论文 (20a)+(20b))
```

### 5.1 接口

```python
def solve_qcqp(H, g, A, b, Hq, gq, uq, verbose=False) -> dict
```

- 输入：Hessian `H`、线性项 `g`、等式约束 `A, b`、二次约束 `Hq(3D), gq(2D), uq`
- 返回：`{'x', 'status', 'status_str', 'obj_value', 'solve_time', 'iters'}`

> **变量名说明**：HPIPM C/Python wrapper 内部决策变量名为 `v`（`qcqp_sol.get('v')`），与论文符号一致；返回 dict 的键名为 `'x'`，与 MATLAB 端 `result{'x'}` 对齐。

### 5.2 HPIPM 维度设置

对应 **HPIPM 论文 Section II-B** 的数据结构（维度对象）：

```python
dim = hpipm_dense_qcqp_dim()
dim.set('nv', n)    # 36  (决策变量数)
dim.set('ne', ne)   # 18  (等式约束数)
dim.set('nb', 0)    # 0   (box bounds, 未启用)
dim.set('ng', 0)    # 0   (一般线性约束, 未启用)
dim.set('nq', nq)   # 72  (二次约束数)
```

### 5.3 求解流程

对应 **HPIPM 论文 Section III** — 求解器模式：
- `balance`：平衡模式（默认，平衡速度与鲁棒性）
- `robust`：鲁棒模式（更稳定但更慢，balance 失败时回退）

```python
# balance 模式 (默认)
arg = hpipm_dense_qcqp_solver_arg(dim, 'balance')
arg.set('iter_max', 300)
arg.set('tol_stat', 1e-6)
arg.set('tol_eq', 1e-6)
arg.set('tol_ineq', 1e-6)
arg.set('tol_comp', 1e-6)

solver = hpipm_dense_qcqp_solver(dim, arg)
solver.solve(qcqp, qcqp_sol)
x = qcqp_sol.get('v').flatten()
status = int(solver.get('status'))  # 0=SUCCESS (HPIPM 论文 Table I)

# balance 失败时回退到 robust 模式
if status != 0:
    arg2 = hpipm_dense_qcqp_solver_arg(dim, 'robust')
    arg2.set('iter_max', 500)
    # ... 同样的 tol 设置
    solver2 = hpipm_dense_qcqp_solver(dim, arg2)
    solver2.solve(qcqp, qcqp_sol)
    x = qcqp_sol.get('v').flatten()
    status = int(solver2.get('status'))
```

### 5.4 DLL 路径自动设置

```python
_hpipm_lib_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', '..', 'third_party', 'hpipm', 'lib')
if sys.platform.startswith('win'):
    os.add_dll_directory(_hpipm_lib_dir)  # Python 3.8+ DLL 搜索路径
    # 同时加入 PATH 作为 fallback
```

---

## 6. `config.m` — 算法参数

```matlab
params.Lx = 0.655;          % 车身长度 (m)
params.Ly = 0.335;          % 车身宽度 (m)
params.a = params.Lx/2;     % 车轮 x 向偏移
params.b = params.Ly/2;     % 车轮 y 向偏移
params.wheel_pos = [a,b; -a,b; -a,-b; a,-b];  % 4 个车轮位置
params.vimax = 5;           % 最大轮速 (m/s)
params.phidotmax = 5*pi;   % 最大转向角速率 (rad/s)
params.k1 = 0.15;           % 输出增益 (被 control_RSS.m 硬编码为 1 覆盖)
params.k2 = 0.15;
params.k3 = 0.1;
params.eps = 0.001;
params.dt = 0.01;           % 时间步长 (s)
params.t_end = 1;           % 仿真时长 (s)
params.num_steps = 100;
params.ctrl_pts = [...];    % 贝塞尔曲线控制点
params.num_path_pts = 100;
```

> **注意**：`config.m` 中的 `k1=0.15` 在 `control_RSS.m` 中被硬编码的 `k1=1` 覆盖。

---

## 7. `build_hpipm_windows.sh` — 编译脚本

Windows 下用 MSYS2 编译 `libhpipm.dll`：

```bash
# 编译 blasfeo 静态库
cd $BLASFEO_DIR
make -j4 static_library TARGET=GENERIC USE_C99_MATH=1 EXT_DEP=1 OS=WINDOWS

# 编译 hpipm 共享库
cd $HPIPM_DIR
make -j4 shared_library TARGET=GENERIC USE_C99_MATH=1 EXT_DEP=1 OS=WINDOWS

# Windows: 复制 .so 为 .dll
cp -f lib/libhpipm.so lib/libhpipm.dll
```

产物：
- `third_party/blasfeo/lib/libblasfeo.a`
- `third_party/hpipm/lib/libhpipm.so`
- `third_party/hpipm/lib/libhpipm.dll`（Windows Python ctypes 用）

---

## 8. 调用逻辑（单步 MPC）

```
control_RSS.m
  ├─ params = config()                          % 读取参数
  ├─ u_hat = zeros(3, K)                        % static init
  ├─ persistent: 加载 hpipm_qp_solver 模块      % 只一次, 防内存泄漏
  │
  └─ for m = 1:max_iter (3次)
       │
       ├─ qp = construct_complete_qp_from_rss(path, step, current_nu, state, u_hat, params)
       │     │
       │     ├─ 1. 构造 H, g (代价: 位置 + 姿态 + 控制正则 + RSS 强凸)
       │     ├─ 2. 构造 A, b (等式: 初始 3 + 递推 15 = 18)
       │     ├─ 3. 构造 Hq, gq, uq (二次约束: 轮速 24 + 转向锥 48 = 72)
       │     │     ├─ 轮速: ||H_n*nu(:,k)|| <= vimax
       │     │     └─ 转向锥: C = A - B - L <= 0 (两组 R 矩阵)
       │     └─ return qp_problem struct
       │
       ├─ cell -> 3D numpy 转换
       │
       ├─ result = py.hpipm_qp_solver.solve_qcqp(H, g, A, b, Hq, gq, uq)
       │     │
       │     ├─ 维度设置 (nv=36, ne=18, nb=0, ng=0, nq=72)
       │     ├─ balance 模式求解
       │     ├─ 失败则回退 robust 模式
       │     └─ return dict('x', 'status', 'obj_value', 'solve_time')
       │
       ├─ u = reshape(x(1:18), 3, K)            % 提取 u
       └─ u_hat = u                              % 更新 u_hat, 下次凸化使用
  │
  └─ 输出: new_state_dot = R(psi)*(state_dot + 1*u(:,1))
          velocity     = state_dot + u(:,1)
```

---

## 9. 关键设计点

1. **MATLAB + Python 分工**：MATLAB 端构造 QP 矩阵 + SQP 外层循环，Python 端仅负责 HPIPM 求解。每次 SQP 迭代有一次跨语言调用。

2. **SQP 无条件更新**：`u_hat = u` 无论求解成功失败都更新（论文 Alg.1 无 break-on-failure）。失败时 `u_sol = zeros(3,K)`。

3. **k1=1 硬编码**：`control_RSS.m` 硬编码 `k1=1`，覆盖 `config.m` 中的 `0.15`。对应论文 Alg.1 line 11: `velocity = state_dot + k1*u(:,1)`（增益为 1）。

4. **persistent 防内存泄漏**：Python 模块 `hpipm_qp_solver` 只 import 一次，避免重复 reload 导致 `libhpipm.dll` 内存泄漏（曾触发 WinError 8）。

5. **DLL 路径自动设置**：`hpipm_qp_solver.py` 在 import 时自动通过 `os.add_dll_directory` 把 `third_party/hpipm/lib/` 加入 DLL 搜索路径，无需手动配置 PATH。

6. **变量命名**：HPIPM C/Python wrapper 内部决策变量名为 `v`（与论文符号一致），返回 dict 的键名为 `'x'`（代码命名），MATLAB 端通过 `result{'x'}` 读取。

7. **硬约束子集**：本代码未启用 HPIPM 的 slack/box/一般线性约束（`nb=0, ng=0, ns=0`），退化为纯 QCQP。所有非线性约束都转为二次形式 `0.5*x'*Hq*x + gq'*x <= uq`。
