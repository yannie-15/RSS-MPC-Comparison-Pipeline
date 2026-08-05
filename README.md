# RSS-MPC-Comparison-Pipeline

基于论文 *RSS2026: Exploit Agile Mobility of Steerable-Wheeled Mobile Robots: A Fast Motion Planning Approach* 的四轮全向底盘 MPC 仿真复现项目，对比四种控制器。

## 目录结构

```
RSS-MPC-Comparison-Pipeline-rss_hpipm/
├── main.py                              # Python 入口, 通过 MATLAB Engine 调用
├── default_python_config.json           # Python 默认配置
│
├── algorithms/                          # 算法包
│   ├── RSS_proposed/                    # proposed 控制器 (普通目录, HPIPM dense QCQP)
│   │   ├── control_RSS.m                #   3 次 SQP 迭代, k1=1
│   │   ├── construct_complete_qp_from_rss.m  # QP 矩阵构造
│   │   ├── config.m                     #   算法参数
│   │   ├── hpipm_qp_solver.py           #   HPIPM Python 求解器接口
│   │   └── build_hpipm_windows.sh       #   Windows MSYS2 编译脚本
│   ├── RSS_sqp/                         # → github.com/serendipitjx/RSS_sqp (submodule, main)
│   ├── RSS_fmincon/                     # → github.com/serendipitjx/RSS_fmincon (submodule, main)
│   ├── RSS_active_set/                  # → github.com/serendipitjx/RSS_fmincon (submodule, active-set 分支)
│   └── .gitattributes
│
├── core/                                # 仿真核心工具
│   ├── setup_paths.m                    # 路径设置 (含 submodule 检查)
│   ├── defaultConfig.m                  # 默认参数
│   ├── generateReference.m              # Bernstein 多项式参考轨迹
│   ├── propagateState.m                 # 状态传播
│   ├── computeWheelOutputs.m            # 轮速 / 轮角计算
│   └── computeMetrics.m                 # RMSE / J / 求解时间 / 约束违反率
│
├── batch_simulation/                    # 批量仿真 (多 seed 随机场景)
│   ├── main.m                           # MATLAB 批量仿真主入口
│   ├── run_batch_simulation.m           # (seed × algorithm) 批量仿真循环
│   ├── run_one_case.m                   # 单场景闭环仿真, 按算法名分发到 submodule
│   ├── comparison_init.m / _load.m / _save.m  # comparison 结构体管理
│   ├── scenario_bank.m / scenario_generator.m # 场景库与采样
│   ├── plot_one_algorithm.m             # 单算法 summary 图
│   ├── plot_paper_comparison.m          # 论文风格对比图 (Fig.3/4/5)
│   ├── replot_per_seed.m                # 逐 seed 轨迹重绘
│   ├── print_table_ii.m                 # Table II 汇总打印
│   ├── save_algorithm_csv.m             # 单算法 CSV 导出
│   ├── matlab_bridge.py                 # MATLAB Engine 启动/调用/关闭
│   ├── config_io.py                     # JSON 配置加载/合并/保存
│   └── result_io.py                     # summary JSON 解析/打印
│
├── paper_reproduction/                  # 论文 Section IV 复现
│   ├── paper_reproduction.m             # 复现入口
│   ├── run_paper_baseline_case.m        # 论文复现专用仿真
│   └── results/                         # 输出目录 (运行时生成)
│
├── verification/                        # HPIPM 约束验证
│   ├── verify_constraints_hpipm.m       # HPIPM 解约束验证
│   └── results/                         # 输出目录 (运行时生成)
│
└── third_party/                         # 第三方求解器源码
    ├── blasfeo/                         # BLASFEO 线性代数库 (submodule)
    └── hpipm/                           # HPIPM QP 求解器
```

## 算法说明

| 算法 | 求解器 | 来源 | 特点 |
|---|---|---|---|
| proposed-3iter | HPIPM (dense QCQP) | RSS_proposed (普通目录) | RSS 凸化 + 3 次 SQP 迭代, K=6 |
| e-lmpc | fmincon SQP | RSS_sqp (submodule) | MaxIter=1, K=6 |
| interior-point | fmincon interior-point | RSS_fmincon (submodule) | K=6 |
| active-set | fmincon active-set | RSS_active_set (submodule) | K=6 |

## 获取项目

```bash
# clone 时带 --recursive (自动初始化全部 submodule)
git clone --recursive <repo-url>

# 或已 clone 后初始化 submodule
git submodule update --init --recursive
```

## 命令速查

| 场景 | 入口 | 命令 |
|------|------|------|
| 4 算法批量对比 (Python) | `main.py` | `python main.py` |
| 4 算法批量对比 (MATLAB) | `batch_simulation/main.m` | `cd batch_simulation; main` |
| 论文 Section IV 复现 | `paper_reproduction/paper_reproduction.m` | `cd paper_reproduction; paper_reproduction` |
| 约束验证 | `verification/verify_constraints_hpipm.m` | `cd verification; verify_constraints_hpipm` |

## Python 入口

通过 `main.py` 调用 MATLAB Engine 执行批量仿真。

### 安装 MATLAB Engine for Python

```bash
cd <matlabroot>/extern/engines/python
python setup.py install
python -c "import matlab.engine; print('OK')"
```

### 运行

```bash
# 默认: 1:10 seeds, 全部 4 种算法
python main.py

# 自定义 seed 范围
python main.py --seeds 1:100

# 只跑一种算法
python main.py --algorithms proposed-3iter

# 组合参数 + 强制重新生成场景
python main.py --seeds 1:50 --algorithms proposed-3iter,e-lmpc --force-regen
```

### CLI 参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--seeds` | 种子范围，格式 `start:end` 或单个整数 | `1:10` |
| `--algorithms` | 算法列表，逗号分隔 | `e-lmpc,active-set,interior-point,proposed-3iter` |
| `--force-regen` | 强制重新生成场景文件 | 关闭 |

### 输出

运行后在 `results/batch/` 生成：
- `comparison_checkpoint.mat` / `comparison_final.mat`
- `csv/{算法名}_results.csv`
- `per_algorithm/{算法名}_summary.png` + `per_algorithm/{算法名}/seed_XXXX.png`
- `comparison/fig3_trajectories.png` 等论文对比图

## MATLAB 入口

### 环境要求

- MATLAB (需 Optimization Toolbox 提供 fmincon)
- Python 3.9+ (HPIPM 接口, 仅 proposed-3iter 需要)
- MSYS2 + gcc + make + bc (编译 HPIPM/BLASFEO, 仅 Windows 需要)

### HPIPM 编译 (仅 proposed-3iter 需要)

```powershell
# 1. MSYS2 UCRT64 安装工具链
pacboy sync:mman-git ucrt64/toolchain msys/make msys/bc

# 2. 运行编译脚本
C:\msys64\usr\bin\env.exe MSYSTEM=UCRT64 /usr/bin/bash -lc "/d/PROJECT/RSS-MPC-Comparison-Pipeline-rss_hpipm/algorithms/RSS_proposed/build_hpipm_windows.sh"

# 3. 验证 Python 加载
python -c "import sys; sys.path.insert(0,'algorithms/RSS_proposed'); import hpipm_qp_solver; print('HPIPM_OK=', hpipm_qp_solver._HPIPM_OK)"
```

编译产物：`third_party/hpipm/lib/libhpipm.dll`、`third_party/blasfeo/lib/libblasfeo.a`

### 运行批量对比

```matlab
cd('batch_simulation');

main                                          % 默认: 1:10 seeds, 4 种算法
main('seeds', 1:100)                          % 自定义 seed 范围
main('algorithms', {'proposed-3iter'})        % 只跑一种算法
main('seeds', 1:50, 'algorithms', {'proposed-3iter'}, 'forceRegen', true)
```

## 论文复现

```matlab
cd('paper_reproduction');

paper_reproduction                                  % 重跑全部 4 种算法
paper_reproduction({'proposed-3iter'})              % 只重跑指定算法, 保留其余
paper_reproduction({'e-lmpc','active-set'})         % 只重跑指定算法, 保留其余
```

输出写入 `paper_reproduction/results/`：
- `per_algorithm/{算法名}/{算法名}_summary.png`
- `csv/{算法名}_results.csv`
- `paper_reproduction.mat`

## 约束验证

```matlab
cd('verification');
verify_constraints_hpipm
```

验证 HPIPM 求解的解是否满足原始非线性约束（轮速 SOC 约束、转向角速率约束）。输出汇总报告 + `verification/results/hpipm_constraint_verification.mat`。

---

## RSS_proposed 算法文件详解

### 文件总览

| 文件 | 作用 |
|------|------|
| `control_RSS.m` | 算法主入口，SQP 外层循环 + 调用 Python 求解器 |
| `construct_complete_qp_from_rss.m` | QP/QCQP 矩阵显式构造（含 Bug1-5 修复） |
| `hpipm_qp_solver.py` | Python 端 HPIPM dense QCQP 求解器封装 |
| `config.m` | 算法参数（机器人几何、控制器增益、仿真参数、路径控制点） |
| `build_hpipm_windows.sh` | Windows MSYS2 编译 blasfeo + hpipm 共享库脚本 |
| `.gitignore` | 忽略 MATLAB 自动备份与验证产物 |

### 调用逻辑（单步 MPC）

```
control_RSS.m
  ├─ 1. 读取 config() 参数
  ├─ 2. 初始化 persistent Python 环境 (首次调用时)
  │     └─ sys.path.append(当前目录) + importlib.reload(hpipm_qp_solver)
  ├─ 3. SQP 外层循环 (m = 1..3, 无条件更新 u_hat)
  │     │
  │     ├─ 3.1 调用 construct_complete_qp_from_rss(path, step, current_nu, state, u_hat, params)
  │     │       │
  │     │       ├─ 构造 Hessian H (n_var × n_var)
  │     │       │     位置代价 w_pos=30 + 姿态代价 w_psi=1 + 控制正则 w_control=0.3 + RSS 松弛 rho=0.01
  │     │       ├─ 构造等式约束 A*x = b (动力学递推, 3K 条)
  │     │       ├─ 构造二次约束 Hq/gq/uq
  │     │       │     ├─ 轮速 SOC 约束 (4 轮 × K 步 = 24 条)
  │     │       │     └─ 转向锥凸化约束 (两组 R 矩阵 × 4 轮 × K 步 = 48 条)
  │     │       └─ 返回 qp 结构体 (H, g, A, b, Hq, gq, uq, 元数据)
  │     │
  │     ├─ 3.2 cell → 3D numpy 数组转换
  │     │       Hq cell (n_var × n_var) → Hq_3d (n_var × n_var × n_qcqp)
  │     │       gq cell (n_var × 1)    → gq_2d (n_var × n_qcqp)
  │     │
  │     ├─ 3.3 调用 Python 求解器
  │     │       result = py.hpipm_qp_solver.solve_qcqp(H, g, A, b, Hq_3d, gq_2d, uq, verbose=False)
  │     │       │
  │     │       └─ [Python 内部] hpipm_qp_solver.solve_qcqp
  │     │             ├─ 数组类型统一 (np.asarray float64)
  │     │             ├─ 设置 hpipm_dense_qcqp_dim (nv, ne, nq)
  │     │             ├─ 设置 hpipm_dense_qcqp (H, g, A, b, Hq, gq, uq)
  │     │             ├─ 求解: balance 模式 + tol=1e-6 + iter_max=300
  │     │             ├─ 若失败回退: robust 模式 + iter_max=500
  │     │             └─ 返回 {x, status, obj_value, solve_time, iters}
  │     │
  │     ├─ 3.4 提取解 x (n_var,), 取前 3K 维 reshape 为 u_sol (3 × K)
  │     ├─ 3.5 记录诊断 (status, optval, solve_time, solver_name)
  │     └─ 3.6 无条件更新 u_hat = u_sol (与原始 CVX 0121 分支一致)
  │
  ├─ 4. 输出状态导数 new_state_dot = R(psi) * (current_nu + k1 * u(:,1))  [k1=1 硬编码]
  ├─ 5. 输出车体系速度 velocity = current_nu + u(:,1)
  └─ 6. 汇总 diagnostics.total_solve_time
```

### 决策变量与约束规模

| 项目 | 维度 | 说明 |
|------|------|------|
| 决策变量 `x` | 36 | `x = [u(:); nu(:)]`，u 为 3×K 控制增量，nu 为 3×K 车体系速度 |
| 等式约束 | 18 | 动力学递推：`nu(:,1)=current_nu+u(:,1)`，`nu(:,k+1)=nu(:,k)+u(:,k+1)` |
| 二次约束 | 24 | 轮速 SOC：`‖H_n · nu(:,k)‖ ≤ vimax`，4 轮 × 6 步 |
| 二次约束 | 48 | 转向锥凸化：两组 R 矩阵线性化，4 轮 × 6 步 × 2 |
| 预测时域 K | 6 | 固定 |

### 各文件目的

#### `control_RSS.m` — 算法主入口
- **接口**：`[u, new_state_dot, velocity, diagnostics] = control_RSS(path, step, state_dot, state)`
- **职责**：
  1. 加载参数（`K=6, rho=0.01, k1=1, max_iter=3`）
  2. 初始化 Python 环境（`persistent` 变量保证只设置一次 `sys.path` 与一次 `importlib.reload`）
  3. SQP 外层循环 3 次迭代：每次调用 `construct_complete_qp_from_rss` 构造 QP → 调用 `py.hpipm_qp_solver.solve_qcqp` 求解
  4. **关键策略**：无论求解成功失败都更新 `u_hat = u_sol`，不 break、不回退（与原始 CVX 0121 分支完全一致）
  5. 输出控制增量 `u`、世界系状态导数、车体系速度、诊断结构体

#### `construct_complete_qp_from_rss.m` — QP 矩阵构造
- **接口**：`qp = construct_complete_qp_from_rss(path, step, current_nu, state, u_hat, params)`
- **职责**：将 RSS 原始非凸问题显式构造为 dense QCQP 标准形式
  - **目标函数** `0.5·x'Hx + g'x`：
    - 位置跟踪代价（`w_pos=30`，k=2..K）
    - 姿态跟踪代价（`w_psi=1`，k=1..K）
    - 控制正则化（`w_control=0.3`，`‖u‖²`）
    - RSS 松弛项（`rho=0.01`，`‖u-u_hat‖²`）
  - **等式约束** `A·x = b`：动力学递推（3K 条）
  - **二次约束** `0.5·x'Hq_i·x + gq_i'·x ≤ uq_i`：
    - 轮速 SOC 约束（4×K=24 条）
    - 转向锥凸化约束（两组 R 矩阵，2×4×K=48 条）
- **Bug 修复记录**（5 个索引/系数 Bug，与原 CVX 公式对齐）：
  - Bug1：nu 索引偏移从 +1/+2/+3 改为 +0/+1/+2
  - Bug2：加入跨预测阶段的 Hessian 交叉项
  - Bug3：位置代价常数部分加入 `current_nu(1:2)·dt` 位移
  - Bug4：Hessian 系数改为 `2·w·dt²`（适配 `0.5·x'Hx` 形式）
  - Bug5：等式约束递推段索引从 `+i` 改为 `+(i-1)`
- **返回**：结构体 `qp`（含 H, g, A, b, Hq cell, gq cell, uq, 元数据）

#### `hpipm_qp_solver.py` — Python 求解器封装
- **接口**：`result = solve_qcqp(H, g, A, b, Hq, gq, uq, verbose=False)`
- **职责**：仅负责求解，不构造 QP
  1. import 时自动设置 `libhpipm.dll` 搜索路径（`os.add_dll_directory` + PATH fallback）
  2. 加载 HPIPM Python wrapper（`hpipm_dense_qcqp_dim/qcqp/sol/solver_arg/solver`）
  3. 将 MATLAB 传入的 numpy 数组统一转为 `float64`
  4. 设置维度 `nv/ne/nb/ng/nq`，填入 QP 数据
  5. 求解：`balance` 模式 + `tol=1e-6` + `iter_max=300`
  6. **失败回退**：`robust` 模式 + `iter_max=500`
  7. 返回 dict：`{x, status, status_str, obj_value, solve_time, iters}`

#### `config.m` — 算法参数
- **机器人几何**：`Lx=0.655m, Ly=0.335m`，4 轮位置（`wheel_pos`）
- **约束上限**：`vimax=5 m/s`（最大轮速），`phidotmax=5π rad/s`（最大转向角速率）
- **控制器增益**：`k1=k2=0.15, k3=0.1, eps=0.001`（注：`control_RSS.m` 中 `k1` 硬编码为 1，覆盖此值）
- **仿真参数**：`dt=0.01s, t_end=1s, num_steps=100`
- **参考路径**：5 个 Bernstein 控制点（生成贝塞尔曲线参考轨迹）

#### `build_hpipm_windows.sh` — 编译脚本
- **目的**：在 Windows MSYS2 环境编译 `blasfeo` 静态库 + `hpipm` 共享库
- **前置**：MSYS2 UCRT64 + gcc + make + bc
- **产物**：
  - `third_party/blasfeo/lib/libblasfeo.a`
  - `third_party/hpipm/lib/libhpipm.so` 与 `libhpipm.dll`（Windows Python ctypes 加载需 `.dll` 扩展名）
- **用法**：`C:\msys64\usr\bin\env.exe MSYSTEM=UCRT64 /usr/bin/bash -lc "<脚本路径>"`

#### `.gitignore`
- 忽略 MATLAB 自动备份（`*.asv, *.m~`）
- 忽略约束验证产物（`constraint_verification_results.mat`）

### 关键设计点

1. **MATLAB + Python 分工**：MATLAB 负责构造 QP 矩阵（修复 Bug1-5）与 SQP 外层循环；Python 仅负责求解（`balance` + `robust` 回退）
2. **SQP 无条件更新**：与原始 CVX 0121 分支一致，无论成功失败都更新 `u_hat`，不 break、不回退
3. **k1=1 硬编码**：`control_RSS.m` 第 25 行 `k1=1`，覆盖 `config.m` 中的 `params.k1=0.15`，匹配原 CVX 0121 分支行为
4. **persistent 环境初始化**：避免每次调用都重新加载 Python 模块，防止 `libhpipm.dll` 内存泄漏（曾导致 WinError 8）
5. **DLL 路径自动设置**：`hpipm_qp_solver.py` 在 import 时自动通过 `os.add_dll_directory` 添加 `third_party/hpipm/lib/` 到搜索路径
