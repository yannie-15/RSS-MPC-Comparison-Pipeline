# RSS: Robotic Steering Safety Controller

基于论文 *RSS: A Reactive Safety Shield for Robotic Systems* 的四轮全向底盘模型预测控制（MPC）仿真复现项目。通过 **git submodule** 以引用方式接入作者的三种算法仓库，加上本地维护的 active-set 实现，共对比四种控制器。

## 目录结构

```
RSS_V2/
├── main.py                        # Python 入口 (推荐), 通过 MATLAB Engine 调用
├── main.m                         # MATLAB 兼容入口 (调用 run_closed_loop)
├── run_one_case.m                 # JSON 接口: 读取 config.json → run_closed_loop → 返回 summary JSON
├── run_closed_loop.m              # 从原 main.m 抽取的 100 步闭环函数, 返回结构化 result
├── default_python_config.json     # Python 默认配置
│
├── python/                        # Python 桥接层
│   ├── __init__.py
│   ├── matlab_bridge.py           # MATLAB Engine 启动/调用/关闭
│   ├── config_io.py               # JSON 配置加载/合并/保存
│   └── result_io.py               # summary JSON 解析/打印
│
├── tests/                         # 测试
│   ├── baseline/
│   │   ├── freeze_baseline.m      # Phase 0: 冻结原 main.m 基线
│   │   └── .gitkeep
│   ├── matlab/
│   │   ├── test_run_one_case.m    # MATLAB smoke test
│   │   └── test_parity.m          # 新旧结果回归对比
│   └── python/
│       └── test_smoke.py          # Python smoke test (不需 MATLAB Engine)
│
├── matlab/                        # 批量对比仿真框架
│   │
│   │ ── 主入口与编排 ──────────────────────────────────────
│   ├── main.m                     # 批量仿真主入口, 按步骤编排 (Step 0-7)
│   ├── paper_reproduction.m       # 论文 Section IV 复现独立脚本
│   │
│   │ ── 状态管理 ────────────────────────────────────────
│   ├── comparison_init.m          # 初始化 comparison 结构体
│   ├── comparison_load.m          # 加载 / 迁移 checkpoint (兼容旧格式)
│   ├── comparison_save.m          # 增量保存 checkpoint
│   │
│   │ ── 核心仿真 ────────────────────────────────────────
│   ├── run_batch_simulation.m     # (seed × algorithm) 批量仿真循环
│   ├── run_one_case.m             # 单场景闭环仿真, 按算法名分发到 submodule
│   ├── run_paper_baseline_case.m  # 论文复现专用仿真 (含 iter_num 跟踪与解有限性检查, 4 种算法)
│   │
│   │ ── 查询辅助 ────────────────────────────────────────
│   ├── pair_is_completed.m        # 判断 (seed, alg) 是否已成功
│   ├── find_pair_index.m          # 查找 (seed, alg) 在 completedPairs 中的行号
│   ├── get_alg_results.m          # 提取某算法的所有 results
│   ├── safe_metric.m              # 安全读取 metrics 字段 (缺失/非数值 → NaN)
│   │
│   │ ── 场景与参考轨迹 ──────────────────────────────────
│   ├── scenario_bank.m            # 场景库 (每个 seed 一个 .mat, 可复现)
│   ├── scenario_generator.m       # Latin 超立方采样
│   ├── generateReference.m        # Bernstein 多项式参考轨迹
│   ├── defaultConfig.m            # 默认参数 (论文固定参数)
│   │
│   │ ── 动力学与指标 ────────────────────────────────────
│   ├── propagateState.m           # 状态传播
│   ├── computeWheelOutputs.m      # 轮速 / 轮角计算
│   ├── computeMetrics.m           # RMSE / J / 求解时间 / 约束违反率
│   ├── verify_constraints_RSS.m   # SQP约束验证 (每次子迭代原始约束重检)
│   │
│   │ ── 可视化与导出 ────────────────────────────────────
│   ├── plot_one_algorithm.m       # 单算法 summary 图
│   ├── save_algorithm_csv.m       # 单算法 CSV 导出
│   ├── print_table_ii.m           # Table II 汇总打印
│   ├── plot_paper_comparison.m    # 论文风格对比图 (Fig.3/4/5)
│   ├── replot_per_seed.m          # 逐 seed 轨迹重绘
│   │
│   │ ── 工具 ────────────────────────────────────────────
│   └── setup_paths.m              # 路径设置 (含 submodule 检查)
│
├── algorithms/                    # 所有算法 (统一目录)
│   │
│   │ ── git submodule (引用外部仓库, 不复制代码) ────────
│   ├── RSS_proposed/              # → github.com/serendipitjx/RSS_proposed
│   │   └── control_RSS.m          #   proposed (CVX+SDPT3, 3 次迭代)
│   ├── RSS_sqp/                   # → github.com/serendipitjx/RSS_sqp
│   │   └── control_RSS.m          #   e-LMPC (fmincon SQP, MaxIter=1)
│   ├── RSS_fmincon/               # → github.com/serendipitjx/RSS_fmincon
│   │   └── control_RSS.m          #   interior-point (fmincon interior-point)
│   │
│   │ ── 本地实现 (原仓库无此版本) ───────────────────────
│   └── control_active_set.m       # active-set (fmincon active-set)
│
├── third_party/                   # 第三方求解器源码 (参考)
│   ├── blasfeo/                   # BLASFEO 线性代数库
│   └── hpipm/                     # HPIPM QP 求解器
```

> 本地运行仿真后还会生成 `results/`（输出结果）、`scenario_bank/`（预生成场景）两个目录。

## 算法引用架构（git submodule）

四种算法**统一放在 `algorithms/` 目录**下，其中三种以 git submodule 引用外部仓库，active-set 为本地实现：

| 算法 | 来源 | 接口 |
|---|---|---|
| proposed-3iter | `algorithms/RSS_proposed/` (submodule) | `[u, new_state_dot, velocity, diagnostics] = control_RSS(path, k, state_dot, state)` |
| e-lmpc | `algorithms/RSS_sqp/` (submodule) | `[new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state)` |
| interior-point | `algorithms/RSS_fmincon/` (submodule) | `[new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state, params)` |
| active-set | `algorithms/control_active_set.m` (本地) | `[u, worldVelocity, bodyVelocity] = control_active_set(path, k, state_dot, state, config)` |

> **Phase 3 修复**：`control_RSS.m`（RSS_proposed）新增第 4 输出 `diagnostics`（可选，旧调用不受影响）。求解失败时不再无条件执行 `u_hat = u`，改为跟踪 `u_last_feasible` 并在失败时回退到最后可行解。求解失败后停止当前步后续 SQP 迭代。
>
> **Phase 6 修复**：日志解析支持 SDPT3 格式（`Total CPU time (secs) = x.xx`），解析失败记为 `NaN`（不再记为 `0`）。移除了 `pause(0.01)` 延迟。

`run_one_case.m` 按算法名动态切换 submodule 路径（`addpath`/`rmpath`），避免三个 `control_RSS.m` 和 `config.m` 同名冲突。`run_paper_baseline_case.m`（论文复现专用）改为在循环外一次性 `addpath` 当前算法的 submodule 并用 `onCleanup` 注册 `rmpath`，循环内不再切换路径，从根本上避免多版本 `control_RSS` 函数缓存冲突（此前 `clear functions` 不足以解决 e-lmpc/interior-point 的"输入/输出参数太多"错误）。

### 获取项目（含 submodule）

```bash
# 方式 1: clone 时带 --recursive
git clone --recursive <repo-url>

# 方式 2: 已 clone, 后续初始化 submodule
git submodule update --init --recursive
```

## 命令速查

三种入口，按使用场景选择：

| 场景 | 入口 | 命令 |
|------|------|------|
| **单次 proposed 仿真（推荐）** | Python | `python main.py --config default_python_config.json` |
| **单次 proposed 仿真（MATLAB）** | 根目录 `main.m` | `cd D:\PROJECT\RSS_V2; main` |
| **JSON 接口（MATLAB）** | `run_one_case.m` | `summary_json = run_one_case('default_python_config.json')` |
| **4 算法批量对比** | `matlab/main.m` | `cd matlab; setup_paths; main` |
| **论文 Section IV 复现** | `paper_reproduction.m` | `cd matlab; setup_paths; paper_reproduction` |
| **约束验证** | `verify_constraints_RSS.m` | `cd matlab; setup_paths; verify_constraints_RSS` |
| **Python smoke test** | `tests/python/` | `python tests/python/test_smoke.py` |
| **MATLAB smoke test** | `tests/matlab/` | `cd D:\PROJECT\RSS_V2; test_run_one_case` |
| **新旧结果回归** | `tests/matlab/` | 先 `freeze_baseline`，再 `test_parity` |

> **注意**：根目录 `main.m` 和 `matlab/main.m` 是两个不同的入口。根目录 `main.m` 调用 `run_closed_loop` 跑单次 proposed-3iter 仿真；`matlab/main.m` 是 4 算法批量对比编排器。

## 快速开始

### 方式一：Python 入口（推荐）

通过 `main.py` 调用 MATLAB Engine，一次性运行完整 100 步闭环仿真。Python 负责配置和结果汇总，MATLAB 负责全部 RSS 数学、CVX+SDPT3 求解和闭环仿真。

#### 安装 MATLAB Engine for Python

```bash
# 1. 在 MATLAB 中获取 matlabroot
#    >> matlabroot
# 2. 安装 MATLAB Engine
cd <matlabroot>/extern/engines/python
python setup.py install

# 3. 验证安装
python -c "import matlab.engine; print('OK')"
```

> MATLAB R2026a 支持 Python 3.9-3.12（以本机 `matlabroot/extern/engines/python/setup.py` 为准）。

#### 运行

```bash
# 使用默认配置 (100步, SDPT3, 无实时绘图)
python main.py --config default_python_config.json

# 自定义 seed 和输出目录
python main.py --seed 3 --output results/seed_0003 --no-plot

# 指定 case ID
python main.py --case-id test_001 --output results/test_001

# 开启实时绘图
python main.py --live-plot

# 查看所有参数
python main.py --help
```

#### CLI 参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--config` | JSON 配置文件路径 | `default_python_config.json` |
| `--case-id` | 案例标识符 | `paper_fixed_001` |
| `--seed` | 随机种子 | `1` |
| `--output` | 输出目录 | `results/<case_id>` |
| `--live-plot` | 开启 MATLAB 实时绘图 | 关闭 |
| `--no-plot` | 关闭绘图（默认） | — |
| `--save-figures` | 保存图片到输出目录 | 关闭 |
| `--force` | 覆盖已有输出目录 | 关闭 |

#### 退出码

| 退出码 | 含义 |
|--------|------|
| 0 | 成功（算法完成所有步） |
| 2 | 算法失败（求解器失败、约束违反等） |
| 3 | 基础设施错误（MATLAB Engine 未安装、配置错误等） |

#### `default_python_config.json` 字段

```json
{
  "case_id": "paper_fixed_001",
  "seed": 1,
  "trajectory_mode": "paper_fixed",
  "solver": "sdpt3",
  "num_steps": 100,
  "live_plot": false,
  "save_figures": false,
  "save_full_log": false,
  "output_dir": "results/paper_fixed_001"
}
```

| 字段 | 类型 | 说明 |
|------|------|------|
| `case_id` | string | 案例标识符，用于输出目录命名 |
| `seed` | int | 随机种子（未来随机场景用，当前固定轨迹不依赖） |
| `trajectory_mode` | string | 轨迹模式（`paper_fixed` = 论文固定贝塞尔轨迹） |
| `solver` | string | CVX 求解器（`sdpt3` 或 `ecos`，当前推荐 `sdpt3`） |
| `num_steps` | int | 仿真步数（默认 100，会同步 `num_path_pts`） |
| `live_plot` | bool | 是否开启 MATLAB 实时绘图 |
| `save_figures` | bool | 是否保存图片到输出目录 |
| `save_full_log` | bool | 是否保存完整日志 |
| `output_dir` | string | 输出目录路径 |

#### 输出

每次运行在 `output_dir` 中生成：
- `resolved_config.json` — 合并后的完整配置
- `result.mat` — 完整仿真数据（状态、速度、控制、轮速、轮角、求解时间等）
- `summary.json` — 仿真摘要（成功/失败、RMSE、代价、约束违反等）

#### Python smoke test

```bash
python tests/python/test_smoke.py
```

不需要 MATLAB Engine，仅测试 Python 模块导入和 JSON 解析。

### 方式二：MATLAB 入口（单次 proposed 仿真）

#### 兼容入口 `main.m`（根目录）

```matlab
cd D:\PROJECT\RSS_V2
main                    % 调用 run_closed_loop, 开实时绘图
```

#### 直接调用 `run_closed_loop`

```matlab
cd D:\PROJECT\RSS_V2
cfg = config();         % 从 algorithms/RSS_proposed/config.m 加载默认参数
cfg.live_plot = false;  % 关闭实时绘图 (批量运行推荐)
result = run_closed_loop(cfg);

% result 包含完整仿真数据:
%   result.state_history         - Nx3 状态 [x, y, psi]
%   result.body_velocity_history - Nx3 车体系速度
%   result.control_history       - Nx3 控制增量
%   result.wheel_speed_history   - Nx4 轮速
%   result.wheel_angle_history   - Nx4 轮角
%   result.steering_rate_history - Nx4 转向率
%   result.solver_time_history   - Nx1 每步求解时间
%   result.cvx_status_history    - {Nx1} CVX 状态
%   result.trajectory_cost       - 总代价 J
%   result.position_rmse         - 位置 RMSE
%   result.completed_steps       - 完成步数
%   result.success               - 是否成功
```

#### JSON 接口 `run_one_case`（与 Python 桥接）

```matlab
cd D:\PROJECT\RSS_V2
summary_json = run_one_case('default_python_config.json');
summary = jsondecode(summary_json);
disp(summary.success);
disp(summary.position_rmse);
```

> `run_one_case` 会自动创建 `output_dir`，保存 `result.mat` 和 `summary.json`。

#### MATLAB 测试

```matlab
cd D:\PROJECT\RSS_V2
test_run_one_case    % smoke test: run_closed_loop + run_one_case + diagnostics

% 新旧结果回归对比 (需先冻结基线)
cd tests/baseline
freeze_baseline      % 冻结原 main.m 基线 (约 20 分钟, 300 次 SDPT3)
cd D:\PROJECT\RSS_V2
test_parity          % 12 项检查: state/velocity/control/cost/RMSE/约束
```

### 方式三：4 算法批量对比（`matlab/main.m`）

批量对比 proposed-3iter / e-lmpc / interior-point / active-set 四种算法，支持多 seed。

#### 环境要求

- MATLAB（需 Optimization Toolbox 提供 fmincon；需 CVX + SDPT3 用于 proposed 算法）
- git（用于 submodule 管理）

> **注（ECOS → SDPT3 切换）**：原 RSS_proposed submodule 0121 分支 `control_RSS.m` 内部使用 `cvx_solver ECOS`，但 ECOS 2.0.10 在 MATLAB R2026a 上对此 problem 触发 segfault（访问冲突，MATLAB 整体崩溃）。已将 submodule 的 `control_RSS.m` 第34行改为 `cvx_solver SDPT3`，并在 [run_paper_baseline_case.m](matlab/run_paper_baseline_case.m) 中针对 proposed-3iter 算法添加全局 `cvx_solver sdpt3` 兜底 + 首步诊断输出（打印 control_RSS.m 实际路径、第34行内容、CVX 可用 solver 列表），便于定位缓存/遮蔽问题。同时修复了日志解析中 `runtime_line` 为空时直接索引 `{1}` 导致报错的问题（改为先判空再输出）。ECOS 长期方案待定（重新编译 ECOS 或降级 MATLAB）。

#### 运行批量对比

```matlab
cd('d:\PROJECT\RSS_V2\matlab'); setup_paths;

% 默认: 1:10 seeds, 4 种算法
main

% 自定义 seed 范围
main('seeds', 1:100)

% 只跑一种算法
main('algorithms', {'proposed-3iter'})

% 组合参数 + 强制重新生成场景
main('seeds', 1:50, 'algorithms', {'proposed-3iter'}, 'forceRegen', true)
```

#### 批量对比编排流程

`matlab/main.m` 按下列步骤顺序编排：

```
Step 0  路径设置             → setup_paths()  (含 submodule 检查)
Step 1  生成 seed 列表        → (参数解析)
Step 2  选择算法              → (参数解析)
Step 3  检查断点续跑          → comparison_load()
Step 4  运行批量仿真          → run_batch_simulation()
                                └─ run_one_case()  按算法名分发到 algorithms/ 下
Step 5  保存最终结果          → save_all_artifacts()
Step 6  生成论文对比图        → plot_paper_comparison()
Step 7  打印 Table II         → print_table_ii()
```

## 已知限制

- **proposed 算法**（RSS_proposed submodule）只输出 `new_state_dot`，缺少 `u` 和 `solve_time`，相关 metrics 记为 NaN，Table II 中求解时间列显示 N/A
- 三个 submodule 的 `control_RSS.m` 接口各不相同，`run_one_case.m` 负责适配
- submodule 内部调用各自的 `config.m`（无参函数），与项目的 `defaultConfig.m` 独立
- **main 批量仿真时**，`run_one_case.m` 通过 `setup_config_override` 用临时目录的 `config.m` 覆盖 submodule 的 `config.m`，使 proposed/e-lmpc 能用到 main 传入的随机 seed 场景参数（ctrl_pts/vimax/phidotmax/Lx/Ly 等）；interior-point/active-set 直接接收外部 config

## 论文结果复现

运行 [paper_reproduction.m](matlab/paper_reproduction.m) 可复现论文 Section IV 的固定轨迹实验。无需 seed，每个算法使用各自 submodule 的 config.m 参数：

```matlab
cd('d:\PROJECT\RSS_V2\matlab'); setup_paths;

paper_reproduction                                  % 默认: 重跑全部 4 种算法
paper_reproduction({'proposed-3iter'})              % 只重跑 proposed-3iter, 保留其余算法已有结果
paper_reproduction({'e-lmpc','active-set'})         % 只重跑指定算法, 保留其余
```

**增量复现模式**：默认不清空整个 `results/paper_reproduction/` 目录，而是：
- 加载已有的 `paper_reproduction.mat`，保留未重跑算法的旧结果
- 仅清理本次将重跑算法对应的 CSV / 汇总图 / 逐 seed 图目录
- 最终合并保存（Table II 汇总覆盖全部 4 种算法）

若需全量重跑，先删除 `results/paper_reproduction/paper_reproduction.mat` 或传入全部算法列表。

**参数来源**（每个算法用各自的 config.m，不在 paper_reproduction 中覆盖）：

| 算法 | 参数来源 |
|---|---|
| proposed-3iter | `algorithms/RSS_proposed/config.m` |
| e-lmpc | `algorithms/RSS_sqp/config.m` |
| interior-point | `algorithms/RSS_fmincon/config.m` |
| active-set | `matlab/defaultConfig.m`（本地，对齐 RSS_fmincon） |

**固定初始状态**（论文无随机性，所有算法统一）：

| 参数 | 值 | 适用算法 |
|---|---|---|
| 初始位姿 | `[0.05, 0.1, 0.2]` | 全部 4 种算法 |
| 初始速度 | `[0.01, 0.01, 0.01]` | 全部 4 种算法 |

> 注：RSS_proposed submodule 已切换到 `0121` 分支（commit d25d299 "proposed method 最终版"），其 `main.m` 初始条件为 `[0.05,0.1,0.2]+[0.01,0.01,0.01]`，与其余三种算法完全一致，无需特判。

**各 submodule config.m 参数对比**（submodule 保持原始参数，不修改）：

| 参数 | RSS_proposed (0121) | RSS_sqp | RSS_fmincon |
|---|---|---|---|
| Lx | 0.655 (m) | 0.655 (m) | 0.655 (m) |
| vimax | 5 (m/s) | 5 (m/s) | 5 (m/s) |
| phidotmax | 5π | 5π | 5π |
| t_end | 1.0 s | 1.0 s | 1.0 s |
| num_steps | 100 | 100 | 100 |
| ctrl_pts | 5 个点 | 5 个点 | 5 个点 |
| K (预测时域) | 6 | 6 | 6 |

> 注：三个 submodule 的 config.m 参数已天然一致（均为 m 单位）。RSS_proposed 的 `0121` 分支为最终版，与 RSS_sqp/RSS_fmincon 使用相同的几何/仿真参数。仅控制器增益 k1 在各算法 control_RSS.m 内部各自设定。

**与 main.m 的区别**：
- `main.m` 用 `scenario_bank(seed)` 生成随机场景，所有算法共用同一套参数（defaultConfig）
- `paper_reproduction.m` 每个算法用各自的 config.m 参数，仅初始状态固定

**输出**（写入 `results/paper_reproduction/`）：
- `per_algorithm/{算法名}/{算法名}_summary.png` — 单算法汇总图
- `per_algorithm/{算法名}/seed_0001.png` — 逐 seed 轨迹图
- `{算法名}_results.csv` — 单算法指标 CSV
- `paper_reproduction.mat` — 完整 comparison 结构体

**额外诊断**（相比 `run_one_case.m`，由 `run_paper_baseline_case.m` 提供）：
- 每步 `iter_num`（求解器迭代数，proposed/active-set 为 NaN）
- 解的有限性检查（NaN/Inf 标记失败步）
- warm-up 排除后的中位数/分位数耗时（P1-4 复现要求）

## 约束验证（verify_constraints_RSS）

[verify_constraints_RSS.m](matlab/verify_constraints_RSS.m) 用于验证 RSS proposed 算法在 SQP 求解过程中是否满足**原始问题约束**（不仅仅是最终输出轨迹）。

### 背景

RSS proposed 算法包含两类约束：

| 约束 | 类型 | CVX 中的处理 | 是否可能违反 |
|------|------|-------------|-------------|
| 轮速 `‖H·ν‖ ≤ vimax` | **凸**（SOC） | CVX 直接强制 | 求解成功 → 不可能违反 |
| 转向角速率 `|Δθ| ≤ dt·phidotmax` | **非凸** | SQP 线性化近似 | **可能违反**（线性化仅局部近似） |

现有的 `computeMetrics.m` 只检查最终执行轨迹（3 次 SQP 迭代后的输出），不检查每次子迭代（m=1,2,3）的中间解。本脚本填补这一盲区。

### 用法

```matlab
cd('D:\PROJECT\RSS_V2\matlab')
setup_paths;
verify_constraints_RSS
```

### 验证逻辑

1. 复现 `control_RSS.m` 的完整 CVX 问题（`quiet` 模式，与原代码一致）
2. 每次 CVX 求解后（100 步 × 3 次 = 300 次），立即用**原始非线性公式**重检：
   - 轮速：`norm(H{n} * nu(:,k))` vs `vimax`
   - 转角变化：`atan2(wv_k) - atan2(wv_{k-1})` vs `delta_theta = dt * phidotmax`
3. 按 SQP 迭代次数分组统计（iter1/iter2/iter3），观察收敛过程中违反是否减少
4. 输出汇总报告 + 保存 `constraint_verification_results.mat`

### 输出

- 命令窗口打印每次违反的详细信息（步号、迭代号、轮编号、超出量）
- 汇总报告：求解器失败率、轮速违反次数/最大超出量、转角违反次数/最大超出量
- 按 SQP 子迭代的分组分析
- `constraint_verification_results.mat`：完整验证数据（保存到 `results/` 目录）

## 仿真步数说明

论文复现（`paper_reproduction.m`）固定 `t_end=1s, dt=0.01s` → `num_steps=100`，且 `num_path_pts = num_steps`（路径离散点数与仿真步数一一对应，第 k 步跟踪路径第 k 个点）。四种算法面对完全相同的 100 个参考点和 100 步控制，确保结果可比性。

随机场景生成器 `scenario_generator.m` 的 `t_end ∈ [0.5, 2.0]` 随机，步数随之变为 50~200，用于压力测试算法在不同路径长度下的鲁棒性。

> 注：`num_path_pts = num_steps` 将路径离散化精度与仿真时长绑定。在论文复现场景下无影响，但若未来需独立调整路径精度或仿真时长，应解耦这两个参数。

## 算法说明

| 算法 | 求解器 | 来源 | 特点 |
|---|---|---|---|
| proposed-3iter | CVX + SDPT3 | RSS_proposed submodule | RSS 凸化 + 3 次迭代，K=6 |
| e-lmpc | fmincon SQP | RSS_sqp submodule | MaxIter=1，K=6 |
| interior-point | fmincon interior-point | RSS_fmincon submodule | K=6，可违反转向锥约束 |
| active-set | fmincon active-set | 本地 `algorithms/` | K=6，对齐 RSS_fmincon |

## 模块职责（速查）

| 模块 | 职责 |
|---|---|
| `main.py` | **Python 入口（推荐）**，CLI 参数解析 → MATLAB Engine → run_one_case |
| `run_closed_loop.m` | **根目录**，从原 main.m 抽取的 100 步闭环函数，返回结构化 result |
| `run_one_case.m` (根目录) | **JSON 接口**，读取 config.json → run_closed_loop → 返回 summary JSON |
| `main.m` (根目录) | MATLAB 兼容入口，调用 run_closed_loop |
| `python/matlab_bridge.py` | MATLAB Engine 启动/路径设置/调用 run_one_case/关闭 |
| `python/config_io.py` | JSON 配置加载/合并/CLI 覆盖/保存 |
| `python/result_io.py` | summary JSON 解析/打印 |
| `matlab/main.m` | 批量仿真主入口，按 Step 0-7 编排全流程 |
| `matlab/run_one_case.m` | 批量仿真单场景闭环，按算法名分发到 submodule |
| `matlab/comparison_init/load/save.m` | comparison 结构体的初始化 / 加载 / 增量保存 |
| `matlab/run_batch_simulation.m` | 核心 (seed × alg) 双层循环 |
| `matlab/setup_paths.m` | 路径设置 + submodule 初始化检查 |
| `matlab/pair_is_completed.m` | 判断 (seed, alg) 是否已成功 |
| `matlab/get_alg_results.m` | 按算法名提取所有 results |
| `matlab/plot_one_algorithm.m` | 单算法 summary 图 |
| `matlab/print_table_ii.m` | Table II 风格汇总表打印 |
| `tests/baseline/freeze_baseline.m` | 冻结原 main.m 基线数据 |
| `tests/matlab/test_run_one_case.m` | MATLAB smoke test |
| `tests/matlab/test_parity.m` | 新旧结果回归对比 |
| `tests/python/test_smoke.py` | Python smoke test（不需 MATLAB Engine） |

## Python → MATLAB 架构

```
用户/CI
  ↓
根目录 main.py
  ├── 读取 JSON 配置
  ├── 启动一次 MATLAB Engine
  ├── addpath(仓库根目录)
  ├── 调用 run_one_case(config_json_path)  ← 一次跨语言调用
  ├── 解析 MATLAB 返回的 summary JSON
  └── 根据结果设置退出码
            ↓
MATLAB run_one_case.m (根目录)
  ├── 合并 config() 默认值 + JSON 覆盖
  ├── 调用 run_closed_loop(cfg)
  │     ├── bezier_path 生成参考轨迹
  │     ├── 100 步完整闭环
  │     │   └── 每步调用 control_RSS (CVX + SDPT3)
  │     ├── 约束检查和指标计算
  │     └── 返回 result 结构体
  ├── 保存 result.mat 和 summary.json
  └── 返回 jsonencode(summary)
```

**关键原则**：Python 只调用一次 MATLAB，不逐 MPC step 调用。整个 100 步闭环在 MATLAB 内部完成，保证计时公平性和 CVX warm start 的完整性。

## 关键约束

- submodule 以引用方式接入，不修改 submodule 内的代码（RSS_proposed 的 control_RSS.m 例外：已修复失败回退和计时解析）
- 算法实现细节与原仓库完全一致（CVX 问题公式、代价函数、约束、K、rho 均未修改）
- active-set 为本地实现（原仓库无此版本）
- checkpoint 文件不可删除，用于断点续跑
- Python 入口不计算机器人状态，不读写 MATLAB Base Workspace，不解析 MATLAB 控制台文本
- MATLAB Engine 只启动一次，整个闭环只跨语言调用一次
