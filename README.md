# RSS: Robotic Steering Safety Controller

基于论文 *RSS: A Reactive Safety Shield for Robotic Systems* 的四轮全向底盘模型预测控制（MPC）仿真复现项目。通过 **git submodule** 以引用方式接入作者的三种算法仓库，加上本地维护的 active-set 实现，共对比四种控制器。

## 目录结构

```
RSS_V2/
├── matlab/                        # 仿真框架（主入口）
│   │
│   │ ── 主入口与编排 ──────────────────────────────────────
│   ├── main.m                     # 主入口, 按步骤编排 (Step 0-7)
│   ├── compare_algorithms.m       # [兼容包装] 旧入口, 转调 main
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
│   ├── run_paper_baseline_case.m  # 基线专用仿真 (含 exitflag 跟踪)
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
├── algorithms/                    # 本地算法实现
│   └── control_active_set.m       # active-set (fmincon active-set)
│
├── third_party/                   # 外部引用 (git submodule)
│   ├── RSS_proposed/              # → github.com/serendipitjx/RSS_proposed
│   │   └── control_RSS.m          #   proposed (CVX+ECOS, 3 次迭代)
│   ├── RSS_sqp/                   # → github.com/serendipitjx/RSS_sqp
│   │   └── control_RSS.m          #   e-LMPC (fmincon SQP, MaxIter=1)
│   └── RSS_fmincon/               # → github.com/serendipitjx/RSS_fmincon
│       └── control_RSS.m          #   interior-point (fmincon interior-point)
│
└── tests/                         # 单元测试
```

> 本地运行仿真后还会生成 `results/`（输出结果）、`scenario_bank/`（预生成场景）两个目录，已在 `.gitignore` 中排除，不入库。

## 算法引用架构（git submodule）

三种算法通过 git submodule 以**引用**方式接入，不复制代码：

| 算法 | 来源 | 接口 |
|---|---|---|
| proposed-3iter | `third_party/RSS_proposed/` | `[new_state_dot] = control_RSS(path, k, state_dot, state)` |
| e-lmpc | `third_party/RSS_sqp/` | `[new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state)` |
| interior-point | `third_party/RSS_fmincon/` | `[new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state, params)` |
| active-set | `algorithms/control_active_set.m` (本地) | `[u, worldVelocity, bodyVelocity] = control_active_set(path, k, state_dot, state, config)` |

`run_one_case.m` 按算法名动态切换 submodule 路径（`addpath`/`rmpath`），避免三个 `control_RSS.m` 和 `config.m` 同名冲突。

### 获取项目（含 submodule）

```bash
# 方式 1: clone 时带 --recursive
git clone --recursive <repo-url>

# 方式 2: 已 clone, 后续初始化 submodule
git submodule update --init --recursive
```

## 主入口流程

[main.m](matlab/main.m) 按下列步骤顺序编排：

```
Step 0  路径设置             → setup_paths()  (含 submodule 检查)
Step 1  生成 seed 列表        → (参数解析)
Step 2  选择算法              → (参数解析)
Step 3  检查断点续跑          → comparison_load()
Step 4  运行批量仿真          → run_batch_simulation()
                                └─ run_one_case()  按算法名分发到 submodule/本地
Step 5  保存最终结果          → save_all_artifacts()
Step 6  生成论文对比图        → plot_paper_comparison()
Step 7  打印 Table II         → print_table_ii()
```

## 快速开始

### 环境要求

- MATLAB（需 Optimization Toolbox 提供 fmincon；需 CVX + ECOS 用于 proposed 算法）
- git（用于 submodule 管理）

### 运行批量对比

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

### 已知限制

- **proposed 算法**（RSS_proposed submodule）只输出 `new_state_dot`，缺少 `u` 和 `solve_time`，相关 metrics 记为 NaN，Table II 中求解时间列显示 N/A
- 三个 submodule 的 `control_RSS.m` 接口各不相同，`run_one_case.m` 负责适配
- submodule 内部调用各自的 `config.m`（无参函数），与项目的 `defaultConfig.m` 独立

## 算法说明

| 算法 | 求解器 | 来源 | 特点 |
|---|---|---|---|
| proposed-3iter | CVX + ECOS | RSS_proposed submodule | RSS 凸化 + 3 次迭代，K=5 |
| e-lmpc | fmincon SQP | RSS_sqp submodule | MaxIter=1，K=6 |
| interior-point | fmincon interior-point | RSS_fmincon submodule | K=6，可违反转向锥约束 |
| active-set | fmincon active-set | 本地 `algorithms/` | K=6，对齐 RSS_fmincon |

## 模块职责（速查）

| 模块 | 职责 |
|---|---|
| `main.m` | 主入口，按 Step 0-7 编排全流程 |
| `run_one_case.m` | 单次闭环仿真，按算法名分发到 submodule 或本地 |
| `comparison_init/load/save.m` | comparison 结构体的初始化 / 加载 / 增量保存 |
| `run_batch_simulation.m` | 核心 (seed × alg) 双层循环 |
| `setup_paths.m` | 路径设置 + submodule 初始化检查 |
| `pair_is_completed.m` | 判断 (seed, alg) 是否已成功 |
| `get_alg_results.m` | 按算法名提取所有 results |
| `plot_one_algorithm.m` | 单算法 summary 图 |
| `print_table_ii.m` | Table II 风格汇总表打印 |

## 关键约束

- submodule 以引用方式接入，不修改 submodule 内的代码
- 算法实现细节与原仓库完全一致
- active-set 为本地实现（原仓库无此版本）
- checkpoint 文件不可删除，用于断点续跑
