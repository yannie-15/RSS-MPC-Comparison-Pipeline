# RSS2026: Exploit Agile Mobility of Steerable-Wheeled Mobile Robots: A Fast Motion Planning Approach

基于论文 *RSS2026: Exploit Agile Mobility of Steerable-Wheeled Mobile Robots: A Fast Motion Planning Approach* 的四轮全向底盘模型预测控制（MPC）仿真复现项目。通过 **git submodule** 以引用方式接入作者的算法仓库（含 active-set 分支），共对比四种控制器。

## 目录结构

按功能划分，不再按语言（MATLAB/Python）分目录。

```
RSS_V2/
├── main.py                              # Python 入口 (推荐), 通过 MATLAB Engine 调用
├── default_python_config.json           # Python 默认配置
│
├── algorithms/                          # 算法包
│   │ ── RSS_proposed 为普通目录 (已内嵌 HPIPM 版), 其余为 git submodule ──
│   ├── RSS_proposed/                    # 普通目录 (原 github.com/serendipitjx/RSS_proposed 0121 分支)
│   │   ├── control_RSS.m                #   proposed 控制器 (HPIPM dense QCQP, 3 次 SQP 迭代, k1=1)
│   │   ├── construct_complete_qp_from_rss.m  # QP 矩阵构造 (修复 Bug1-5)
│   │   ├── config.m                     #   算法参数
│   │   ├── hpipm_qp_solver.py           #   HPIPM dense QCQP 求解器 (Python 接口, 含 DLL 路径自动设置)
│   │   └── build_hpipm_windows.sh       #   Windows MSYS2 编译脚本
│   ├── RSS_sqp/                         # → github.com/serendipitjx/RSS_sqp (main 分支)
│   │   └── control_RSS.m                #   e-LMPC (fmincon SQP, MaxIter=1)
│   ├── RSS_fmincon/                     # → github.com/serendipitjx/RSS_fmincon (main 分支)
│   │   └── control_RSS.m                #   interior-point (fmincon interior-point)
│   ├── RSS_active_set/                  # → github.com/serendipitjx/RSS_fmincon (active-set 分支)
│   │   └── control_RSS.m                #   active-set (fmincon active-set)
│   └── .gitattributes
│
├── core/                                # 仿真核心工具 (跨功能共享)
│   ├── setup_paths.m                    # 路径设置 (含 submodule 检查)
│   ├── defaultConfig.m                  # 默认参数 (论文固定参数)
│   ├── generateReference.m              # Bernstein 多项式参考轨迹
│   ├── propagateState.m                 # 状态传播
│   ├── computeWheelOutputs.m            # 轮速 / 轮角计算
│   └── computeMetrics.m                 # RMSE / J / 求解时间 / 约束违反率
│
├── batch_simulation/                    # 批量仿真 (多 seed 随机场景)
│   │ ── MATLAB ──
│   ├── main.m                           # 批量仿真主入口, 按步骤编排 (Step 0-7)
│   ├── run_batch_simulation.m           # (seed × algorithm) 批量仿真循环
│   ├── run_one_case.m                   # 单场景闭环仿真, 按算法名分发到 submodule
│   ├── comparison_init.m                # 初始化 comparison 结构体
│   ├── comparison_load.m                # 加载 / 迁移 checkpoint (兼容旧格式)
│   ├── comparison_save.m                # 增量保存 checkpoint
│   ├── scenario_bank.m                  # 场景库 (每个 seed 一个 .mat, 可复现)
│   ├── scenario_generator.m             # Latin 超立方采样
│   ├── pair_is_completed.m              # 判断 (seed, alg) 是否已成功
│   ├── find_pair_index.m                # 查找 (seed, alg) 在 completedPairs 中的行号
│   ├── get_alg_results.m                # 提取某算法的所有 results
│   ├── safe_metric.m                    # 安全读取 metrics 字段 (缺失/非数值 → NaN)
│   ├── plot_one_algorithm.m             # 单算法 summary 图
│   ├── plot_paper_comparison.m          # 论文风格对比图 (Fig.3/4/5, 含 boxplot 回退)
│   ├── replot_per_seed.m                # 逐 seed 轨迹重绘
│   ├── print_table_ii.m                 # Table II 汇总打印
│   ├── save_algorithm_csv.m             # 单算法 CSV 导出
│   │ ── Python ──
│   ├── __init__.py
│   ├── matlab_bridge.py                 # MATLAB Engine 启动/调用/关闭
│   ├── config_io.py                     # JSON 配置加载/合并/保存
│   └── result_io.py                     # summary JSON 解析/打印
│
├── paper_reproduction/                  # 论文 Section IV 复现
│   ├── paper_reproduction.m             # 复现入口
│   ├── run_paper_baseline_case.m        # 论文复现专用仿真 (含 iter_num 跟踪与解有限性检查)
│   └── results/                         # 输出目录 (运行时生成)
│       ├── paper_reproduction.mat
│       ├── csv/
│       └── per_algorithm/
│
├── verification/                        # HPIPM 约束验证
│   ├── verify_constraints_hpipm.m       # HPIPM 解约束验证 (原始非线性约束重检)
│   └── results/                         # 输出目录 (运行时生成)
│       └── hpipm_constraint_verification.mat
│
└── third_party/                         # 第三方求解器源码
    ├── blasfeo/                         # BLASFEO 线性代数库
    └── hpipm/                           # HPIPM QP 求解器
```

> 本地运行仿真后还会生成 `results/batch/`（批量仿真输出）、`scenario_bank/`（预生成场景）两个目录。
>
> `results/batch/` 目录结构：
> ```
> results/batch/
> ├── comparison_final.mat     # 完整 comparison 结构体
> ├── comparison_checkpoint.mat # 断点续跑 checkpoint
> ├── csv/                     # {算法名}_results.csv
> ├── per_algorithm/           # {算法名}_summary.png + {算法名}/seed_XXXX.png
> └── comparison/              # fig3/4/5 论文对比图
> ```

## 算法引用架构（git submodule）

四种算法**统一放在 `algorithms/` 目录**下。其中 `RSS_proposed` 已转为普通目录（HPIPM 版 `control_RSS.m` 直接内嵌于主仓库），其余三个仍以 git submodule 引用外部仓库：

| 算法 | 来源 | 接口 |
|---|---|---|
| proposed-3iter | `algorithms/RSS_proposed/` (普通目录, 已内嵌 HPIPM 版 control_RSS.m) | `[u, new_state_dot, velocity, diagnostics] = control_RSS(path, k, state_dot, state)` |
| e-lmpc | `algorithms/RSS_sqp/` (submodule) | `[new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state)` |
| interior-point | `algorithms/RSS_fmincon/` (submodule) | `[new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state, params)` |
| active-set | `algorithms/RSS_active_set/` (submodule, active-set 分支) | `[new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state, params)` |

> **proposed-3iter 求解器**：`algorithms/RSS_proposed/control_RSS.m`（已转为普通目录，随主仓库直接 checkout，无需 submodule update）通过 MATLAB `py.*` 接口调用同目录下的 [hpipm_qp_solver.py](algorithms/RSS_proposed/hpipm_qp_solver.py)，将原 CVX 建模的 RSS 子问题转为 dense QCQP 由 HPIPM 求解。SQP 外层循环（3 次迭代）保留在 MATLAB 端，每次迭代调用一次 Python 求解器。变量 `x = [u(:); nu(:)]`（36 维），含 18 条等式约束（动力学）、24 条凸 SOC 约束（轮速，转为二次形式）、48 条非凸线性化约束（grip，两组 R 矩阵）。

`run_one_case.m` 按算法名动态切换 submodule 路径（`addpath`/`rmpath`），避免四个 `control_RSS.m` 和 `config.m` 同名冲突。`run_paper_baseline_case.m`（论文复现专用）改为在循环外一次性 `addpath` 当前算法的 submodule 并用 `onCleanup` 注册 `rmpath`，循环内不再切换路径，从根本上避免多版本 `control_RSS` 函数缓存冲突。

### 获取项目（含 submodule）

```bash
# 方式 1: clone 时带 --recursive (自动初始化 RSS_sqp/RSS_fmincon/RSS_active_set 三个 submodule)
git clone --recursive <repo-url>

# 方式 2: 已 clone, 后续初始化 submodule
git submodule update --init --recursive
```

> **注**：`algorithms/RSS_proposed/` 已转为普通目录，其文件（含 HPIPM 版 `control_RSS.m` 和 `hpipm_qp_solver.py`）随主仓库直接 checkout，**无需 submodule update**。只有 `RSS_sqp`、`RSS_fmincon`、`RSS_active_set` 三个仍为 submodule，需要上述命令初始化。

## 命令速查

四种入口，按使用场景选择：

| 场景 | 入口 | 命令 |
|------|------|------|
| **4 算法批量对比（Python 入口）** | 根目录 `main.py` | `python main.py` |
| **4 算法批量对比（MATLAB）** | `batch_simulation/main.m` | `cd batch_simulation; main` |
| **论文 Section IV 复现** | `paper_reproduction/paper_reproduction.m` | `cd paper_reproduction; paper_reproduction` |
| **约束验证** | `verification/verify_constraints_hpipm.m` | `cd verification; verify_constraints_hpipm` |

> **main.py 与 batch_simulation/main.m 的关系**：`main.py` 通过 MATLAB Engine 调用 `batch_simulation/main.m`，两者行为完全一致（Step 0-7 全流程）。`main.py` 只是 Python 薄包装层，方便从命令行批量运行。参数对应：`--seeds` → `'seeds'`，`--algorithms` → `'algorithms'`，`--force-regen` → `'forceRegen'`。

## 快速开始

### 方式一：Python 入口

通过 `main.py` 调用 MATLAB Engine，执行与 `batch_simulation/main.m` 完全相同的批量仿真流程（Step 0-7）。Python 负责参数解析，MATLAB 负责闭环仿真；proposed-3iter 算法在 MATLAB 中通过 `py.*` 回调同目录下的 HPIPM 求解器。

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
# 默认: 1:10 seeds, 全部 4 种算法
python main.py

# 自定义 seed 范围
python main.py --seeds 1:100

# 只跑一种算法
python main.py --algorithms proposed-3iter

# 组合参数 + 强制重新生成场景
python main.py --seeds 1:50 --algorithms proposed-3iter,e-lmpc --force-regen

# 查看所有参数
python main.py --help
```

#### CLI 参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--seeds` | 种子范围，格式 `start:end` 或单个整数 | `1:10` |
| `--algorithms` | 算法列表，逗号分隔 | `e-lmpc,active-set,interior-point,proposed-3iter` |
| `--force-regen` | 强制重新生成场景文件 | 关闭 |

#### 退出码

| 退出码 | 含义 |
|--------|------|
| 0 | 成功（全部 (seed, algorithm) 组合完成） |
| 2 | 仿真失败 |
| 3 | 基础设施错误（MATLAB Engine 未安装等） |

#### 输出

每次运行在 `results/batch/` 目录生成（与 `batch_simulation/main.m` 完全一致）：
- `comparison_checkpoint.mat` — 断点续跑 checkpoint（增量保存）
- `comparison_final.mat` — 完整 comparison 结构体
- `csv/{算法名}_results.csv` — 单算法指标 CSV
- `per_algorithm/{算法名}_summary.png` — 单算法汇总图
- `per_algorithm/{算法名}/seed_XXXX.png` — 逐 seed 轨迹图（由 `replot_per_seed` 生成）
- `comparison/fig3_trajectories.png` 等 — 论文对比图（Step 6）
- Table II 汇总打印（Step 7）

### 方式二：MATLAB 入口

批量对比 proposed-3iter / e-lmpc / interior-point / active-set 四种算法，支持多 seed。

#### 环境要求

- MATLAB（需 Optimization Toolbox 提供 fmincon）
- **Statistics and Machine Learning Toolbox**（可选，仅 Fig.4 求解时间箱线图用到 `boxplot`）。未安装时 [plot_paper_comparison.m](batch_simulation/plot_paper_comparison.m) 会自动回退到基础 MATLAB 的 `boxchart`（R2020a+ 自带），不影响仿真主流程
- git（用于 submodule 管理）
- Python 3.9+（HPIPM Python 接口，用于 proposed-3iter 算法）
- MSYS2 + gcc + make + bc（编译 HPIPM/BLASFEO 共享库，仅 Windows 需要）

#### HPIPM 求解器编译（仅 proposed-3iter 算法需要）

`proposed-3iter` 算法通过 MATLAB `py.*` 接口调用 Python 端 [hpipm_qp_solver.py](algorithms/RSS_proposed/hpipm_qp_solver.py)，使用 HPIPM dense QCQP 求解器。需要先编译 `libhpipm.dll`。

**Windows (MSYS2)**：

1. 安装 MSYS2：https://www.msys2.org/
2. 在 MSYS2 UCRT64 终端安装工具链：
   ```bash
   pacboy sync:mman-git ucrt64/toolchain msys/make msys/bc
   ```
3. 运行编译脚本（PowerShell 中）：
   ```powershell
   C:\msys64\usr\bin\env.exe MSYSTEM=UCRT64 /usr/bin/bash -lc "/d/PROJECT/RSS_V2/algorithms/RSS_proposed/build_hpipm_windows.sh"
   ```
   编译产物：`third_party/hpipm/lib/libhpipm.dll`（+ `libhpipm.so`）、`third_party/blasfeo/lib/libblasfeo.a`
4. 验证 Python 加载：
   ```powershell
   python -c "import sys; sys.path.insert(0,'algorithms/RSS_proposed'); import hpipm_qp_solver; print('HPIPM_OK=', hpipm_qp_solver._HPIPM_OK)"
   ```
   输出 `HPIPM_OK= True` 即成功。

> **DLL 路径自动设置**：[hpipm_qp_solver.py](algorithms/RSS_proposed/hpipm_qp_solver.py) 在 import 时自动通过 `os.add_dll_directory` 把 `third_party/hpipm/lib/` 加入 DLL 搜索路径，无需手动配置 PATH。
>
> **MATLAB Python 环境**：需在 MATLAB 中配置 Python：`pyenv('Version', '<python.exe 路径>')`。MATLAB Engine 调用 Python 时会继承上述 DLL 路径设置。
>
> **注（CVX → HPIPM 切换历史）**：原 RSS_proposed submodule 0121 分支 `control_RSS.m` 使用 CVX+SDPT3/ECOS 求解。ECOS 2.0.10 在 MATLAB R2026a 上触发 segfault，曾改用 SDPT3。现已完全替换为 HPIPM dense QCQP 求解器，不再依赖 CVX/ECOS/SDPT3。HPIPM 版本与原始 CVX 0121 分支对齐的关键点：`k1=1` 硬编码（非 `params.k1=0.15`）、SQP 外层循环无条件更新 `u_hat=u`（无失败回退/break）、QP 构造由 [construct_complete_qp_from_rss.m](algorithms/RSS_proposed/construct_complete_qp_from_rss.m) 完成（修复了 5 个数值 Bug）、Python 端仅负责求解（`balance` 模式 + tol=1e-6 + robust 回退）。

#### 运行批量对比

```matlab
cd('d:\PROJECT\RSS_V2\batch_simulation');

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

`batch_simulation/main.m` 按下列步骤顺序编排：

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

- **proposed 算法**（RSS_proposed）通过 MATLAB `py.*` 调用同目录下的 Python HPIPM 求解器，每次 SQP 迭代有跨语言调用开销（约 1-5ms/次 QCQP 求解），但求解时间由 `diagnostics` 结构体准确记录
- 三个 submodule 的 `control_RSS.m` 接口各不相同，`run_one_case.m` 负责适配
- submodule 内部调用各自的 `config.m`（无参函数），与项目的 `defaultConfig.m` 独立
- **main 批量仿真时**，`run_one_case.m` 通过 `setup_config_override` 用临时目录的 `config.m` 覆盖 submodule 的 `config.m`，使 proposed/e-lmpc 能用到 main 传入的随机 seed 场景参数（ctrl_pts/vimax/phidotmax/Lx/Ly 等）；interior-point/active-set 直接接收外部 config

## 论文结果复现

运行 [paper_reproduction.m](paper_reproduction/paper_reproduction.m) 可复现论文 Section IV 的固定轨迹实验。无需 seed，每个算法使用各自 submodule 的 config.m 参数：

```matlab
cd('d:\PROJECT\RSS_V2\paper_reproduction');

paper_reproduction                                  % 默认: 重跑全部 4 种算法
paper_reproduction({'proposed-3iter'})              % 只重跑 proposed-3iter, 保留其余算法已有结果
paper_reproduction({'e-lmpc','active-set'})         % 只重跑指定算法, 保留其余
```

**增量复现模式**：默认不清空整个 `paper_reproduction/results/` 目录，而是：
- 加载已有的 `paper_reproduction.mat`，保留未重跑算法的旧结果
- 仅清理本次将重跑算法对应的 CSV / 汇总图 / 逐 seed 图目录
- 最终合并保存（Table II 汇总覆盖全部 4 种算法）

若需全量重跑，先删除 `paper_reproduction/results/paper_reproduction.mat` 或传入全部算法列表。

**参数来源**（每个算法用各自的 config.m，不在 paper_reproduction 中覆盖）：

| 算法 | 参数来源 |
|---|---|
| proposed-3iter | `algorithms/RSS_proposed/config.m` |
| e-lmpc | `algorithms/RSS_sqp/config.m` |
| interior-point | `algorithms/RSS_fmincon/config.m` |
| active-set | `algorithms/RSS_active_set/config.m` |

**固定初始状态**（论文无随机性，所有算法统一）：

| 参数 | 值 | 适用算法 |
|---|---|---|
| 初始位姿 | `[0.05, 0.1, 0.2]` | 全部 4 种算法 |
| 初始速度 | `[0.01, 0.01, 0.01]` | 全部 4 种算法 |

> 注：RSS_proposed 已切换到 `0121` 分支（commit d25d299 "proposed method 最终版"），其 `main.m` 初始条件为 `[0.05,0.1,0.2]+[0.01,0.01,0.01]`，与其余三种算法完全一致，无需特判。

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

**与 batch_simulation/main.m 的区别**：
- `batch_simulation/main.m` 用 `scenario_bank(seed)` 生成随机场景，所有算法共用同一套参数（defaultConfig）
- `paper_reproduction.m` 每个算法用各自的 config.m 参数，仅初始状态固定

**输出**（写入 `paper_reproduction/results/`）：
- `per_algorithm/{算法名}/{算法名}_summary.png` — 单算法汇总图
- `per_algorithm/{算法名}/seed_0001.png` — 逐 seed 轨迹图
- `csv/{算法名}_results.csv` — 单算法指标 CSV
- `paper_reproduction.mat` — 完整 comparison 结构体

**额外诊断**（相比 `run_one_case.m`，由 `run_paper_baseline_case.m` 提供）：
- 每步 `iter_num`（求解器迭代数，proposed 为 NaN）
- 解的有限性检查（NaN/Inf 标记失败步）
- warm-up 排除后的中位数/分位数耗时（P1-4 复现要求）

## 约束验证（verify_constraints_hpipm）

[verify_constraints_hpipm.m](verification/verify_constraints_hpipm.m) 用于验证 HPIPM 求解的 RSS proposed 解是否满足**原始问题约束**（不仅仅是 QCQP 线性化形式）。

### 背景

rss_hpipm 分支的 RSS proposed 算法采用 MATLAB 构造 QCQP + Python HPIPM 求解的架构。HPIPM 严格满足 QCQP 形式的约束，但其中转向角速率约束是原始非凸约束的线性化近似，因此需要用原始非线性公式重新检验。

| 约束 | 类型 | QCQP 中的处理 | 是否可能违反 |
|------|------|-------------|-------------|
| 轮速 `‖H·ν‖ ≤ vimax` | **凸**（SOC） | HPIPM 直接强制 | 求解成功 → 不可能违反 |
| 转向角速率 `|Δθ| ≤ dt·phidotmax` | **非凸** | SQP 线性化近似（两组转向锥凸化约束） | **可能违反**（线性化仅局部近似） |

### 用法

```matlab
cd('D:\PROJECT\RSS_V2\verification')
verify_constraints_hpipm
```

### 验证逻辑

1. 调用 `control_RSS(path, step, last_vel, state)` 获取 HPIPM 实际输出的 `u` (3×K) 和 `velocity`
2. 由动力学等式约束重构完整 `nu` 序列：`nu(:,1)=current_nu+u(:,1)`，`nu(:,k+1)=nu(:,k)+u(:,k+1)`
3. 用**原始非线性公式**重检（100 步）：
   - 轮速：`norm(H{n} * nu(:,k))` vs `vimax`
   - 转角变化：`atan2(vy_k, vx_k) - atan2(vy_{k-1}, vx_{k-1})` vs `delta_theta = dt * phidotmax`
4. 输出汇总报告 + 保存 `hpipm_constraint_verification.mat`

### 输出

- 命令窗口打印每次违反的详细信息（步号、轮编号、超出量）
- 汇总报告：HPIPM 求解失败率、轮速违反步数/最大超出量、转角违反步数/最大超出量
- `hpipm_constraint_verification.mat`：完整验证数据（保存到 `verification/results/` 目录）

## 仿真步数说明

论文复现（`paper_reproduction.m`）固定 `t_end=1s, dt=0.01s` → `num_steps=100`，且 `num_path_pts = num_steps`（路径离散点数与仿真步数一一对应，第 k 步跟踪路径第 k 个点）。四种算法面对完全相同的 100 个参考点和 100 步控制，确保结果可比性。

随机场景生成器 `scenario_generator.m` 的 `t_end ∈ [0.5, 2.0]` 随机，步数随之变为 50~200，用于压力测试算法在不同路径长度下的鲁棒性。

> 注：`num_path_pts = num_steps` 将路径离散化精度与仿真时长绑定。在论文复现场景下无影响，但若未来需独立调整路径精度或仿真时长，应解耦这两个参数。

## 算法说明

| 算法 | 求解器 | 来源 | 特点 |
|---|---|---|---|
| proposed-3iter | HPIPM (dense QCQP, Python 接口) | RSS_proposed (普通目录) | RSS 凸化 + 3 次迭代，K=6，MATLAB 通过 `py.*` 调用 Python |
| e-lmpc | fmincon SQP | RSS_sqp submodule | MaxIter=1，K=6 |
| interior-point | fmincon interior-point | RSS_fmincon submodule | K=6，可违反转向锥约束 |
| active-set | fmincon active-set | RSS_active_set submodule | K=6，对齐 RSS_fmincon |

## 模块职责（速查）

| 模块 | 职责 |
|---|---|
| `main.py` | **Python 入口（推荐）**，CLI 参数解析 → MATLAB Engine → batch_simulation/main.m |
| `batch_simulation/main.m` | 批量仿真主入口，按 Step 0-7 编排全流程 |
| `batch_simulation/matlab_bridge.py` | MATLAB Engine 启动/路径设置/调用 main.m/关闭 |
| `batch_simulation/config_io.py` | JSON 配置加载/合并/保存 |
| `batch_simulation/result_io.py` | summary JSON 解析/打印 |
| `batch_simulation/run_one_case.m` | 批量仿真单场景闭环，按算法名分发到 submodule |
| `batch_simulation/comparison_init/load/save.m` | comparison 结构体的初始化 / 加载 / 增量保存 |
| `batch_simulation/run_batch_simulation.m` | 核心 (seed × alg) 双层循环 |
| `core/setup_paths.m` | 路径设置 + submodule 初始化检查 |
| `core/defaultConfig.m` | 默认参数（论文固定参数） |
| `core/generateReference.m` | Bernstein 多项式参考轨迹 |
| `core/propagateState.m` | 状态传播 |
| `core/computeWheelOutputs.m` | 轮速 / 轮角计算 |
| `core/computeMetrics.m` | RMSE / J / 求解时间 / 约束违反率 |
| `algorithms/RSS_proposed/control_RSS.m` | proposed-3iter 控制器（HPIPM dense QCQP, 3 次 SQP 迭代） |
| `algorithms/RSS_proposed/construct_complete_qp_from_rss.m` | QP 矩阵构造（修复 Bug1-5） |
| `algorithms/RSS_proposed/hpipm_qp_solver.py` | HPIPM dense QCQP 求解器（Python 接口） |
| `algorithms/RSS_proposed/build_hpipm_windows.sh` | Windows MSYS2 编译脚本 |
| `paper_reproduction/paper_reproduction.m` | 论文 Section IV 复现入口 |
| `paper_reproduction/run_paper_baseline_case.m` | 论文复现专用仿真（含 iter_num 跟踪与解有限性检查） |
| `verification/verify_constraints_hpipm.m` | HPIPM 解约束验证（原始非线性约束重检） |

## Python → MATLAB 架构

```
用户/CI
  ↓
根目录 main.py
  ├── 解析 CLI 参数 (seeds, algorithms, forceRegen)
  ├── 启动一次 MATLAB Engine
  ├── addpath(仓库根目录 + core/ + batch_simulation/ + algorithms/)
  ├── cd(batch_simulation/)
  ├── 调用 main('seeds', ..., 'algorithms', ..., 'forceRegen', ...)  ← 一次跨语言调用
  └── 根据结果设置退出码
            ↓
MATLAB batch_simulation/main.m
  ├── Step 0: setup_paths()
  ├── Step 3: comparison_load() (断点续跑)
  ├── Step 4: run_batch_simulation()
  │     └─ run_one_case()  按算法名分发到 algorithms/ 下
  │           └─ 每步调用 control_RSS (HPIPM dense QCQP, 3 次 SQP 迭代)
  ├── Step 5: save_all_artifacts()
  ├── Step 6: plot_paper_comparison()
  └── Step 7: print_table_ii()
```

**关键原则**：Python 只调用一次 MATLAB，不逐 MPC step 调用。整个 100 步闭环在 MATLAB 内部完成，保证计时公平性和 CVX warm start 的完整性。

## 关键约束

- submodule 以引用方式接入，不修改 submodule 内的代码（RSS_proposed 的 control_RSS.m 例外：已修复失败回退和计时解析）
- 算法实现细节与原仓库完全一致（CVX 问题公式、代价函数、约束、K、rho 均未修改）
- active-set 为 RSS_fmincon 仓库的 active-set 分支（同仓库不同分支，独立 submodule）
- checkpoint 文件不可删除，用于断点续跑
- Python 入口不计算机器人状态，不读写 MATLAB Base Workspace，不解析 MATLAB 控制台文本
- MATLAB Engine 只启动一次，整个闭环只跨语言调用一次
