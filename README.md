# RSS-MPC-Comparison-Pipeline

基于论文 *RSS2026: Exploit Agile Mobility of Steerable-Wheeled Mobile Robots: A Fast Motion Planning Approach* 的四轮全向底盘 MPC 仿真复现项目，对比四种控制器。

> 当前分支 `hpipm_ocp_qcqp`：proposed-3iter 的默认求解器已从 legacy OCP QP（切平面线性化 + SCP）升级为 **HPIPM 原生 OCP QCQP**（凸二次约束 + SCP），与 Dense QCQP golden oracle 在 100 步闭环测试中数值一致。

## 目录结构

```
RSS-MPC-Comparison-Pipeline-rss_hpipm/
├── main.py                              # Python 入口, 通过 MATLAB Engine 调用
├── default_python_config.json           # Python 默认配置
│
├── algorithms/                          # 算法包
│   ├── RSS_proposed/                    # proposed 控制器 (HPIPM OCP QCQP + SCP)
│   │   ├── control_RSS_ocpqcqp.m        #   默认: OCP QCQP 主控制器 (SCP 外层循环)
│   │   ├── control_RSS_denseqcqp.m      #   Dense QCQP golden oracle (离线对照)
│   │   ├── construct_ocp_qcqp_from_rss.m    #   OCP QCQP 矩阵构造 (原生二次约束)
│   │   ├── construct_complete_qp_from_rss.m #   Dense QCQP 矩阵构造 (golden oracle)
│   │   ├── check_qk_violation.m         #   约束违反诊断工具
│   │   ├── config.m                     #   算法参数
│   │   ├── hpipm_qp_solver.py           #   HPIPM Python 接口 (solve_ocp_qcqp / solve_qcqp)
│   │   ├── build_hpipm_windows.sh       #   Windows MSYS2 编译脚本
│   │   ├── tests/                       #   测试脚本 (Test A/B/C/D)
│   │   │   ├── test_ocp_qcqp_construction_equivalence.m  # Test A: 构造等价性
│   │   │   ├── test_hpipm_ocp_qcqp_smoke.m              # Test B: wrapper 冒烟测试
│   │   │   ├── test_ocp_qcqp_fixed_qk_alignment.m       # Test C: 固定 QK 对齐
│   │   │   └── test_ocp_qcqp_golden_100steps.m          # Test D: 100 步闭环 Golden 对比
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
│   ├── run_paper_baseline_case.m        # 论文复现专用仿真 (默认 solver=ocpqcqp)
│   └── results/                         # 输出目录 (运行时生成)
│
├── verification/                        # HPIPM 约束验证
│   ├── verify_constraints_hpipm.m       # HPIPM 解约束验证
│   └── results/                         # 输出目录 (运行时生成)
│
└── third_party/                         # 第三方求解器源码
    ├── blasfeo/                         # BLASFEO 线性代数库 (submodule)
    └── hpipm/                           # HPIPM QP/QCQP 求解器
```

## 算法说明

| 算法 | 求解器 | 来源 | 特点 |
|---|---|---|---|
| proposed-3iter | HPIPM (OCP QCQP + SCP) | RSS_proposed (本地目录) | 原生凸二次约束 + SCP 迭代, K=6, 默认 solver |
| e-lmpc | fmincon SQP | RSS_sqp (submodule) | MaxIter=1, K=6 |
| interior-point | fmincon interior-point | RSS_fmincon (submodule) | K=6 |
| active-set | fmincon active-set | RSS_active_set (submodule) | K=6 |

### proposed-3iter 求解策略 (OCP QCQP + SCP)

论文 Algorithm 1 的凸子问题 (公式 17) 含二次约束 (转向锥、轮速 SOC)。本分支采用 **HPIPM 原生 OCP QCQP** 求解，无需将二次约束线性化为切平面：

| 组件 | 说明 |
|---|---|
| `solve_ocp_qcqp` | HPIPM 的 `ocp_qcqp` IPM 求解器，原生处理凸二次约束 (转向锥 + 轮速 SOC) |
| `construct_ocp_qcqp_from_rss` | 构造 OCP QCQP 矩阵 (A/B/Bb/Q/R/S/q/r/Qq/Sq/Rq/qq/rq/uq)，二次约束直接以二次形式给出 |
| **SCP 外层** | 固定锚点 `u_hat` 后原非凸问题转化为凸 QCQP；每个 MPC step 执行 3 次 SCP 迭代 |
| **warm_start=2** | 选中 HPIPM "全量裁剪初始化"分支，绕过 DLL C 层 `d_ocp_qcqp_ipm_arg_set_t0_init` 误写 `t_lam_min` 的 bug (step 94 剧烈机动段触发) |
| **求解器参数** | `mode=speed`、`mu0=10`、`warm_start=2` (通过环境变量 `HPIPM_OCP_QCQP_MODE` 切换) |

#### 与 Dense QCQP golden oracle 的对齐验证

100 步闭环 Golden 对比 (`test_ocp_qcqp_golden_100steps.m`)：

| 指标 | Dense QCQP (oracle) | OCP QCQP (默认) | 差异 |
|---|---|---|---|
| RMSE | 0.036793 | 0.036793 | < 1e-6 |
| J_total | 13.3838 | 13.3838 | < 1e-6 |
| validSteps | 100/100 | 100/100 | — |
| medianSolveTime | ~0.008s | 0.0016s | OCP QCQP 快约 5× |

> 结论：PASS (100 步 Golden 对齐)。OCP QCQP 在数值精度上与 Dense QCQP 完全等价，且因利用 OCP 结构稀疏性而显著更快。

### Dense QCQP golden oracle

`control_RSS_denseqcqp.m` + `construct_complete_qp_from_rss.m` 保留作为**离线对照基准**，用于验证 OCP QCQP 实现的正确性。运行时禁止使用 Dense QCQP 作为 fallback（仅作为测试基准）。

### Legacy OCP QP (已删除)

早期实现 `control_RSS.m` + `construct_ocp_qp_from_rss.m` 通过切平面线性化将二次约束转为线性约束，数学上不等价于原问题，已被 OCP QCQP 实现取代并删除。

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
| OCP QCQP Golden 对比 | `test_ocp_qcqp_golden_100steps.m` | 见下方"测试"章节 |

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

### 指定 solver (proposed-3iter 专用)

`run_paper_baseline_case.m` 默认使用 `solver='ocpqcqp'`。可通过 `config.solver` 切换：

```matlab
cfg = defaultConfig();
cfg.algorithm = 'proposed-3iter';
cfg.solver = 'ocpqcqp';   % 默认: 原生 OCP QCQP (推荐)
% cfg.solver = 'denseqcqp'; % 可选: Dense QCQP golden oracle (离线对照)
run_paper_baseline_case(cfg);
```

## 测试 (OCP QCQP 验证)

位于 `algorithms/RSS_proposed/`，按以下顺序执行：

| 测试 | 文件 | 验证内容 |
|---|---|---|
| Test A | `test_ocp_qcqp_construction_equivalence.m` | OCP QCQP 与 Dense QCQP 矩阵构造等价性 |
| Test B | `test_hpipm_ocp_qcqp_smoke.m` | HPIPM OCP QCQP wrapper 冒烟测试 (从最小问题到 72 条约束) |
| Test C | `test_ocp_qcqp_fixed_qk_alignment.m` | 固定 QK 时 OCP QCQP 与 Dense QCQP 目标函数一致性 |
| Test D | `test_ocp_qcqp_golden_100steps.m` | 100 步闭环 Golden 对比 (RMSE/J_total/validSteps) |

运行 Test D：

```matlab
cd('d:/PROJECT/RSS-MPC-Comparison-Pipeline-rss_hpipm');
addpath('core'); setup_paths;
addpath('algorithms/RSS_proposed');
addpath('algorithms/RSS_proposed/tests');
addpath('paper_reproduction');
setenv('HPIPM_OCP_QCQP_MODE', 'speed');   % 必须设为 speed 以绕过 DLL bug
test_ocp_qcqp_golden_100steps;
```

## 约束验证

```matlab
cd('verification');
verify_constraints_hpipm
```

验证 HPIPM 求解的解是否满足原始非线性约束（轮速 SOC 约束、转向角速率约束）。输出汇总报告 + `verification/results/hpigm_constraint_verification.mat`。

---

## RSS_proposed 算法文件详解

详见 [RSS_proposed_算法详解.md](RSS_proposed_算法详解.md)（含文件总览、物理量参数、`construct_complete_qp_from_rss.m` 的 Dense QCQP 矩阵构造、`construct_ocp_qcqp_from_rss.m` 的 OCP QCQP 矩阵构造、调用逻辑、关键设计点）。
