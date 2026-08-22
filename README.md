# RSS-MPC-Comparison-Pipeline

基于论文 *RSS2026: Exploit Agile Mobility of Steerable-Wheeled Mobile Robots: A Fast Motion Planning Approach* 的四轮全向底盘 MPC 仿真复现项目，对比四种控制器。

> 项目主干为 `pipeline/` 统一单场景仿真入口（dt (tau) 固定 0.01 s，4 算法可选）；批量仿真 / 论文复现 / 约束验证等非主干内容统一放在 `others/`。

## 目录结构

```
RSS-MPC-Comparison-Pipeline-rss_hpipm/
├── pipeline/                            # 统一单场景 pipeline (4 算法可选, 项目主干)
│   ├── main.py                          #   入口: --seed/--algorithm/--K/--iters
│   ├── params.py                        #   参数定义 (dt=0.01 常量/车辆/权重/求解器)
│   ├── trajectory_generator.py          #   轨迹生成 (seed_id, K -> Bezier 参考 + 场景)
│   ├── dynamics.py                      #   原始动力学 (状态传播 + 轮子正运动学)
│   ├── simulator.py                     #   闭环仿真器 (4 算法统一闭环)
│   ├── metrics.py                       #   评估 (RMSE/J 分解/求解时长/约束违反)
│   ├── plotting.py                      #   结果画图 (轨迹/误差/轮速/求解时长...)
│   ├── matlab_algorithm.py              #   MATLAB Engine 算法桥 (e-lmpc 等 3 算法每步求解)
│   └── matlab_control_bridge.m          #   MATLAB 侧每步求解桥 (Engine 调 control_RSS)
│
├── algorithms/                          # 算法包
│   ├── RSS_proposed/                    # proposed 生产实现 (纯 Python) + golden oracle
│   │   ├── control_rss_ocpqcqp.py       #   proposed 控制器 (SCP 外层, 迭代数可指定)
│   │   ├── construct_ocp_qcqp.py        #   OCP QCQP 矩阵构造
│   │   ├── hpipm_qp_solver.py           #   HPIPM Python 接口 (solve_ocp_qcqp / solve_qcqp)
│   │   ├── control_RSS_denseqcqp.m      #   Dense QCQP golden oracle (离线对照)
│   │   ├── construct_complete_qp_from_rss.m #   Dense QCQP 矩阵构造 (golden oracle)
│   │   ├── config.m                     #   算法参数
│   │   └── build_hpipm_windows.sh       #   Windows MSYS2 编译脚本
│   ├── RSS_sqp/                         # e-lmpc 算法 (fmincon SQP, 普通目录, 原 submodule 已并入)
│   ├── RSS_fmincon/                     # interior-point 算法 (fmincon IPM, 普通目录, 原 submodule 已并入)
│   ├── RSS_active_set/                  # active-set 算法 (fmincon SQP, 普通目录, 原 submodule 已并入)
│   └── .gitattributes
│
├── others/                              # 与 pipeline 主干无关的批量/复现/校验内容
│   ├── setup_paths.m                    # 路径设置 (含算法目录检查)
│   ├── batch_simulation/                # MATLAB 批量仿真 (多 seed 随机场景)
│   │   ├── main.m / main.py             # 批量入口 (MATLAB / Python)
│   │   ├── run_one_case.m               # 单场景闭环仿真, 按算法名分发到算法目录
│   │   ├── run_batch_simulation.m       # (seed × algorithm) 批量仿真循环
│   │   ├── comparison_init.m / _load.m / _save.m  # comparison 结构体管理
│   │   ├── plot_one_algorithm.m         # 单算法 summary 图
│   │   ├── plot_paper_comparison.m      # 论文风格对比图 (Fig.3/4/5)
│   │   ├── replot_per_seed.m            # 逐 seed 轨迹重绘
│   │   ├── print_table_ii.m             # Table II 汇总打印
│   │   ├── save_algorithm_csv.m         # 单算法 CSV 导出
│   │   └── matlab_bridge.py / config_io.py / result_io.py  # Python 桥接
│   ├── paper_reproduction/              # 论文 Section IV 复现 (MATLAB)
│   │   ├── paper_reproduction.m         # 复现入口
│   │   └── run_paper_baseline_case.m    #   论文复现专用仿真 (proposed 经 Dense QCQP golden oracle)
│   └── verification/                    # HPIPM 约束验证
│       └── verify_constraints_hpipm.m   #   HPIPM 解约束验证脚本
│
├── scenario_bank/                       # 场景库 (scenario_seed{1..150}.mat, 运行时生成)
│
├── third_party/                         # 第三方求解器源码
│   ├── blasfeo/                         # BLASFEO 线性代数库 (submodule)
│   └── hpipm/                           # HPIPM QP/QCQP 求解器
│
├── README.md                            # 本文件
└── RSS_proposed_算法详解.md             # proposed 算法推导与实现说明
```

## 算法说明

| 算法 | 求解器 | 来源 | 特点 |
|---|---|---|---|
| proposed-3iter | HPIPM (OCP QCQP + SCP) | RSS_proposed/ (纯 Python) | 原生凸二次约束 + SCP 迭代，默认 K=6，支持 `--K` 指定；Dense QCQP golden oracle 同目录 |
| e-lmpc | fmincon SQP | RSS_sqp/ (普通目录) | MaxIter=1，预测时域 K 可通过 `--K` 或 `config.K` 配置 |
| interior-point | fmincon interior-point | RSS_fmincon/ (普通目录) | 预测时域 K 可通过 `--K` 或 `config.K` 配置 |
| active-set | fmincon active-set | RSS_active_set/ (普通目录) | 预测时域 K 可通过 `--K` 或 `config.K` 配置 |

> MATLAB 三算法逐步日志统一口径：`exitflag` 为 fmincon 原生码（正值 1~5 均为收敛，仅判据不同；0=达迭代上限；负值=失败），另打印 `status = sign(exitflag)`（1=收敛 / 0=达 MaxIterations / -1=失败）。e-lmpc 为 "1 iteration edition" 设计（每步单次 SQP 迭代），其 status=0 属预期而非失败。

### proposed-3iter 求解策略 (OCP QCQP + SCP)

论文 Algorithm 1 的凸子问题 (公式 17) 含二次约束 (转向锥、轮速 SOC)。本分支采用 **HPIPM 原生 OCP QCQP** 求解，无需将二次约束线性化为切平面：

| 组件 | 说明 |
|---|---|
| `solve_ocp_qcqp` | HPIPM 的 `ocp_qcqp` IPM 求解器，原生处理凸二次约束 (转向锥 + 轮速 SOC) |
| `algorithms/RSS_proposed/construct_ocp_qcqp.py` | 构造 OCP QCQP 矩阵 (A/B/Bb/Q/R/S/q/r/Qq/Sq/Rq/qq/rq/uq)，二次约束直接以二次形式给出 |
| **SCP 外层** | 固定锚点 `u_hat` 后原非凸问题转化为凸 QCQP；每个 MPC step 执行 3 次 SCP 迭代 (论文基准; 纯 Python pipeline 可用 `--iters` 指定) |
| **warm_start=2** | 选中 HPIPM "全量裁剪初始化"分支，绕过 DLL C 层 `d_ocp_qcqp_ipm_arg_set_t0_init` 误写 `t_lam_min` 的 bug (step 94 剧烈机动段触发) |
| **求解器参数** | `mode=speed`、`mu0=10`、`warm_start=2` (通过环境变量 `HPIPM_OCP_QCQP_MODE` 切换) |

与 Dense QCQP golden oracle（离线对照基准，运行时禁止作为 fallback）的 100 步闭环对齐验证：

| 指标 | Dense QCQP (oracle) | OCP QCQP (默认) | 差异 |
|---|---|---|---|
| RMSE | 0.036793 | 0.036793 | < 1e-6 |
| J_total | 13.3838 | 13.3838 | < 1e-6 |
| validSteps | 100/100 | 100/100 | — |
| medianSolveTime | ~0.008s | 0.0016s | OCP QCQP 快约 5× |

> 结论：PASS (100 步 Golden 对齐)。OCP QCQP 在数值精度上与 Dense QCQP 完全等价，且因利用 OCP 结构稀疏性而显著更快。

## pipeline 命令行调用

```bash
# proposed-3iter: 论文固定场景 (seed=0), K=6, 默认 3 次 SCP 迭代 (纯 Python)
python pipeline/main.py --seed 0 --algorithm proposed-3iter --K 6

# proposed-3iter: scenario_bank 随机场景 (seed>=1), 指定 5 次 SCP 外层迭代
python pipeline/main.py --seed 1 --algorithm proposed-3iter --K 6 --iters 5

# MATLAB 算法: 经 MATLAB Engine 每步求解 (引擎启动一次约 30s, 需 matlabengine 包)
python pipeline/main.py --seed 0 --algorithm e-lmpc --K 6
python pipeline/main.py --seed 1 --algorithm active-set --K 6
python pipeline/main.py --seed 0 --algorithm interior-point --K 6

# 覆盖约束值 (vimax 轮速上限 m/s / phidotmax 转向速率上限 rad/s; 默认用场景值)
python pipeline/main.py --seed 0 --algorithm proposed-3iter --K 6 --vimax 8 --phidotmax 10

# 控制正则化 rho: 单值 (常数) 或逗号分隔逐步序列 (第 k 步取第 k 个值,
# 序列短于总步数时保持末值; 如 1,0 = 第 1 步 rho=1, 之后 rho=0)
python pipeline/main.py --seed 0 --algorithm proposed-3iter --K 6 --rho 1,0
python pipeline/main.py --seed 0 --algorithm proposed-3iter --K 6 --rho 0.05

# 模块方式启动
python -m pipeline.main --seed 0 --algorithm proposed-3iter --K 6
```

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--seed` | 场景 seed_id (0=paper_fixed 固定场景, >=1=scenario_bank) | `0` |
| `--algorithm` | 算法名, 4 选 1: `proposed-3iter` / `e-lmpc` / `active-set` / `interior-point`。后三者经 MATLAB Engine 每步求解 | `proposed-3iter` |
| `--K` | MPC 预测时域步长 (所有算法均可指定，默认 `6`) | `6` |
| `--iters` | proposed 的 SCP 外层迭代数 (每步 HPIPM 求解次数)。默认 3 = 论文基准；非 3 时结果目录名带 `_iters{N}` 后缀 | `3` |
| `--vimax` | 覆盖最大轮速约束 vimax (m/s)，4 算法统一生效；不传用场景值 (seed=0 为论文基准 5)。结果目录名带 `_vimax{N}` 后缀 | 场景值 |
| `--phidotmax` | 覆盖最大转向角速率约束 phidotmax (rad/s)，4 算法统一生效；不传用场景值 (seed=0 为论文基准 5π)。结果目录名带 `_phidot{N}` 后缀 | 场景值 |
| `--rho` | 控制正则化 rho：单值 (如 `0.01`) 或逗号分隔逐步序列 (如 `1,0` = 第 1 步 rho=1 之后 rho=0)，序列短于总步数时保持末值。真正生效于 proposed-3iter 的 RSS 锚点正则化 (论文式 17)；3 个 fmincon 基线算法的 rho 为目标函数签名中的遗留参数 (未参与代价)，传入仅贯通管路、不改变其基线结果。结果目录名带 `_rho...` 后缀 | 各算法默认 |
| `--out` | 结果输出根目录 | `results/pipeline` |
| `--quiet` | 关闭算法 per-iteration 打印 | 关闭 |

运行后在 `results/pipeline/{algorithm}/seed{N}_K{K}/` 生成结果三件套：

1. **`run_config.json`** — 记录传入参数 (seed_id / 算法 / K / iters / dt / 场景 / 权重 / 求解器)
2. **`metrics.json`** — 评估指标 (RMSE / J_total / medianSolveTime / 约束违反量 / 成功率...)
3. **`figures/*.png`** — 画图 (轨迹 / 跟踪误差 / 轮速 / 转向速率 / 求解时长 / 控制输入)
4. **`simulation_data.npz`** — 原始数组 (复画图 / golden 对比用)

环境要求：proposed-3iter 需 HPIPM DLL（Windows 编译脚本 `algorithms/RSS_proposed/build_hpipm_windows.sh`，构建依赖 `third_party/blasfeo` submodule，clone 后执行 `git submodule update --init third_party/blasfeo`）；MATLAB 算法需 `matlabengine` 包（版本须与 MATLAB 发行版匹配，安装见 `pipeline/matlab_algorithm.py` 模块注释）+ MATLAB 在 PATH。三个算法目录（RSS_sqp/RSS_fmincon/RSS_active_set）已作为普通目录并入主仓库，clone 即得，无需 submodule 操作。

## 各板块调用逻辑与流程图

### 总体流程图

```
python pipeline/main.py --seed S --algorithm A --K K [--iters N]
        │  (dt = 0.01 s 写死)
        ▼
trajectory_generator.py + params.py
 (seed,K) → Bezier 参考 + 场景; 车辆/权重/求解器参数
        │
        ▼
simulator.py (原始动力学闭环, 按步长 K 循环调用算法 A)
        │
        │  ◄── 唯一分支点: 仿真器每步调用算法时按 A 分发
        │
        ├─ A = proposed-3iter (纯 Python 路径)
        │      ▼
        │  algorithms/RSS_proposed/control_rss_ocpqcqp.py (SCP 外层 × iters)
        │      │ 每次迭代
        │      ▼
        │  algorithms/RSS_proposed/construct_ocp_qcqp.py (OCP QCQP 矩阵)
        │      ▼
        │  algorithms/RSS_proposed/hpipm_qp_solver.py → third_party/hpipm DLL
        │      │ u* 返回 simulator (状态推进见下方 dynamics.py)
        │
        └─ A ∈ {e-lmpc, active-set, interior-point} (MATLAB 算法, 每步借 MATLAB 求解)
               ▼
           pipeline/matlab_algorithm.py —— Python 侧的"传话人"
           (仿真开始时把 MATLAB 启动一次, 约 30 秒, 之后一直开着不重复启动;
            参考轨迹、K、车辆和权重参数也在开头一次性交给 MATLAB 存好,
            每步只传"当前是第几步、现在的速度和位置"这几个数)
               ▼
           pipeline/matlab_control_bridge.m —— MATLAB 侧的"接应"
           (第一步被调用时, 把 Python 传来的参数整理成算法认识的 config
            结构体, 之后每步直接复用; 参数怎么交给算法, 取决于算法的
            函数签名: interior-point / active-set 的第 5 个参数就是
            config, 直接递进去; e-lmpc 的签名不收参数, 它内部自己调
            config() 从文件里读, 所以接应把整理好的参数写进一个临时
            config.m 放到路径最前面"顶替"原文件)
               ▼
           algorithms/{RSS_sqp | RSS_fmincon | RSS_active_set}/control_RSS.m
           (真正干活的算法: 用 fmincon 求解这一步的最优控制, 解完打印
            一行 exitflag+status 方便排查收敛情况)
               ▼
           解出的控制量/速度传回 Python
        │
        ▼  (两条分支在此汇合: 算法解出的控制量都交回仿真器)
dynamics.py —— 机器人的"真身", 按原始动力学推进一步
 (propagate_state: 用世界速度把位姿 [x;y;ψ] 前推一个 dt, 不做任何
  线性化——算法内部的近似模型只管算控制量, 机器人按真模型动,
  两者的差距才是真实跟踪误差; compute_wheel_outputs 顺带算出
  各轮速度/转向角, 记录下来供画图)
        │
        │  步数没走完 且 没有失败 (step_failed/异常)?
        │      是 ──► 回到上面的"仿真器每步调用算法"处, 带着
        │             新位姿和新速度走第 k+1 步 (就这样闭环转满
        │             num_steps 圈; seed=0 是 100 步)
        ▼ 否
统一 SimResult (两条路径同构)
        │
   ┌────┼─────────┐
   ▼    ▼         ▼
metrics.py  plotting.py  结果三件套落盘
(RMSE/J/耗时) (6 张图) (run_config.json / metrics.json / figures/ / npz)
```

> 说明：仿真器对 4 种算法一视同仁（同一条"循环调用算法"主干）：轨迹生成、参数配置、原始动力学闭环推进、评估与画图全部在 Python 侧；MATLAB 系算法仅在每步求解时经 MATLAB Engine 调用算法目录的 `control_RSS.m`（引擎会话随仿真开启/关闭，启动一次约 30s，之后每步毫秒级调用）。

### 各板块调用逻辑

| 板块 | 文件 | 职责 | 上游调用者 | 下游被调对象 |
|---|---|---|---|---|
| 入口 | `pipeline/main.py` | 解析 CLI (seed/algorithm/K/iters)，4 算法统一走 Python 主干 (轨迹→仿真→评估→落盘)，统一落盘结果三件套 | 用户命令行 | `trajectory_generator` / `simulator` / `metrics` / `plotting` |
| 参数定义 | `pipeline/params.py` | `DT=0.01` 常量、车辆参数（轮位/vimax/phidotmax）、权重、求解器参数、RunConfig 记录结构 | 所有 pipeline 模块 | — |
| 轨迹生成 | `pipeline/trajectory_generator.py` | (seed_id, K) → Bernstein/Bezier 参考轨迹 + 场景参数；seed=0 为 paper_fixed 固定场景，seed≥1 读 `scenario_bank/scenario_seed{N}.mat`（与 MATLAB 同源） | `main.py` | `simulator.py` |
| 原始动力学 | `pipeline/dynamics.py` | `propagate_state`（状态传播）+ `compute_wheel_outputs`（轮速/轮角正运动学）；使用**原始动力学**而非算法内部展开/线性化误差动力学 | `simulator.py` | — |
| 闭环仿真器 | `pipeline/simulator.py` | 原始动力学闭环：每步 `control(path, step, state_dot, state, params)` → u/世界速度/车体速度/诊断，推进状态并逐步记录；step_failed/NaN 检查失败即终止；4 算法均在此闭环（MATLAB 算法经 `MatlabAlgorithmBridge` 每步求解） | `main.py` | `dynamics` / `control_rss_ocpqcqp` / `matlab_algorithm` |
| proposed 控制器 | `algorithms/RSS_proposed/control_rss_ocpqcqp.py` | RSS 控制律：SCP 外层循环 `--iters` 次（默认 3），每次构造并求解凸 QCQP 子问题，失败时保留最近可行 incumbent | `simulator.py` | `construct_ocp_qcqp` |
| QCQP 构造 | `algorithms/RSS_proposed/construct_ocp_qcqp.py` | RSS 模型 → OCP QCQP 矩阵（A/B/Bb/Q/R/S/q/r/Qq/Sq/Rq/qq/rq/uq），转向锥/轮速 SOC 以原生二次约束给出 | `control_rss_ocpqcqp.py` | `hpipm_qp_solver` |
| HPIPM 接口 | `algorithms/RSS_proposed/hpipm_qp_solver.py` | ctypes 封装，加载 `third_party/hpipm/lib/libhpipm.dll`，调 `ocp_qcqp` IPM 求解器 | `construct_ocp_qcqp.py` | HPIPM/BLASFEO DLL |
| MATLAB Engine 桥 | `pipeline/matlab_algorithm.py` | `MatlabAlgorithmBridge`：启动常驻 MATLAB Engine 会话（一次）；轨迹/算法/seed 经 base workspace 传入，完整 run_config（K/dt/车辆/权重）写临时 .mat 加载为 `PIPELINE_CONFIG`，K 另写 `PIPELINE_K`；每步 `control()` 调 `matlab_control_bridge` 求解并回传，签名与 `control_rss_ocpqcqp` 统一 | `simulator.py` | `matlab_control_bridge.m` |
| MATLAB 侧求解桥 | `pipeline/matlab_control_bridge.m` | persistent 一次性初始化（从 `PIPELINE_CONFIG` 重建算法 config 结构体：车辆几何→Lx/Ly/wheel_pos/vimax/phidotmax、权重→k1/k2/k3/eps、dt/num_steps/K，`PIPELINE_K` 覆盖预测时域，保证与 Python 侧参数一致）；e-lmpc/active-set 写临时 config.m 覆盖算法目录 config()，interior-point 直接传 config 第 5 参；每步按算法分发 `control_RSS`，evalc 捕获 exitflag/status 日志回传 | `matlab_algorithm.py` | 算法目录 `control_RSS.m` |
| 单 case 仿真 | `others/batch_simulation/run_one_case.m` | MATLAB 版闭环（批量仿真用）：按算法名 addpath 对应算法目录并临时覆盖 config，循环调用 `control_RSS.m`，原始动力学推进 | `others/batch_simulation/main.m` | 算法目录 `control_RSS.m` |
| 对比算法 | `algorithms/{RSS_sqp,RSS_fmincon,RSS_active_set}/control_RSS.m` | fmincon 系 NLP 求解（SQP / interior-point / active-set），统一 exitflag+status 日志口径 | `matlab_control_bridge.m` / `run_one_case.m` | fmincon |
| 评估 | `pipeline/metrics.py` | RMSE / J_total 与 J 分解 / medianSolveTime（排除 warm-up，论文 P1-4 口径）/ 约束违反量 / 成功率 | `main.py` | — |
| 画图 | `pipeline/plotting.py` | 6 张图：轨迹跟踪 / 跟踪误差 / 轮速（含约束线）/ 转向速率 / 每步求解时长 / 控制输入 | `main.py` | — |

### 与 MATLAB 的一致性

- proposed-3iter golden 基准 (seed=0, K=6, iters=3) 与 MATLAB 完全对齐：RMSE=0.036793, J_total=13.3838, validSteps=100/100
- MATLAB Engine 桥接算法 (Python 主干 + MATLAB 每步求解) 与 MATLAB 侧闭环结果完全一致：e-lmpc seed=0: RMSE=0.037776, J_total=44.4226, validSteps=100/100（两侧完全一致）
