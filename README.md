# RSS: Robotic Steering Safety Controller

基于论文 *RSS: A Reactive Safety Shield for Robotic Systems* 的四轮全向底盘模型预测控制（MPC）仿真复现项目。实现并对比四种控制器：proposed（HPIPM 凸 QP 迭代）、e-LMPC（SQP 单步）、active-set、interior-point。

## 目录结构

```
RSS/
├── matlab/                     # 仿真框架（主入口）
│   ├── setup_paths.m           # 路径设置
│   ├── defaultConfig.m         # 默认参数（论文固定参数）
│   ├── scenario_bank.m         # 场景库（可复现随机场景）
│   ├── scenario_generator.m    # Latin 超立方采样
│   ├── generateReference.m     # Bernstein 多项式参考轨迹
│   ├── run_one_case.m          # 单场景闭环仿真
│   ├── run_paper_baseline_case.m  # 基线专用仿真（含 exitflag 跟踪）
│   ├── propagateState.m        # 状态传播
│   ├── computeWheelOutputs.m   # 轮速/轮角计算
│   ├── computeMetrics.m        # 指标计算（RMSE/J/求解时间/约束违反率）
│   ├── compare_algorithms.m    # 批量对比主入口（含断点续跑）
│   ├── paper_reproduction.m    # 论文 Section IV 复现独立脚本
│   ├── plot_paper_comparison.m # 论文对比图
│   ├── replot_per_seed.m       # 逐 seed 轨迹重绘
│   └── clean_baseline_checkpoint.m  # 清理基线 checkpoint
├── RSS_proposed/               # 控制器与求解器
│   ├── control_RSS_v2.m        # proposed（HPIPM 3 次迭代）
│   ├── control_eLMPC.m         # e-LMPC（SQP MaxIter=1）
│   ├── control_active_set.m    # active-set
│   ├── control_interior_point.m # interior-point
│   ├── construct_qp_from_rss_control.m  # QP 构造
│   ├── construct_complete_qp_from_rss.m # 完整 QP 构造
│   ├── solve_qp_with_python_hpipm.m     # Python HPIPM 桥接
│   ├── hpipm_solver_wrapper.m            # 求解器封装
│   └── config_hpipm.m                    # HPIPM 配置
├── python/
│   └── hpipm_solver.py         # HPIPM Python 接口
└── tests/                      # Python 单元测试
    ├── test_matlab_baseline.py
    └── test_python_rss_controller.py
```

> 本地运行仿真后还会生成 `results/`（输出结果）、`scenario_bank/`（预生成场景）、`third_party/`（HPIPM/BLASFEO 源码）三个目录，已在 `.gitignore` 中排除，不入库。

## 仿真流水线

```
场景生成 → 轨迹生成 → 控制器求解 → 状态传播 → 轮速计算 → 数据记录 → 指标计算 → 可视化 → 输出
```

1. **场景生成**：Latin 超立方采样，每个 seed 可复现
2. **轨迹生成**：Bernstein 多项式（Bézier）参考轨迹
3. **控制器求解**（瓶颈环节）：4 种算法
4. **状态传播**：离散动力学推进一步
5. **轮速计算**：雅可比映射到各轮
6. **数据记录**：闭环 100 步，每完成一个 (seed, alg) 存 checkpoint
7. **指标计算**：RMSE、轨迹代价 J、求解时间统计、约束违反率
8. **可视化**：逐 seed 轨迹图（3 子图）+ 汇总柱状图 + 论文对比图
9. **输出**：.mat / .csv / .png + Table II 汇总

## 快速开始

### 环境要求

- MATLAB（需 Optimization Toolbox，提供 fmincon）
- Python 3.x + hpipm_python（用于 proposed 算法）

### 运行批量对比

```matlab
cd('d:\Projects\RSS\matlab'); setup_paths;
comparison = compare_algorithms(1:20, {'proposed-3iter', 'e-lmpc', 'active-set', 'interior-point'});
```

支持断点续跑：扩大 seed 范围时已完成的 (seed, alg) 自动跳过。

### 运行论文复现

```matlab
cd('d:\Projects\RSS\matlab'); setup_paths;
comparison = paper_reproduction();
```

使用论文固定参数，仅运行三种基线算法，输出 Table II 格式汇总。

## 算法说明

| 算法 | 求解器 | 特点 |
|---|---|---|
| proposed-3iter | HPIPM（凸 QP） | RSS 凸化 + 3 次迭代更新，严格满足约束 |
| e-lmpc | fmincon SQP | 同 (P^K) 问题，MaxIter=1，对齐 RSS_sqp |
| active-set | fmincon active-set | 同 (P^K) 问题，对齐 RSS_fmincon |
| interior-point | fmincon interior-point | 同 (P^K) 问题，可违反转向锥约束 |

三种基线严格对齐原始仓库（RSS_sqp / RSS_fmincon）实现，求解相同的 (P^K) 问题，确保对比公平。

## 关键约束

- 基线控制器不含 RSS 正则项（rho）
- 转向锥约束采用论文式 (14) 双线性形式
- 状态递推：nu^{k+1} = nu^k + u^{k+1}
- 位置预测使用固定 R_psi0
- checkpoint 文件不可删除，用于断点续跑
