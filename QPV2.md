# QPV2 — H 矩阵新构造方法测试

> 本文档记录一种**基于选择矩阵的 H 构造思路**，作为当前 `construct_complete_qp_from_rss.m` 的对照方案。
> **仅文档，不改代码。**

## 1. 动机

用**选择矩阵**把所有 $\mathbf{e}_k$ 统一表达为 $\mathbf{e}_k = \mathbf{B}\,\mathbf{I}_k\,\boldsymbol{v} + \mathbf{c}_k$，然后直接做矩阵乘法 $\mathbf{e}_k^\top \mathbf{Q} \mathbf{e}_k = (\mathbf{B}\,\mathbf{I}_k\,\boldsymbol{v} + \mathbf{c}_k)^\top \mathbf{Q} (\mathbf{B}\,\mathbf{I}_k\,\boldsymbol{v} + \mathbf{c}_k)$，一次性得到 H/g/const，避免逐项展开。

$\boldsymbol{v}$ 含完整预测时域 $v_0..v_K$（含 $v_K$），从而控制正则和 RSS 正则的求和范围 $k=1..K$ 与 V1 完全一致，**V1 与 V2 严格等价**（见第 3 节）。

## 2. 新方法思路（V2）— 选择矩阵法

### 2.1 变量定义

将预测时域内的速度变量（**含 $v_0$ 和 $v_K$**）堆叠成长向量 $\boldsymbol{v}$（$3(K+1) \times 1$）：

$$
\boldsymbol{v} = \begin{bmatrix} v_0 \\ v_1 \\ v_2 \\ \vdots \\ v_K \end{bmatrix} \quad (3(K+1) \times 1)
$$

- $v_0 = v_0^{\text{cur}}$（当前车体系速度 `state_dot`，**已知量**，靠等式约束 $v_0 = v_0^{\text{cur}}$ 锁定）
- $v_l = \nu(:,l)$ for $l=0..K$（共 $K+1$ 块）
- **含 $v_K$**：对应 V1 中的 $\nu(:,K)$，使控制正则和 RSS 正则的求和范围 $k=1..K$ 与 V1 一致（含 $u_K$）
- $v_K$ 在**跟踪代价**中不出现（$\mathbf{e}_k$ 最大用到 $v_0..v_{K-1}$），但出现在**控制/RSS 正则**中（$u_K = v_K - v_{K-1}$）

```matlab
% 对应代码 (§8):
n_var = 3 * (K+1);    % 21 (V1 是 36, V2 消去了 u 但保留 v_0 和 v_K)
% v_0 = v0 作为等式约束锁定, 不分离
```

### 2.2 误差递推（$\mathbf{c}_k$ 随 $k$ 变化）

论文公式 (19) 的离散化形式（求和从 $l=0$ 开始）：

$$
\mathbf{e}_k = \boldsymbol{\xi}_{\text{cur}} - \boldsymbol{\xi}_k^{\text{ref}} + \mathbf{R}(\psi_0) \cdot \tau \cdot \sum_{l=0}^{k-1} v_l
$$

定义 $\mathbf{c}_k$（随 $k$ 变化的常数部分，**不含 $v_0$**，因为 $v_0$ 已在 $\boldsymbol{v}$ 中）：

$$
\mathbf{c}_k = \boldsymbol{\xi}_{\text{cur}} - \boldsymbol{\xi}_k^{\text{ref}} \quad (3 \times 1,\ \text{随}\ k\ \text{变化})
$$

- $\boldsymbol{\xi}_{\text{cur}} - \boldsymbol{\xi}_k^{\text{ref}}$：随 $k$ 变化（因为 $\boldsymbol{\xi}_k^{\text{ref}}$ 每步不同）
- $v_0$ 的贡献由选择矩阵 $\mathbf{I}_k$ 从 $\boldsymbol{v}$ 中提取（见 2.3）
- **$v_K$ 不出现在任何 $\mathbf{e}_k$ 中**：因求和上限为 $k-1 \le K-1$，$v_K$ 块在 $\mathbf{I}_k$ 中恒为 $\mathbf{0}$

因此误差可写成：

$$
\mathbf{e}_k = \underbrace{\mathbf{c}_k}_{\text{随}\ k\ \text{变化}} + \underbrace{\mathbf{B}\,\mathbf{I}_k\,\boldsymbol{v}}_{\text{含}\ v_0,\ \text{不含}\ v_K}
$$

```matlab
% 对应代码 (§8):
ref_idx = min(size(path, 2), step + k);
ref_xy_k  = path(1:2, ref_idx);
ref_psi_k = path(3,   ref_idx);

% c_k: 3×1, 随 k 变化 (不含 v_0, 因为 v_0 在 v 中)
%   c_k = [ current_xy - ref_xy_k ]   ← 位置误差常数 (2×1)
%         [ psi0      - ref_psi_k ]   ← 姿态误差常数 (1×1)
c_k = [current_xy - ref_xy_k;
            psi0      - ref_psi_k];
```

### 2.3 选择矩阵 $\mathbf{I}_k$

定义**选择矩阵 $\mathbf{I}_k$**（$3 \times 3(K+1)$）：

$$
\mathbf{I}_k = \begin{bmatrix} \mathbf{I} & \mathbf{I} & \cdots & \mathbf{I} & \mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} \end{bmatrix}
$$

- $\mathbf{I}$ 是 $3 \times 3$ 单位矩阵
- 前 $k$ 个块是 $\mathbf{I}$（含 $v_0$ 块，对应 $l=0..k-1$），后 $K+1-k$ 个块是 $\mathbf{0}$（**含 $v_K$ 块**）
- 作用：从 $\boldsymbol{v}$ 中选出前 $k$ 个速度块（含 $v_0$，**不含 $v_K$）并求和

$$
\mathbf{I}_k \cdot \boldsymbol{v} = v_0 + v_1 + \cdots + v_{k-1} = \sum_{l=0}^{k-1} v_l \quad (3 \times 1)
$$

```matlab
% 对应代码 (§8):
I_k = zeros(3, n_var);                          % 3×3(K+1)
for l = 0:k-1
    I_k(:, l*3+1 : (l+1)*3) = eye(3);           % 前 k 块为 I (从 v_0 开始)
end
```

### 2.4 HPIPM 求解器要求的 QP 形式

HPIPM 求解的 dense QP 标准形式（无 slack 的硬约束子集）：

$$
\begin{aligned}
\min_{\boldsymbol{v}} \quad & \tfrac{1}{2} \boldsymbol{v}^\top \mathbf{H} \boldsymbol{v} + \mathbf{g}^\top \boldsymbol{v} + \text{const} \\
\text{s.t.} \quad & \mathbf{A}\boldsymbol{v} = \mathbf{b} \\
& \tfrac{1}{2} \boldsymbol{v}^\top \mathbf{H}_{q,i} \boldsymbol{v} + \mathbf{g}_{q,i}^\top \boldsymbol{v} \le u_{q,i}
\end{aligned}
$$

- $\mathbf{H}$：Hessian 矩阵（$n_{\text{var}} \times n_{\text{var}}$，对称半正定）
- $\mathbf{g}$：一次项系数向量（$n_{\text{var}} \times 1$）
- $\text{const}$：常数项（标量，不影响最优解，仅用于目标值 obj 比较）
- $\boldsymbol{v}$：决策变量 $= [v_0; v_1; \ldots; v_K]$（$3(K+1) \times 1 = 21$ 维）

**关键约定**：HPIPM 目标中的二次项系数是 $\tfrac{1}{2}$，而我们代价展开得到的是 $\boldsymbol{v}^\top \mathbf{A} \boldsymbol{v}$（系数为 1）。为了让两者匹配，**所有二次项贡献在累加到 $\mathbf{H}$ 时都乘以 2**，使得：

$$
\boldsymbol{v}^\top \mathbf{A} \boldsymbol{v} = \tfrac{1}{2} \boldsymbol{v}^\top (2\mathbf{A}) \boldsymbol{v} \quad \Leftarrow \quad \mathbf{H} \text{ 中存的是 } 2\mathbf{A}
$$

### 2.5 RSS 论文的代价函数

论文公式 (17)-(18)，代价函数分为三项（求和范围 $k=1..K$，**含 $u_K$**，与 V1 一致）：

$$
J = \underbrace{\sum_{k=1}^{K} \mathbf{e}_k^\top \mathbf{Q} \mathbf{e}_k}_{\text{① 跟踪代价}} + \underbrace{\sum_{k=1}^{K} \mathbf{u}_k^\top \mathbf{R} \mathbf{u}_k}_{\text{② 控制正则化}} + \underbrace{\rho \sum_{k=1}^{K} \|\mathbf{u}_k - \hat{\mathbf{u}}_k\|^2}_{\text{③ RSS 强凸正则化}}
$$

参数：
- $\mathbf{Q} = \text{diag}(30, 30, 1)$：跟踪权重（位置 30，姿态 1）
- $\mathbf{R} = 0.3 \cdot \mathbf{I}_3$：控制权重
- $\rho = 0.01$：RSS 强凸参数
- $\hat{\mathbf{u}}_k$：上次迭代求得的 $\mathbf{u}_k$（已知常数）

对应到 HPIPM 形式，三项分别贡献到 H/g/const：

$$
\mathbf{H} = \mathbf{H}_{\text{track}} + \mathbf{H}_u + \mathbf{H}_\rho, \quad \mathbf{g} = \mathbf{g}_{\text{track}} + \mathbf{g}_u + \mathbf{g}_\rho, \quad \text{const} = \text{const}_{\text{track}} + \text{const}_u + \text{const}_\rho
$$

下面分别展开三项。

### 2.6 ① 跟踪代价 $\mathbf{e}_k^\top \mathbf{Q} \mathbf{e}_k$

#### 误差的仿射表达

由 §2.2、§2.3，误差 $\mathbf{e}_k$ 可写为决策变量 $\boldsymbol{v}$ 的仿射函数：

$$
\mathbf{e}_k = \underbrace{\mathbf{c}_k}_{\text{常数}} + \underbrace{\mathbf{B}\,\mathbf{I}_k\,\boldsymbol{v}}_{\text{线性项}}
$$

- $\mathbf{c}_k$：$3 \times 1$ 常数向量（随 $k$ 变化，见 §2.6.1）
- $\mathbf{B}$：$3 \times 3$ 单步转移矩阵（见 §2.6.2）
- $\mathbf{I}_k$：$3 \times 3(K+1)$ 选择矩阵（前 $k$ 块为 $\mathbf{I}$，含 $v_0$，**不含 $v_K$**）

#### 2.6.1 $\mathbf{c}_k$ 是什么

**$\mathbf{c}_k$ 是跟踪误差的"常数偏移部分"**（$3 \times 1$ 向量，随 $k$ 变化），对应论文公式 (19) 中**不含决策变量**的那一项。

**物理含义**：$\mathbf{c}_k$ 表示**当前状态与第 $k$ 步参考轨迹之间的初始误差**，由两部分组成：

$$
\mathbf{c}_k = \begin{bmatrix} \boldsymbol{\xi}_{\text{cur}}^{xy} - \boldsymbol{\xi}_k^{\text{ref},xy} \\ \psi_0 - \psi_k^{\text{ref}} \end{bmatrix} \quad \begin{matrix} \leftarrow \text{位置误差常数}\ (2 \times 1,\ \text{世界系}) \\ \leftarrow \text{姿态误差常数}\ (1 \times 1,\ \text{标量}) \end{matrix}
$$

- $\boldsymbol{\xi}_{\text{cur}}^{xy}$（$2 \times 1$）：当前世界系位置 $[\text{state}(1); \text{state}(2)]$，**已知量**（不随 $k$ 变化）
- $\boldsymbol{\xi}_k^{\text{ref},xy}$（$2 \times 1$）：第 $k$ 步参考位置 `path(1:2, ref_idx)`，**随 $k$ 变化**（贝塞尔曲线采样点）
- $\psi_0$（标量）：当前航向 `state(3)`，**已知量**（不随 $k$ 变化）
- $\psi_k^{\text{ref}}$（标量）：第 $k$ 步参考航向 `path(3, ref_idx)`，**随 $k$ 变化**

**为什么 $\mathbf{c}_k$ 随 $k$ 变化**：因为 $\boldsymbol{\xi}_k^{\text{ref},xy}$ 和 $\psi_k^{\text{ref}}$ 是参考轨迹上第 $k$ 个采样点，每一步都不同（机器人沿轨迹前进，参考点也在移动）。

**为什么 $\mathbf{c}_k$ 不含 $v_0$**：在 V2 中 $v_0$ 放进了决策变量 $\boldsymbol{v}$ 中，$v_0$ 的贡献由选择矩阵 $\mathbf{I}_k$ 从 $\boldsymbol{v}$ 中提取（见 §2.3）。$\mathbf{c}_k$ 只包含**纯已知常数**（当前状态 - 参考轨迹），不含任何决策变量。

**$\mathbf{c}_k$ 在跟踪代价中的角色**：跟踪误差 $\mathbf{e}_k = \mathbf{c}_k + \mathbf{B}\,\mathbf{I}_k\,\boldsymbol{v}$，其中：
- $\mathbf{c}_k$：常数偏移（产生常数项和一次项）
- $\mathbf{B}\,\mathbf{I}_k\,\boldsymbol{v}$：决策变量的线性组合（产生一次项和二次项）

**参考轨迹索引 `ref_idx`**：

```matlab
ref_idx = min(size(path, 2), step + k);
```

- `step`：当前全局步数（机器人在轨迹上的位置）
- `k`：预测时域内的步数（$1..K$）
- `step + k`：第 $k$ 步预测对应的参考轨迹索引
- `min(..., size(path, 2))`：防止索引越界（轨迹末端用最后一个点填充）

```matlab
% 对应代码 (§8):
ref_idx = min(size(path, 2), step + k);
ref_xy_k  = path(1:2, ref_idx);   % 第 k 步参考位置 (2×1)
ref_psi_k = path(3,   ref_idx);   % 第 k 步参考航向 (标量)

% c_k: 3×1, 随 k 变化 (不含 v_0, 因为 v_0 在 v 中)
c_k = [current_xy - ref_xy_k;
            psi0  - ref_psi_k];
```

#### 代入代价展开

利用 $\mathbf{e}_k = \mathbf{c}_k + \mathbf{B}\,\mathbf{I}_k\,\boldsymbol{v}$：

$$
\mathbf{e}_k^\top \mathbf{Q} \mathbf{e}_k = (\mathbf{c}_k + \mathbf{B}\,\mathbf{I}_k\,\boldsymbol{v})^\top \mathbf{Q} (\mathbf{c}_k + \mathbf{B}\,\mathbf{I}_k\,\boldsymbol{v})
$$

$$
= \underbrace{\mathbf{c}_k^\top \mathbf{Q} \mathbf{c}_k}_{(A)\ \text{常数项}} + \underbrace{2 \mathbf{c}_k^\top \mathbf{Q} (\mathbf{B}\,\mathbf{I}_k) \boldsymbol{v}}_{(B)\ \text{一次项}} + \underbrace{\boldsymbol{v}^\top (\mathbf{B}\,\mathbf{I}_k)^\top \mathbf{Q} (\mathbf{B}\,\mathbf{I}_k) \boldsymbol{v}}_{(C)\ \text{二次项}}
$$

#### (A) 常数项

$$
\text{const}_{\text{track}}(k) = \mathbf{c}_k^\top \mathbf{Q} \mathbf{c}_k \quad \text{(标量，随}\ k\ \text{变化)}
$$

累加到 `objective_constant`：

$$
\text{const}_{\text{track}} = \sum_{k=1}^{K} \mathbf{c}_k^\top \mathbf{Q} \mathbf{c}_k
$$

#### (B) 一次项

$$
2 \mathbf{c}_k^\top \mathbf{Q} (\mathbf{B}\,\mathbf{I}_k) \boldsymbol{v} = \big[ \underbrace{2 (\mathbf{B}\,\mathbf{I}_k)^\top \mathbf{Q} \mathbf{c}_k}_{\mathbf{g}_{\text{track}}(k)} \big]^\top \boldsymbol{v}
$$

$\mathbf{g}_{\text{track}}(k) = 2(\mathbf{B}\,\mathbf{I}_k)^\top \mathbf{Q} \mathbf{c}_k$（$3(K+1) \times 1$ 向量，随 $k$ 变化），累加：

$$
\mathbf{g}_{\text{track}} = \sum_{k=1}^{K} 2 (\mathbf{B}\,\mathbf{I}_k)^\top \mathbf{Q} \mathbf{c}_k = \sum_{k=1}^{K} 2 \mathbf{I}_k^\top \mathbf{B}^\top \mathbf{Q} \mathbf{c}_k
$$

#### (C) 二次项

$$
\boldsymbol{v}^\top (\mathbf{B}\,\mathbf{I}_k)^\top \mathbf{Q} (\mathbf{B}\,\mathbf{I}_k) \boldsymbol{v} = \tfrac{1}{2} \boldsymbol{v}^\top \underbrace{[2 (\mathbf{B}\,\mathbf{I}_k)^\top \mathbf{Q} (\mathbf{B}\,\mathbf{I}_k)]}_{\mathbf{H}_{\text{track}}(k)} \boldsymbol{v}
$$

$\mathbf{H}_{\text{track}}(k) = 2 (\mathbf{B}\,\mathbf{I}_k)^\top \mathbf{Q} (\mathbf{B}\,\mathbf{I}_k)$（$3(K+1) \times 3(K+1)$，前乘 2 抵消 HPIPM 的 0.5），累加：

$$
\mathbf{H}_{\text{track}} = 2 \sum_{k=1}^{K} (\mathbf{B}\,\mathbf{I}_k)^\top \mathbf{Q} (\mathbf{B}\,\mathbf{I}_k) = 2 \sum_{k=1}^{K} \mathbf{I}_k^\top \underbrace{(\mathbf{B}^\top \mathbf{Q} \mathbf{B})}_{\mathbf{Q}_B} \mathbf{I}_k = 2 \sum_{k=1}^{K} \mathbf{I}_k^\top \mathbf{Q}_B \mathbf{I}_k
$$

其中 $\mathbf{Q}_B = \mathbf{B}^\top \mathbf{Q} \mathbf{B}$（$3 \times 3$ 等效权重，见 §2.6.2）。

#### 2.6.2 $\mathbf{B}$ 矩阵的定义

$\mathbf{B}$ 是 $3 \times 3$ 单步转移矩阵，编码"车体系速度→世界系位移"的物理含义：

- **位置部分**：$\mathbf{B}_{\text{pos}} = \mathbf{R}(\psi_0) \cdot \tau$（车体系速度→世界系位移，$2 \times 2$）
- **姿态部分**：$B_\psi = \tau$（标量，无需旋转）

$\mathbf{B}$ 的形式（$3 \times 3$ 分块对角）：

$$
\mathbf{B} = \begin{bmatrix} \mathbf{R}(\psi_0) \cdot \tau & \mathbf{0} \\ \mathbf{0} & \tau \end{bmatrix}
$$

代入 $\mathbf{Q}_B = \mathbf{B}^\top \mathbf{Q} \mathbf{B}$：

$$
\mathbf{Q}_B = \mathbf{B}^\top \mathbf{Q} \mathbf{B} = \begin{bmatrix} (\mathbf{R}\tau)^\top & \mathbf{0} \\ \mathbf{0} & \tau \end{bmatrix} \begin{bmatrix} 30\mathbf{I} & \mathbf{0} \\ \mathbf{0} & 1 \end{bmatrix} \begin{bmatrix} \mathbf{R}\tau & \mathbf{0} \\ \mathbf{0} & \tau \end{bmatrix} = \begin{bmatrix} 30 \tau^2 \cdot \mathbf{I} & \mathbf{0} \\ \mathbf{0} & \tau^2 \end{bmatrix}
$$

（利用 $\mathbf{R}^\top \mathbf{R} = \mathbf{I}$）

- $\mathbf{Q}_B$ 的位置块 $= 30 \tau^2 \cdot \mathbf{I}$（$2 \times 2$）
- $\mathbf{Q}_B$ 的姿态块 $= \tau^2$（$1 \times 1$）

```matlab
% 对应代码 (§8):
B = [R_psi0 * dt,   zeros(2,1);           % 3×3 分块对角
         zeros(1,2),       dt    ];
Q = diag([w_pos, w_pos, w_psi]);          % diag(30, 30, 1)
Q_B   = B' * Q * B;               % [30*dt^2*I, 0; 0, dt^2]
```

#### $\mathbf{H}_{\text{track}}$ 的块结构

$\mathbf{I}_k^\top \mathbf{Q}_B \mathbf{I}_k$ 是 $3(K+1) \times 3(K+1)$ 矩阵，**非零块只在前 $k$ 个 $v$ 块**（块 0 到块 $k-1$）：

$$
\mathbf{I}_k^\top \mathbf{Q}_B \mathbf{I}_k = \begin{bmatrix} \mathbf{Q}_B & \cdots & \mathbf{Q}_B & \mathbf{0} & \cdots & \mathbf{0} \\ \vdots & \ddots & \vdots & \vdots & & \vdots \\ \mathbf{Q}_B & \cdots & \mathbf{Q}_B & \mathbf{0} & \cdots & \mathbf{0} \\ \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} \end{bmatrix} \begin{matrix} \left.\begin{matrix} v_0 \\ \vdots \\ v_{k-1} \end{matrix}\right\} k \times k\ \text{个}\ \mathbf{Q}_B\ \text{块} \\ \left.\begin{matrix} v_k \\ \vdots \\ v_K \end{matrix}\right\} K+1-k\ \text{块}\ \mathbf{0} \end{matrix}
$$

对 $k=1..K$ 求和后，跟踪代价 $\mathbf{H}$ 块结构为（$\boldsymbol{v}$ 有 $K+1$ 块，索引 $i,j \in 0..K$）：

$$
\frac{\mathbf{H}_{\text{track}}}{2} = \begin{bmatrix} K\mathbf{Q}_B & (K-1)\mathbf{Q}_B & \cdots & 1\cdot\mathbf{Q}_B & \mathbf{0} \\ (K-1)\mathbf{Q}_B & (K-1)\mathbf{Q}_B & \cdots & 1\cdot\mathbf{Q}_B & \mathbf{0} \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 1\cdot\mathbf{Q}_B & 1\cdot\mathbf{Q}_B & \cdots & 1\cdot\mathbf{Q}_B & \mathbf{0} \\ \mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} \end{bmatrix} \begin{matrix} v_0 \\ v_1 \\ \vdots \\ v_{K-1} \\ v_K \end{matrix}
$$

- **$v_0$ 行/列**：$K \cdot \mathbf{Q}_B$（被所有 $K$ 个 $\mathbf{e}_k$ 使用）
- **$v_{K-1}$ 行/列**：$1 \cdot \mathbf{Q}_B$（仅被 $\mathbf{e}_K$ 使用）
- **$v_K$ 行/列**：$\mathbf{0}$（不被任何 $\mathbf{e}_k$ 使用，跟踪代价不依赖 $v_K$）

> $v_K$ 在 $\mathbf{H}_{\text{track}}$ 中全为 $\mathbf{0}$，但其行/列会在 §2.7 控制正则和 §2.8 RSS 正则中获得非零贡献。

#### 代码对应

```matlab
% 对应代码 (§8) — 三部分一次性计算:
const_track = const_track + c_k' * Q * c_k;                       % (A) 常数
g_track     = g_track     + 2 * (B*I_k)' * Q * c_k;                % (B) 一次
H_track     = H_track     + 2 * (B*I_k)' * Q * (B*I_k);            % (C) 二次
```

### 2.7 ② 控制正则化 $\mathbf{u}_k^\top \mathbf{R} \mathbf{u}_k$

#### $\mathbf{u}_k$ 的差分表达

由论文公式 (9) $\nu_{k+1} = \nu_k + u_{k+1}$ 反推：

$$
\mathbf{u}_k = v_k - v_{k-1} \quad \text{(差分关系)}
$$

定义**差分矩阵 $\mathbf{D}_k$**（$3 \times 3(K+1)$），使得 $\mathbf{u}_k = \mathbf{D}_k \boldsymbol{v}$：

$$
\mathbf{D}_k = \begin{bmatrix} \mathbf{0} & \mathbf{0} & \mathbf{0} & \cdots & -\mathbf{I} & \mathbf{I} & \cdots & \mathbf{0} \end{bmatrix}
$$

- 第 $k-1$ 块为 $-\mathbf{I}$，第 $k$ 块为 $+\mathbf{I}$，其余为 $\mathbf{0}$
- 因 $v_0$ 在 $\boldsymbol{v}$ 中，$\mathbf{u}_k = \mathbf{D}_k \boldsymbol{v}$ **完全由决策变量线性表达，无常数偏移**

```matlab
% 对应代码 (§8):
D_k = zeros(3, n_var);
D_k(:, k*3+1 : (k+1)*3) =  eye(3);   % v_k     块 = +I
D_k(:, (k-1)*3+1 : k*3) = -eye(3);   % v_{k-1} 块 = -I
```

#### 代价展开

因 $\mathbf{u}_k$ 完全由 $\boldsymbol{v}$ 线性表达（无常数偏移），代入 $\mathbf{u}_k^\top \mathbf{R} \mathbf{u}_k$：

$$
\mathbf{u}_k^\top \mathbf{R} \mathbf{u}_k = (\mathbf{D}_k \boldsymbol{v})^\top \mathbf{R} (\mathbf{D}_k \boldsymbol{v}) = \boldsymbol{v}^\top \mathbf{D}_k^\top \mathbf{R} \mathbf{D}_k \boldsymbol{v} \quad \Leftarrow \text{只有二次项}
$$

#### (C) 二次项

$\mathbf{D}_k^\top \mathbf{R} \mathbf{D}_k$ 是 $3(K+1) \times 3(K+1)$ 矩阵，**只在 4 个 $3 \times 3$ 块上非零**：

$$
\begin{aligned}
\text{块}\ (k, k) &: \quad \mathbf{I}^\top \mathbf{R} \mathbf{I} = +\mathbf{R} \quad \Leftarrow v_k \text{对角} \\
\text{块}\ (k, k-1) &: \quad \mathbf{I}^\top \mathbf{R} (-\mathbf{I}) = -\mathbf{R} \quad \Leftarrow v_k \text{与}\ v_{k-1} \text{交叉} \\
\text{块}\ (k-1, k) &: \quad (-\mathbf{I})^\top \mathbf{R} \mathbf{I} = -\mathbf{R} \quad \Leftarrow v_{k-1} \text{与}\ v_k \text{交叉} \\
\text{块}\ (k-1, k-1) &: \quad (-\mathbf{I})^\top \mathbf{R} (-\mathbf{I}) = +\mathbf{R} \quad \Leftarrow v_{k-1} \text{对角}
\end{aligned}
$$

对 $k=1..K$ 求和（前乘 2 抵消 HPIPM 的 0.5），得到 $(K+1) \times (K+1)$ **三对角块矩阵**：

$$
\frac{\mathbf{H}_u}{2} = \begin{bmatrix} \mathbf{R} & -\mathbf{R} & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} \\ -\mathbf{R} & 2\mathbf{R} & -\mathbf{R} & \cdots & \mathbf{0} & \mathbf{0} \\ \mathbf{0} & -\mathbf{R} & 2\mathbf{R} & \cdots & \mathbf{0} & \mathbf{0} \\ \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\ \mathbf{0} & \mathbf{0} & \mathbf{0} & \cdots & 2\mathbf{R} & -\mathbf{R} \\ \mathbf{0} & \mathbf{0} & \mathbf{0} & \cdots & -\mathbf{R} & \mathbf{R} \end{bmatrix} \begin{matrix} v_0 \\ v_1 \\ v_2 \\ \vdots \\ v_{K-1} \\ v_K \end{matrix}
$$

- **对角块**：首尾（$v_0$、$v_K$）为 $\mathbf{R}$（只被一个 $\mathbf{D}_k$ 覆盖），中间为 $2\mathbf{R}$（被两个相邻 $\mathbf{D}_k$ 覆盖）
- **相邻非对角块**：$-\mathbf{R}$（差分交叉项）
- **非相邻块**：$\mathbf{0}$

#### (B) 一次项 / (A) 常数项

**均为 0**（因 $\mathbf{u}_k$ 完全由 $\boldsymbol{v}$ 线性表达，无常数偏移）：

$$
\mathbf{g}_u = \mathbf{0}, \quad \text{const}_u = 0
$$

#### 代码对应

```matlab
% 对应代码 (§8) — 只有二次项:
H_u = H_u + 2 * D_k' * R * D_k;    % (C) 二次
% g_u = 0, const_u = 0 (无常数偏移)
```

### 2.8 ③ RSS 强凸正则化 $\rho \|\mathbf{u}_k - \hat{\mathbf{u}}_k\|^2$

#### 展开形式

由 §2.7 知 $\mathbf{u}_k = \mathbf{D}_k \boldsymbol{v}$，所以：

$$
\mathbf{u}_k - \hat{\mathbf{u}}_k = \mathbf{D}_k \boldsymbol{v} - \hat{\mathbf{u}}_k
$$

其中 $\hat{\mathbf{u}}_k$ 是上次迭代解（**已知常数向量**）。代入代价：

$$
\rho \|\mathbf{u}_k - \hat{\mathbf{u}}_k\|^2 = \rho (\mathbf{D}_k \boldsymbol{v} - \hat{\mathbf{u}}_k)^\top (\mathbf{D}_k \boldsymbol{v} - \hat{\mathbf{u}}_k) = \rho \big[ \underbrace{\hat{\mathbf{u}}_k^\top \hat{\mathbf{u}}_k}_{\text{常数}} - \underbrace{2 \hat{\mathbf{u}}_k^\top \mathbf{D}_k \boldsymbol{v}}_{\text{一次项}} + \underbrace{\boldsymbol{v}^\top \mathbf{D}_k^\top \mathbf{D}_k \boldsymbol{v}}_{\text{二次项}} \big]
$$

#### (C) 二次项

注意此处权重是 $\rho \mathbf{I}$（不是 $\mathbf{R}$），所以 Hessian 贡献与 §2.7 同构但用 $\mathbf{I}$ 替换 $\mathbf{R}$：

$$
\mathbf{H}_\rho = 2 \rho \sum_{k=1}^{K} \mathbf{D}_k^\top \mathbf{D}_k
$$

三对角块结构（与 $\mathbf{H}_u$ 同形，把 $\mathbf{R}$ 换成 $\rho \mathbf{I}$）：

$$
\frac{\mathbf{H}_\rho}{2} = \rho \begin{bmatrix} \mathbf{I} & -\mathbf{I} & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} \\ -\mathbf{I} & 2\mathbf{I} & -\mathbf{I} & \cdots & \mathbf{0} & \mathbf{0} \\ \mathbf{0} & -\mathbf{I} & 2\mathbf{I} & \cdots & \mathbf{0} & \mathbf{0} \\ \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\ \mathbf{0} & \mathbf{0} & \mathbf{0} & \cdots & 2\mathbf{I} & -\mathbf{I} \\ \mathbf{0} & \mathbf{0} & \mathbf{0} & \cdots & -\mathbf{I} & \mathbf{I} \end{bmatrix} \begin{matrix} v_0 \\ v_1 \\ v_2 \\ \vdots \\ v_{K-1} \\ v_K \end{matrix}
$$

#### (B) 一次项

$$
\mathbf{g}_\rho = -2 \rho \sum_{k=1}^{K} \mathbf{D}_k^\top \hat{\mathbf{u}}_k
$$

逐 $k$ 展开（$\mathbf{D}_k^\top \hat{\mathbf{u}}_k$ 在每个块上的贡献）：

$$
\begin{aligned}
k=1 &: \quad \mathbf{D}_1^\top \hat{\mathbf{u}}_1 \to v_0 \text{块}: -\hat{\mathbf{u}}_1,\ v_1 \text{块}: +\hat{\mathbf{u}}_1 \\
k=2..K-1 &: \quad \mathbf{D}_k^\top \hat{\mathbf{u}}_k \to v_{k-1} \text{块}: -\hat{\mathbf{u}}_k,\ v_k \text{块}: +\hat{\mathbf{u}}_k \\
k=K &: \quad \mathbf{D}_K^\top \hat{\mathbf{u}}_K \to v_{K-1} \text{块}: -\hat{\mathbf{u}}_K,\ v_K \text{块}: +\hat{\mathbf{u}}_K
\end{aligned}
$$

即 $\mathbf{g}_\rho$ 在每个块 $j$ 上的总贡献（$j=0..K$）：

$$
\begin{aligned}
\text{块}\ 0\ (v_0) &: \quad -2\rho \hat{\mathbf{u}}_1 \quad \text{(仅}\ k=1\ \text{贡献)} \\
\text{块}\ j\ (1..K-1) &: \quad +2\rho \hat{\mathbf{u}}_j - 2\rho \hat{\mathbf{u}}_{j+1} \quad \text{(}k=j \text{的} +\hat{\mathbf{u}}_j \text{和} k=j+1 \text{的} -\hat{\mathbf{u}}_{j+1}\text{)} \\
\text{块}\ K\ (v_K) &: \quad +2\rho \hat{\mathbf{u}}_K \quad \text{(仅}\ k=K\ \text{贡献)}
\end{aligned}
$$

#### (A) 常数项

$$
\text{const}_\rho = \rho \sum_{k=1}^{K} \|\hat{\mathbf{u}}_k\|^2
$$

#### 代码对应

```matlab
% 对应代码 (§8) — 三部分:
H_rho     = H_rho     + 2 * rho * D_k' * D_k;                    % (C) 二次
g_rho     = g_rho     - 2 * rho * D_k' * u_hat(:, k);             % (B) 一次
const_rho = const_rho +       rho * u_hat(:, k)' * u_hat(:, k);   % (A) 常数
```

## 3. V1 vs V2 理论等价性分析

### 3.1 核心结论

**V1 与 V2 严格等价**：V2 的 $\boldsymbol{v}$ 含 $v_0..v_K$（$3(K+1)$ 维），控制正则和 RSS 正则的求和范围 $k=1..K$（含 $u_K$）与 V1 完全一致。两者代价函数在 V1 等式约束消元后**完全相同**，最优解、目标值、KKT 系统均严格一致。

### 3.2 代价函数逐项对比

$$
\begin{aligned}
\text{V1: } J &= \sum_{k=1}^{K} [\mathbf{e}_k^\top \mathbf{Q} \mathbf{e}_k] + \sum_{k=1}^{K} [w_{\text{ctrl}} \cdot \|\mathbf{u}_k\|^2] + \sum_{k=1}^{K} [\rho \cdot \|\mathbf{u}_k - \hat{\mathbf{u}}_k\|^2] \\
& \hspace{3.5cm} \uparrow_{k=1..K\ (\text{含}\ u_K)} \hspace{1.2cm} \uparrow_{k=1..K\ (\text{含}\ u_K)} \\
\text{V2: } J &= \sum_{k=1}^{K} [\mathbf{e}_k^\top \mathbf{Q} \mathbf{e}_k] + \sum_{k=1}^{K} [w_{\text{ctrl}} \cdot \|\mathbf{u}_k\|^2] + \sum_{k=1}^{K} [\rho \cdot \|\mathbf{u}_k - \hat{\mathbf{u}}_k\|^2] \\
& \hspace{3.5cm} \uparrow_{k=1..K\ (\text{含}\ u_K)} \hspace{1.2cm} \uparrow_{k=1..K\ (\text{含}\ u_K)}
\end{aligned}
$$

- **跟踪代价**：完全相同（$k=1..K$，$\mathbf{e}_k$ 不依赖 $v_K$）
- **控制正则**：完全相同（$k=1..K$，含 $u_K = v_K - v_{K-1}$）
- **RSS 正则**：完全相同（$k=1..K$，含 $u_K - \hat{\mathbf{u}}_K$）

> 与早期版本（$\boldsymbol{v}$ 不含 $v_K$）的差异：早期 V2 缺 $u_K$ 的正则贡献，仅为近似等价；当前版本 $\boldsymbol{v}$ 含 $v_K$ 后，V1 与 V2 严格等价。

### 3.3 约束消元关系（完全等价）

V1 的等式约束 $\nu_{k+1} = \nu_k + u_{k+1}$ 可以显式求解 $\mathbf{u}$：

$$
\begin{aligned}
\mathbf{u}(:,1)   &= \nu(:,1) - v_0^{\text{cur}} \quad &\text{(初始条件)} \\
\mathbf{u}(:,k+1) &= \nu(:,k+1) - \nu(:,k) \quad &\text{(递推, } k=1..K-1\text{)}
\end{aligned}
$$

代入 V1 代价，消去 $\mathbf{u}$ 后只剩 $\nu$ 作为决策变量。**V2 的 $\boldsymbol{v}$ 正是消元后的 $\nu$**：
- V2 的 $v_0$ ↔ V1 的 $v_0^{\text{cur}}$（已知常数，V2 靠等式约束锁定）
- V2 的 $v_k$ ↔ V1 的 $\nu(:,k)$ for $k=1..K$（V2 含 $v_K$ ↔ V1 的 $\nu(:,K)$）
- V2 的 $\mathbf{u}_k = v_k - v_{k-1}$ ↔ V1 的 $\mathbf{u}(:,k) = \nu(:,k) - \nu(:,k-1)$

消元后 V1 的决策变量维度 18（$\nu_1..\nu_K$），V2 的自由变量维度 18（$v_1..v_K$，$v_0$ 被等式约束锁定）。**两者完全等价**。

### 3.4 $\mathbf{H}$ 矩阵的关系（非逐元素相等，但经约束消元后严格等价）

V1 的 $\mathbf{H}$（$36 \times 36$，$\mathbf{u}$ 和 $\nu$ 独立）：

$$
\mathbf{H}_{V1} = \begin{bmatrix} 2(w_{\text{ctrl}}+\rho) \cdot \mathbf{I} & \mathbf{0} \\ \mathbf{0} & \mathbf{H}_{\text{track,V1}} \end{bmatrix} \begin{matrix} \leftarrow \mathbf{u} \text{块纯对角}\ (18 \times 18) \\ \leftarrow \nu \text{块跨阶段}\ (18 \times 18, \nu_1..\nu_K) \end{matrix}
$$

V2 的 $\mathbf{H}$（$21 \times 21$，$\boldsymbol{v}$ 含 $v_0..v_K = 3(K+1) = 21$ 维）：

$$
\mathbf{H}_{V2} = \underbrace{\mathbf{H}_{\text{track}}}_{\text{跟踪项}\ (v_0 \text{行/列} = K\mathbf{Q}_B,\ v_K \text{行/列} = \mathbf{0})} + \underbrace{2 \sum_{k=1}^{K} \mathbf{D}_k^\top \mathbf{R} \mathbf{D}_k}_{\text{控制项消元后贡献}\ ((K+1) \times (K+1)\ \text{三对角块})} + \underbrace{2\rho \sum_{k=1}^{K} \mathbf{D}_k^\top \mathbf{D}_k}_{\text{RSS 项消元后贡献}\ ((K+1) \times (K+1)\ \text{三对角块})}
$$

**关键点**：
- $\mathbf{H}_{V2} \ne$ $\mathbf{H}_{V1}$ 的任何子矩阵（两者维度不同：$21 \times 21$ vs $36 \times 36$）
- $\mathbf{H}_{V2} = \mathbf{T}^\top \mathbf{H}_{V1} \mathbf{T}$ **严格成立**，其中 $\mathbf{T}$ 是 $36 \times 21$ 的消元矩阵（见 §3.5）
- $v_0$ 行/列在 $\mathbf{H}_{V2}$ 中有非零贡献（$K\mathbf{Q}_B + \mathbf{R} + \rho \mathbf{I}$），但被等式约束 $v_0 = v_0^{\text{cur}}$ 锁定
- $v_K$ 行/列在 $\mathbf{H}_{\text{track}}$ 中为 $\mathbf{0}$，但在 $\mathbf{H}_u$ 和 $\mathbf{H}_\rho$ 中为 $\mathbf{R}$ 和 $\rho \mathbf{I}$（被 $\mathbf{D}_K$ 覆盖一次）
- 两者**不能逐元素对比**，需通过约束消元矩阵 $\mathbf{T}$ 验证严格等价性

### 3.5 等价性验证方法

由于 $\mathbf{H}$ 维度不同，采用以下三种验证：

**(1) 最优解一致性**（最直接）：

$$
\begin{aligned}
\text{V1 解: } \boldsymbol{x}^* &= [\mathbf{u}^*(1..18);\ \nu^*(1..18)] \quad (36\text{维}) \\
\text{V2 解: } \boldsymbol{v}^* &= [v_0^*;\ v_1^*;\ \ldots;\ v_K^*] \quad (21\text{维}, \boldsymbol{v} \text{含}\ v_0 \text{和}\ v_K) \\
\text{验证: } v_0^* &\stackrel{!}{=} v_0^{\text{cur}} \quad \text{(等式约束锁定)} \\
v_k^* &\stackrel{!}{=} \nu^*(:,k) \quad \text{for}\ k=1..K \quad \text{(严格相等)} \\
\mathbf{u}_k^* &\stackrel{!}{=} v_k^* - v_{k-1}^* \quad \text{for}\ k=1..K \quad \text{(差分关系)}
\end{aligned}
$$

**(2) 目标值一致性**：

$$
\begin{aligned}
\text{obj}_{V1} &= \tfrac{1}{2} \boldsymbol{x}^{*\top} \mathbf{H}_{V1} \boldsymbol{x}^* + \mathbf{g}_{V1}^\top \boldsymbol{x}^* + \text{const}_{V1} \\
\text{obj}_{V2} &= \tfrac{1}{2} \boldsymbol{v}^{*\top} \mathbf{H}_{V2} \boldsymbol{v}^* + \mathbf{g}_{V2}^\top \boldsymbol{v}^* + \text{const}_{V2} \\
\text{验证: } &\text{obj}_{V1} \stackrel{!}{=} \text{obj}_{V2} \quad \text{(严格相等, 无差异)}
\end{aligned}
$$

**(3) KKT 系统一致性**（最严格）：

$$
\text{V1 的 KKT:} \quad \begin{bmatrix} \mathbf{H}_{V1} & \mathbf{A}_{V1}^\top \\ \mathbf{A}_{V1} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \boldsymbol{x}^* \\ \boldsymbol{\lambda}^* \end{bmatrix} = \begin{bmatrix} -\mathbf{g}_{V1} \\ \mathbf{b}_{V1} \end{bmatrix}
$$

$$
\text{V2 的 KKT:} \quad \begin{bmatrix} \mathbf{H}_{V2} & \mathbf{A}_{V2}^\top \\ \mathbf{A}_{V2} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \boldsymbol{v}^* \\ \boldsymbol{\mu}^* \end{bmatrix} = \begin{bmatrix} -\mathbf{g}_{V2} \\ \mathbf{b}_{V2} \end{bmatrix}, \quad (\mathbf{A}_{V2} = [\mathbf{I}, \mathbf{0}, \ldots, \mathbf{0}],\ \mathbf{b}_{V2} = v_0^{\text{cur}},\ \text{锁定}\ v_0)
$$

V1 经约束消元后**严格退化为 V2**。消元矩阵 $\mathbf{T}$（$36 \times 21$）满足：

$$
\mathbf{H}_{V2} = \mathbf{T}^\top \mathbf{H}_{V1} \mathbf{T} \quad \text{(严格相等)}, \quad \mathbf{g}_{V2} = \mathbf{T}^\top \mathbf{g}_{V1} \quad \text{(严格相等)}, \quad \text{const}_{V2} = \text{const}_{V1} \quad \text{(严格相等)}
$$

其中 $\mathbf{T}$ 的构造（V2 的 $\boldsymbol{v}$ ↔ V1 的 $\boldsymbol{x} = [\mathbf{u}; \nu]$）：
- $\mathbf{u}$ 块（$k=1..K$）：$\mathbf{u}_k = v_k - v_{k-1}$，即 $\mathbf{T}$ 的 $\mathbf{u}$ 行在 $\boldsymbol{v}$ 的第 $k$ 块为 $\mathbf{I}$、第 $k-1$ 块为 $-\mathbf{I}$
- $\nu$ 块（$k=1..K$）：$\nu_k = v_k$，即 $\mathbf{T}$ 的 $\nu$ 行在 $\boldsymbol{v}$ 的第 $k$ 块为 $\mathbf{I}$

### 3.6 二次约束的等价性

V1 的二次约束（轮速 SOC + 转向锥）依赖 $\nu(:,k)$，V2 的 $\boldsymbol{v}$ 含 $v_0..v_K$：

- **轮速约束** $\|\mathbf{H}_n \cdot \nu(:,k)\| \le v_{\max}$：在 $\boldsymbol{v}$ 的第 $k$ 块（$v_k$，$k=1..K$）上构造，**V2 与 V1 完全一致**（V1 $k=1..K$，V2 $k=1..K$，均含 $v_K$ ↔ $\nu_K$）
- **转向锥约束** $C^k_{i,n}(\mathbf{u}, \hat{\mathbf{u}}) \le 0$：V1 中依赖 $\mathbf{u}$ 和 $\hat{\nu}$，V2 中 $\mathbf{u} = \mathbf{D}\boldsymbol{v}$，需将 $\mathbf{u}$ 替换为 $\boldsymbol{v}$ 的差分，约束结构稍变但**数学等价**

### 3.7 数值预期

如果 V2 实现正确（$\boldsymbol{v}$ 含 $v_0..v_K$，$3(K+1)=21$ 维）：
- **最优解**：$v_0^*$ 应严格等于 $v_0^{\text{cur}}$（等式约束锁定），$v_k^*$（$k=1..K$）应**严格等于** V1 的 $\nu^*(:,k)$
- **目标值**：$\text{obj}_{V2}$ 应**严格等于** V1 的 $\text{obj}_{V1}$（无 $u_K$ 差异）
- **求解时间**：V2 消去了 15 条等式约束（$18 \to 3$），决策变量维度从 36 降为 21，**预期求解更快**
- **RMSE/J_total**：应与基准 `[proposed-3iter: paper_fixed]` 严格一致（RMSE=0.036793, J_total=13.3838）

## 4. 待讨论的问题

1. ~~**$\mathbf{c}_k$ 是否随 $k$ 变化**~~：**已解决**。见 §2.2，$\mathbf{c}_k = \boldsymbol{\xi}_{\text{cur}} - \boldsymbol{\xi}_k^{\text{ref}}$，含 $\boldsymbol{\xi}_k^{\text{ref}}$，随 $k$ 变化（$v_0$ 在 $\boldsymbol{v}$ 中，不进 $\mathbf{c}_k$）。
2. ~~**$\mathbf{B}$ 的确切形式**~~：**已解决**。见 §2.6.2，$\mathbf{B} = \begin{bmatrix} \mathbf{R}(\psi_0) \cdot \tau & \mathbf{0} \\ \mathbf{0} & \tau \end{bmatrix}$（$3 \times 3$ 分块对角），对应位置旋转 + 姿态标量积分。
3. ~~**$v_0$ 的处理**~~：**已解决**。见 §2.1，$v_0 = v_0^{\text{cur}}$ 放进 $\boldsymbol{v}$ 中（$\boldsymbol{v}$ 含 $v_0..v_K$，$3(K+1)=21$ 维），靠 3 条等式约束 $v_0 = v_0^{\text{cur}}$ 锁定。
4. ~~**$\mathbf{u}_k^\top \mathbf{R} \mathbf{u}_k$ 和 $\rho \|\mathbf{u} - \hat{\mathbf{u}}\|^2$ 的处理**~~：**已解决**。见 §2.7、§2.8，通过差分矩阵 $\mathbf{D}_k$ 表达 $\mathbf{u}_k = \mathbf{D}_k \boldsymbol{v}$（$v_0$ 在 $\boldsymbol{v}$ 中，无常数偏移），Hessian 退化为 $\boldsymbol{v}$ 块三对角块结构。
5. ~~**$u_K$ 的处理**~~：**已解决**。$\boldsymbol{v}$ 含 $v_K$，$\mathbf{u}_K = v_K - v_{K-1}$ 自然进入差分表达，$k=1..K$ 与 V1 一致。
6. ~~**数值验证**~~：**已分析**。见第 3 节，V1 与 V2 严格等价（约束消元关系），$\mathbf{H}$ 不可逐元素对比，需通过最优解/目标值/KKT 系统验证。

## 5. 预期测试结果

基准为 `[proposed-3iter: paper_fixed]`：
- RMSE = 0.036793
- medianSolveTime = NaNs
- J_total = 13.3838
- validSteps = 100/100

如果 V2 正确，H/g/const 在**经等式约束消元后**应与 V1 严格等价（V1 的 $\mathbf{u}$ 块对角 ↔ V2 的 $\boldsymbol{v}$ 块三对角是消元前后的不同形态）。求解结果（最优解 $\boldsymbol{x}$、目标值 obj_value）应**完全一致**。

## 6. 待办

- [x] 明确 $\mathbf{B}$ 矩阵的确切形式（见 §2.6.2）
- [x] 明确 $\mathbf{e}_0$ 的定义（随 $k$ 变化，含 $\boldsymbol{\xi}_k^{\text{ref}}$，见 §2.2）
- [x] 处理 $v_0 = v_0^{\text{cur}}$ 的已知量分离（见 §2.1、§2.3）
- [x] 展开 $\mathbf{u}_k^\top \mathbf{R} \mathbf{u}_k$ 和 $\rho \|\mathbf{u} - \hat{\mathbf{u}}\|^2$ 的选择矩阵表达（见 §2.7、§2.8，用差分矩阵 $\mathbf{D}_k$）
- [x] 决定 $u_K$ 的处理方案（$v_K$ 放进 $\boldsymbol{v}$ 中）
- [x] 同步修订 §2.1/§2.3/§2.7/§2.8 的维度与 $k$ 范围（已完成，含 $v_K$）
- [ ] 实现 V2 代码（仅测试，不替换 V1）
- [ ] 对比 V1 和 V2 的 H/g/const（需先做等式约束消元，不能直接逐元素对比）
- [ ] 对照基准结果

## 7. V2 代码实现

> 以下为 V2 的完整 MATLAB 伪代码（仅文档，不替换 V1）。
> 接口与 `construct_complete_qp_from_rss.m` 对齐，便于对照测试。
> 转向锥约束暂标注 TODO，待后续实现。

```matlab
function qp_problem = construct_complete_qp_from_rss_v2(path, step, v0, state, u_hat, params)
% CONSTRUCT_COMPLETE_QP_FROM_RSS_V2
% V2: 基于选择矩阵的 H/g/const 构造方法 (与 V1 严格等价)
%
% 特点:
%   - 决策变量 v = [v_0; v_1; ...; v_K] (3(K+1) 维), u 通过差分消元
%   - v_0 靠等式约束 v_0 = v0 锁定 (3 条)
%   - v_K 在 v 中 (参与控制/RSS 正则 k=1..K, 不参与跟踪代价)
%   - H 用选择矩阵 I_k 和差分矩阵 D_k 构造, 无逐 k 双重循环

    %% =====================================================
    % 参数提取
    % ======================================================
    K = 6;
    dt = params.dt;
    phidotmax = params.phidotmax;
    vimax = params.vimax;
    wheel_pos = params.wheel_pos;
    num_wheels = size(wheel_pos, 1);

    current_xy = [state(1); state(2)];
    psi0 = state(3);
    R_psi0 = [cos(psi0), -sin(psi0); sin(psi0), cos(psi0)];

    w_pos = 30;  w_psi = 1;  w_control = 0.3;  rho = 0.01;

    % 特征矩阵 H_n (论文公式 3)
    Hn = cell(1, num_wheels);
    for n = 1:num_wheels
        Hn{n} = [1, 0, -wheel_pos(n,2); 0, 1, wheel_pos(n,1)];
    end
    % 注: 论文无 M 显式符号, H_n'·H_n 在需要处直接计算

    delta_theta = dt * phidotmax;

    %% =====================================================
    % 决策变量 (2.1): v = [v_0; v_1; ...; v_K]  (3(K+1) 维)
    % v_0 放进 v 中, 靠等式约束锁定; v_K 放进 v 中, 参与控制/RSS 正则
    % ======================================================
    n_var = 3 * (K+1);    % 21 (V1 是 36, V2 消去了 u 但保留 v_0 和 v_K)

    %% =====================================================
    % 1. 跟踪代价 (2.2-2.7): H_track, g_track, const_track
    % =====================================================

    % B 矩阵 (2.8): 3×3, 编码车体系→世界系位移
    B = [R_psi0 * dt,   zeros(2,1);
             zeros(1,2),       dt    ];

    % Q 矩阵 (3×3)
    Q = diag([w_pos, w_pos, w_psi]);

    % Q_B = B'·Q·B (2.8): 3×3 等效权重
    Q_B = B' * Q * B;
    % 结果: [30*dt^2*I, 0; 0, dt^2]

    H_track = zeros(n_var, n_var);
    g_track = zeros(n_var, 1);
    const_track = 0;

    for k = 1:K    % k=1 时 I_1 只选 v_0 块 (被等式约束锁定), e_1 为纯常数, 不影响最优解
        % 参考轨迹 (随 k 变化)
        ref_idx = min(size(path, 2), step + k);
        ref_xy_k  = path(1:2, ref_idx);
        ref_psi_k = path(3,   ref_idx);

        % c_k (2.2): 随 k 变化的常数部分, 3×1 (不含 v_0, 因为 v_0 在 v 中)
        %   c_k = [ current_xy - ref_xy_k ]   ← 位置误差常数 (2×1)
        %         [ psi0      - ref_psi_k ]   ← 姿态误差常数 (1×1)
        c_k = [current_xy - ref_xy_k;
                  psi0    - ref_psi_k];

        % 选择矩阵 I_k (2.3): 3×3(K+1), 前 k 块为 I (含 v_0, 不含 v_K), 后 K+1-k 块为 0
        I_k = zeros(3, n_var);
        for l = 0:k-1
            I_k(:, l*3+1 : (l+1)*3) = eye(3);   % 从 v_0 块开始
        end

        % B·I_k (2.4): 直接在需要处用 B*I_k, 不引入中间变量
        % v_K 列恒为 0 (因 I_k 的 v_K 列为 0)

        % 三部分分离 (2.5)
        const_track = const_track + c_k' * Q * c_k;                       % (A) 常数
        g_track     = g_track     + 2 * (B*I_k)' * Q * c_k;           % (B) 一次
        H_track     = H_track     + 2 * (B*I_k)' * Q * (B*I_k);  % (C) 二次
    end

    %% =====================================================
    % 2. 控制正则化 u'Ru (2.9): H_u, g_u, const_u
    % u_k = D_k·v (v_0 在 v 中, 无常数偏移)
    % k=1..K (含 u_K = v_K - v_{K-1}, 与 V1 一致)
    % =====================================================
    R = w_control * eye(3);

    H_u = zeros(n_var, n_var);
    g_u = zeros(n_var, 1);      % 全为 0 (无常数偏移)
    const_u = 0;                % 全为 0 (无常数偏移)

    for k = 1:K
        % 差分矩阵 D_k (2.9): 3×3(K+1)
        D_k = zeros(3, n_var);
        D_k(:, k*3+1 : (k+1)*3) =  eye(3);   % v_k     块 = +I
        D_k(:, (k-1)*3+1 : k*3) = -eye(3);   % v_{k-1} 块 = -I

        % 只有二次项 (u_k 完全由 D_k·v 表达)
        H_u = H_u + 2 * D_k' * R * D_k;    % (C) 二次
    end

    %% =====================================================
    % 3. RSS 强凸正则化 ρ·‖u-û‖² (2.10): H_rho, g_rho, const_rho
    % u_k - û_k = D_k·v - û_k (û_k 是上次迭代解, 已知常数)
    % k=1..K (含 u_K - û_K, 与 V1 一致)
    % =====================================================
    H_rho = zeros(n_var, n_var);
    g_rho = zeros(n_var, 1);
    const_rho = 0;

    for k = 1:K
        % D_k 同 2.9
        D_k = zeros(3, n_var);
        D_k(:, k*3+1 : (k+1)*3) =  eye(3);   % v_k     块 = +I
        D_k(:, (k-1)*3+1 : k*3) = -eye(3);   % v_{k-1} 块 = -I

        % 三部分 (u_k - û_k = D_k·v - û_k)
        H_rho     = H_rho     + 2 * rho * D_k' * D_k;                    % (C) 二次
        g_rho     = g_rho     - 2 * rho * D_k' * u_hat(:, k);            % (B) 一次
        const_rho = const_rho +       rho * u_hat(:, k)' * u_hat(:, k);  % (A) 常数
    end

    %% =====================================================
    % 4. 合并 H, g, const
    % =====================================================
    H = H_track + H_u + H_rho;
    g = g_track + g_u + g_rho;
    objective_constant = const_track + const_u + const_rho;

    H = 0.5 * (H + H');  % 对称化

    %% =====================================================
    % 5. 等式约束: v_0 = v0 (3 条, 锁定 v_0)
    % =====================================================
    A_eq = zeros(3, n_var);
    A_eq(:, 1:3) = eye(3);          % v_0 块 = I
    b_eq = v0;
    n_eq = 3;

    %% =====================================================
    % 6. 二次约束
    % =====================================================
    Hq_list = {};  gq_list = {};  uq_list = [];

    % 6.1 轮速 SOC 约束 (论文 20b): ||H_n·v_k|| <= vimax
    %     v_k 对应 v 的第 k 块 (k=1..K, 含 v_K, v_0 不参与)
    for k = 1:K
        for n = 1:num_wheels
            Hq_k = zeros(n_var, n_var);
            gq_k = zeros(n_var, 1);
            v_k_start = k*3 + 1;          % v 的第 k 块 (v_0 是第 0 块, 从 1 开始)
            v_k_end   = (k+1)*3;
            Hq_k(v_k_start:v_k_end, v_k_start:v_k_end) = 2 * (Hn{n}' * Hn{n});
            Hq_list{end+1} = Hq_k;
            gq_list{end+1} = gq_k;
            uq_list(end+1) = vimax^2;
        end
    end

    % 6.2 转向锥凸化约束 (论文 20a)
    % TODO: V2 中 u = D·v (d=0), 需将 V1 的 u 依赖替换为 v 的差分
    %       约束结构稍变但数学等价, 待后续实现

    n_qcqp = length(uq_list);

    %% =====================================================
    % 7. 返回 (字段对齐 HPIPM API)
    % =====================================================
    qp_problem.H = H;
    qp_problem.g = g;
    qp_problem.A = A_eq;
    qp_problem.b = b_eq;
    qp_problem.C = [];
    qp_problem.d = [];
    qp_problem.lb = [];
    qp_problem.ub = [];
    qp_problem.Hq = Hq_list;
    qp_problem.gq = gq_list;
    qp_problem.uq = uq_list;
    qp_problem.n_var = n_var;       % 21 (V1 是 36)
    qp_problem.K = K;
    qp_problem.n_eq = n_eq;         % 3  (V1 是 18)
    qp_problem.n_qcqp = n_qcqp;
    qp_problem.objective_constant = objective_constant;
end
```

### 7.1 V1 vs V2 关键差异

| 项 | V1 | V2 |
|---|---|---|
| 决策变量 | $\boldsymbol{x} = [\mathbf{u}(18); \nu(18)] = 36$ 维 | $\boldsymbol{v} = [v_0; v_1..v_K] = 21$ 维 |
| 等式约束 | 18 条 (动力学递推) | 3 条 (仅 $v_0 = v_0^{\text{cur}}$) |
| $\mathbf{H}$ 填充方式 | 逐 $k$ 双重循环填交叉项 | 矩阵乘法 $(\mathbf{B}\,\mathbf{I}_k)^\top \mathbf{Q} (\mathbf{B}\,\mathbf{I}_k)$ |
| $\mathbf{u}^\top \mathbf{R} \mathbf{u}$ 的 $\mathbf{H}$ | $\mathbf{u}$ 块纯对角 | $\boldsymbol{v}$ 块三对角 (差分消元), $\mathbf{g}_u = \mathbf{0}$, $\text{const}_u = 0$ |
| 转向锥约束 | 直接用 $\mathbf{u}$ | TODO: 需重写为 $\boldsymbol{v}$ 差分 |
| 等价性 | — | 与 V1 严格等价 (经约束消元后 H/g/const 完全一致) |
