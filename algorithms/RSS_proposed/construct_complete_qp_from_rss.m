function qp_problem = construct_complete_qp_from_rss(path, step, current_nu, state, u_hat, params)
% CONSTRUCT_COMPLETE_QP_FROM_RSS
% 将论文 RSS26 的凸子问题 Q_K(û) (公式 17) 构造为 dense QCQP 标准形式。
%
% 论文: RSS26 "Exploit Agile Mobility of Steerable-Wheeled Mobile Robots:
%        A Fast Motion Planning Approach"
%
% 论文中的凸子问题 Q_K(û) (公式 17):
%   min  f(u, û) = g(u) + ρ·||u - û||^2
%   s.t. C^k_{i,n}(u, û) ≤ 0     (凸化后的转向锥约束, 公式 15-16)
%        u ∈ ∩_j U_j              (凸约束集, 这里取轮速约束 (20b))
%
% 代价函数 g(u) (论文公式 18):
%   g(u) = Σ_{k=1}^K [ e_k^T Q e_k + u_k^T R u_k ]
%   其中 Q = diag(30, 30, 1), R = diag(0.3, 0.3, 0.3)
%   e_k = ξ_w(t0+kτ) - ξ^k_r 为跟踪误差
%
% 跟踪误差 e_k 的一阶展开 (论文公式 19, SO(2) 上的一阶展开):
%   e_k = ξ_w(t0) + [R(ψ_w(t0)), 0; 0, 1]·Σ_{l=0}^{k-1} ν_l·τ - ξ^k_r + O(τ^2)
%
% 其中 O(τ^2) 是 SO(2) 一阶展开的高阶余项 (τ=0.01s 时 O(τ^2)=1e-4, 可忽略)
% 代码中略去 O(τ^2), 即:
%   位置误差 = current_xy - ref_xy + R(ψ0)·Σ ν_l·τ (旋转矩阵线性化)
%   姿态误差 = ψ0 + current_nu(3)·τ + Σ ν_3·τ - ref_psi
%
% 约束:
%   (1) 动力学等式 (论文公式 20c): ν_{k+1} = ν_k + u_{k+1}
%   (2) 轮速 SOC 约束 (论文公式 20b): ||H_n·ν_k|| ≤ vimax (凸, 直接二次形式)
%   (3) 转向锥凸化约束 (论文公式 15-16, Prop.1):
%       C^k_{i,n}(u, û) = A^k_{i,n}(u) - B^k_{i,n}(û) - L^k_{i,n}(u, û) ≤ 0
%       其中 A, B, L 定义见公式 (16) 及 Appendix A
%
% dense QCQP 标准形式 (HPIPM 论文 Section 2.1):
%
% [HPIPM 论文公式 (1) — 完整 dense QP (线性约束, 含 slack)]
%   min_{v,s}  1/2 [v;1]^T [H g; g^T 0] [v;1]
%             + 1/2 [s^l;s^u;1]^T [Z^l 0 z^l; 0 Z^u z^u; (z^l)^T (z^u)^T 0] [s^l;s^u;1]
%   s.t. A v = b                                                          (等式)
%        [v_; d_] <= [J^{b,v}; C] v + [J^{s,v}; J^{s,g}] s^l              (下界+slack)
%        [J^{b,v}; C] v - [J^{s,v}; J^{s,g}] s^u <= [v^; d^]             (上界+slack)
%        s^l >= s^l_lb,  s^u >= s^u_lb                                     (slack 非负)
%
% [dense QCQP 扩展 — 在 dense QP 基础上增加二次约束]
%   0.5 v^T Hq_i v + gq_i^T v <= uq_i    (二次不等式, 亦可带 slack)
%
% [本代码使用硬约束子集 (nb=0, ng=0, ns=0, 无 slack)]
%   min  0.5 x^T H x + g^T x
%   s.t. A x = b                                  (动力学等式, 论文 (20c))
%        0.5 x^T Hq_i x + gq_i^T x <= uq_i        (二次不等式: 轮速 (20b) + 转向锥 (20a))
% 即 HPIPM 的 slack/box/一般线性约束均未启用, 退化为纯 QCQP (硬约束)
%
% 变量排列: x = [u(:); nu(:)]  (36维, K=6)
%   u(i,k) 的全局索引 = (k-1)*3 + i         (i=1,2,3; k=1,...,K)
%   nu(i,k) 的全局索引 = nu_start + (k-1)*3 + (i-1)
%   其中 nu_start = 3*K + 1 = 19 (K=6时)

    %% =====================================================
    % 参数提取 (论文 IV-A 实验设置)
    % ======================================================
    % 论文 IV-A: "prediction horizon K=6", "discretization interval 0.01s"
    K = 6;
    dt = params.dt;              % τ: 离散化间隔 (论文 IV-A: 0.01s)
    phidotmax = params.phidotmax;  % ω_max: 最大转向角速率 (论文 (5): 5π rad/s)
    vimax = params.vimax;        % z_max: 最大轮速 (论文 (20b): 5 m/s)
    wheel_pos = params.wheel_pos;  % d_n: 轮子位置 (论文 (3): [dx, dy])
    num_wheels = size(wheel_pos, 1);  % N: 轮数 (论文: 4)

    current_xy = [state(1); state(2)];  % 当前世界系位置 ξ_w(t0) 的 xy 分量
    psi0 = state(3);                     % 当前航向 ψ_w(t0)
    % 论文 (1) 式的旋转矩阵 R(ψ_w) (2D 版本, 用于位置误差展开)
    R_psi0 = [cos(psi0), -sin(psi0); sin(psi0), cos(psi0)];

    % 代价函数权重 (论文公式 18: Q=diag(30,30,1), R=diag(0.3,0.3,0.3))
    w_pos = 30;      % Q 的位置分量 (论文 Q=diag(30,30,1))
    w_psi = 1;       % Q 的姿态分量
    w_control = 0.3; % R 的对角元 (论文 R=diag(0.3,0.3,0.3))
    rho = 0.01;       % ρ: 强凸正则化参数 (论文 (17) 式)

    % =====================================================
    % 论文公式 (3): 特征矩阵 H_n (2×3)
    %   H_n = [1, 0, -dy_n; 0, 1, dx_n]
    % 论文公式 (4): z_n = H_n · ν̇_c (alignment 条件)
    % =====================================================
    H_cell = cell(1, num_wheels);
    for n = 1:num_wheels
        H_cell{n} = [
            1, 0, -wheel_pos(n, 2);   % -dy_n
            0, 1,  wheel_pos(n, 1)    %  dx_n
        ];
    end

    % M_n = H_n^T · H_n (3×3), 用于轮速二次约束 (||H_n·ν||^2 = ν^T·M_n·ν)
    M_cell = cell(1, num_wheels);
    for n = 1:num_wheels
        M_cell{n} = H_cell{n}' * H_cell{n};
    end

    % 论文公式 (11): δ_ω = ω_max · τ (单步最大转向角变化)
    % 论文 IV-A: "steering rate 5π rad/s, dt=0.01s → δ_ω = 5π·0.01 = π/20"
    delta_theta = dt * phidotmax;

    %% =====================================================
    % 决策变量排列：x = [u(:); nu(:)]
    % u: 控制增量 (论文 (9): ν_{k+1}=ν_k+u_{k+1})
    % nu: 车体系速度序列 ν_k (论文 (8): ν_k = ξ̇_c(t_k))
    % ======================================================

    n_var = 6 * K;            % 36维: u(18) + nu(18)
    nu_start = 3*K + 1;       % nu(1,1) 在 x 中的起始索引 = 19

    %% =====================================================
    % 1. 代价函数: 0.5·x^T·H·x + g^T·x
    %    对应论文公式 (17) 的 f(u,û) = g(u) + ρ·||u-û||^2
    %    其中 g(u) = Σ e_k^T·Q·e_k + u_k^T·R·u_k (论文公式 18)
    % ======================================================
    % 论文公式 (19): e_k 的一阶展开 (位置部分)
    %   position_error_k = current_xy - ref_xy_k + R(ψ0)·Σ_{l=0}^{k-1} ν_l·τ
    % 其中 ν_0 = current_nu, ν_k = ν_{k-1} + u_k (论文 (9))
    % => Σ_{l=0}^{k-1} ν_l·τ = current_nu·τ + Σ_{j=1}^{k-1} nu(1:2,j)·τ
    %
    % ||v||^2 = v^T·v, 在 0.5·x^T·H·x 形式中 Hessian 系数需乘 2

    % 命名约定: _mat=矩阵, _vec=向量, _eq=等式约束, _ineq=不等式约束
    % 这些中间变量在最后统一赋值到 qp_problem 结构体 (字段名对齐 HPIPM API)
    H_mat = zeros(n_var, n_var);  % 目标 Hessian (对应 HPIPM H)
    g_vec = zeros(n_var, 1);      % 目标线性项 (对应 HPIPM g)

    % -------------------------------------------------------
    % 1.1 位置跟踪代价: Σ_{k=2}^K w_pos · ||position_error_k||^2
    %     论文公式 (18) 位置部分 + 公式 (19) 展开
    % -------------------------------------------------------
    % position_error_k = c_k + R_psi0·dt·S_k
    %   c_k = current_xy - ref_xy + R_psi0·current_nu(1:2)·dt  (常数)
    %   S_k = Σ_{j=1}^{k-1} nu(1:2,j)                           (决策变量)
    %
    % 代价 = w_pos·(||c_k||^2 + 2·c_k^T·R_psi0·dt·S_k + dt^2·||S_k||^2)
    %
    % 因 R_psi0^T·R_psi0 = I (旋转矩阵正交性):
    %   ||R_psi0·dt·S_k||^2 = dt^2·||S_k||^2
    %
    % 在 0.5·x^T·H·x 形式:
    %   二次项 H: H(nu_x(i), nu_x(j)) += 2·w_pos·dt^2, 同理 nu_y
    %             (nu_x 与 nu_y 无交叉, 因 ||S||^2 = Σx^2 + Σy^2)
    %   一次项 g: g(nu_x(j)) += 2·w_pos·dt·(R_psi0^T·c_k)(1)

    for k = 2:K
        % 参考位置 ξ^k_r (论文 IV-A: 贝塞尔曲线采样点)
        ref_idx = min(size(path, 2), step + k);
        ref_xy = path(1:2, ref_idx);

        % c_k = current_xy - ref_xy + R_psi0·current_nu(1:2)·dt
        % (常数部分含 current_nu·dt 位移, 来自 ν_0·τ 的贡献)
        c_k = current_xy - ref_xy + R_psi0 * current_nu(1:2) * dt;

        % 梯度方向 R_psi0^T · c_k (2×1)
        grad_dir = R_psi0' * c_k;

        % 二次项: dt^2·(Σ_{j=1}^{k-1} nu_x(j))^2 = dt^2·Σ_{i,j} nu_x(i)·nu_x(j)
        % H(nu_x(i), nu_x(j)) += 2·w_pos·dt^2 (对所有 i,j ∈ [1, k-1])
        for i = 1:k-1
            for j = 1:k-1
                nu_x_i = nu_start + (i-1)*3;       % nu(1,i) 全局索引
                nu_x_j = nu_start + (j-1)*3;
                nu_y_i = nu_start + (i-1)*3 + 1;   % nu(2,i) 全局索引
                nu_y_j = nu_start + (j-1)*3 + 1;

                H_mat(nu_x_i, nu_x_j) = H_mat(nu_x_i, nu_x_j) + 2 * w_pos * dt^2;
                H_mat(nu_y_i, nu_y_j) = H_mat(nu_y_i, nu_y_j) + 2 * w_pos * dt^2;
            end
        end

        % 一次项: 2·w_pos·dt·(R_psi0^T·c_k)^T · S_k
        for j = 1:k-1
            nu_x_j = nu_start + (j-1)*3;
            nu_y_j = nu_start + (j-1)*3 + 1;

            g_vec(nu_x_j) = g_vec(nu_x_j) + 2 * w_pos * dt * grad_dir(1);
            g_vec(nu_y_j) = g_vec(nu_y_j) + 2 * w_pos * dt * grad_dir(2);
        end
    end

    % -------------------------------------------------------
    % 1.2 姿态跟踪代价: Σ_{k=1}^K w_psi · (psi_k - ref_psi_k)^2
    %     论文公式 (18) 姿态部分
    % -------------------------------------------------------
    % psi(k) = ψ0 + current_nu(3)·dt + dt·Σ_{j=1}^{k-1} nu(3,j)
    % 姿态误差 = psi_c_k + dt·Σ_{j=1}^{k-1} nu(3,j)
    %   其中 psi_c_k = ψ0 + current_nu(3)·dt - ref_psi_k  (常数)
    %
    % 在 0.5·x^T·H·x 形式:
    %   二次项: H(nu_psi(i), nu_psi(j)) += 2·w_psi·dt^2
    %   一次项: g(nu_psi(j)) += 2·w_psi·dt·psi_c_k

    for k = 1:K
        ref_idx = min(size(path, 2), step + k);
        ref_psi = path(3, ref_idx);  % 参考航向

        % psi_c_k = ψ0 + current_nu(3)·dt - ref_psi_k (常数部分)
        psi_c_k = psi0 + current_nu(3)*dt - ref_psi;

        % 二次项: dt^2·(Σ nu_psi)^2 → 跨阶段交叉项
        for i = 1:k-1
            for j = 1:k-1
                nu_psi_i = nu_start + (i-1)*3 + 2;  % nu(3,i) 全局索引
                nu_psi_j = nu_start + (j-1)*3 + 2;

                H_mat(nu_psi_i, nu_psi_j) = H_mat(nu_psi_i, nu_psi_j) + 2 * w_psi * dt^2;
            end
        end

        % 一次项
        for j = 1:k-1
            nu_psi_j = nu_start + (j-1)*3 + 2;

            g_vec(nu_psi_j) = g_vec(nu_psi_j) + 2 * w_psi * dt * psi_c_k;
        end
    end

    % -------------------------------------------------------
    % 1.3 控制正则化: w_control · Σ ||u_k||^2
    %     论文公式 (18) 第二项: u_k^T·R·u_k, R=0.3·I
    % -------------------------------------------------------
    % ||u||^2 在 0.5·x^T·H·x 形式: H 对角块 = 2·w_control·I
    for k = 1:K
        for i = 1:3
            u_idx = (k-1)*3 + i;
            H_mat(u_idx, u_idx) = H_mat(u_idx, u_idx) + 2 * w_control;
        end
    end

    % -------------------------------------------------------
    % 1.4 RSS 强凸正则化: ρ · ||u - û||^2
    %     论文公式 (17): f(u,û) = g(u) + ρ·||u-û||^2
    %     提供强凸性, 保证 S(û) 唯一解 (论文 Theorem 1)
    % -------------------------------------------------------
    % ρ·||u-û||^2 = ρ·(||u||^2 - 2·u^T·û + ||û||^2)
    % H += 2·ρ·I (在 u 块), g -= 2·ρ·û
    for k = 1:K
        for i = 1:3
            u_idx = (k-1)*3 + i;
            H_mat(u_idx, u_idx) = H_mat(u_idx, u_idx) + 2 * rho;
            g_vec(u_idx) = g_vec(u_idx) - 2 * rho * u_hat(i, k);
        end
    end

    % 确保 H 对称 (数值稳定性)
    H_mat = 0.5 * (H_mat + H_mat');

    %% =====================================================
    % 1.5 目标函数常数项 (不影响最优解, 仅用于 obj_value 比较)
    % ======================================================
    % J = 0.5·x^T·H·x + g^T·x + c
    % c = Σ w_pos·||c_k||^2 + Σ w_psi·psi_c_k^2 + ρ·||û||^2

    objective_constant = 0;

    % (a) 位置跟踪常数
    for k = 2:K
        ref_idx = min(size(path, 2), step + k);
        ref_xy = path(1:2, ref_idx);
        c_k = current_xy - ref_xy + R_psi0 * current_nu(1:2) * dt;
        objective_constant = objective_constant + w_pos * (c_k' * c_k);
    end

    % (b) 姿态跟踪常数
    for k = 1:K
        ref_idx = min(size(path, 2), step + k);
        ref_psi = path(3, ref_idx);
        psi_c_k = psi0 + current_nu(3)*dt - ref_psi;
        objective_constant = objective_constant + w_psi * psi_c_k^2;
    end

    % (c) RSS 正则常数: ρ·||û||^2
    objective_constant = objective_constant + rho * sum(u_hat(:).^2);

    %% =====================================================
    % 2. 等式约束 (论文公式 20c: 动力学递推)
    % ======================================================
    % 论文 (9) 式: ν_{k+1} = ν_k + u_{k+1}
    % 离散化为: nu(:,1) = current_nu + u(:,1)
    %           nu(:,k+1) = nu(:,k) + u(:,k+1)
    % => A·x = b 形式

    n_eq = 3*K;  % 3 (初始) + 3·(K-1) (递推) = 18
    A_eq = zeros(n_eq, n_var);  % 等式约束矩阵 (对应 HPIPM A)
    b_eq = zeros(n_eq, 1);     % 等式约束右端 (对应 HPIPM b)

    eq_row = 1;

    % 2.1 初始条件: nu(:,1) - u(:,1) = current_nu
    %     即 ν_1 = ν_0 + u_1 (论文 (9), k=0)
    for i = 1:3
        u_idx = i;                         % u(i,1) 全局索引
        nu_idx = nu_start + i - 1;          % nu(i,1) 全局索引

        A_eq(eq_row, u_idx) = -1;          % -u(i,1)
        A_eq(eq_row, nu_idx) = 1;           % +nu(i,1)
        b_eq(eq_row) = current_nu(i);      % = current_nu(i)

        eq_row = eq_row + 1;
    end

    % 2.2 递推约束: nu(:,k+1) - nu(:,k) - u(:,k+1) = 0
    %     即 ν_{k+1} = ν_k + u_{k+1} (论文 (9), k=1..K-1)
    for k = 1:K-1
        for i = 1:3
            u_idx = k*3 + i;                        % u(i,k+1)
            nu_k_idx = nu_start + (k-1)*3 + i - 1;  % nu(i,k)
            nu_kp1_idx = nu_start + k*3 + i - 1;    % nu(i,k+1)

            A_eq(eq_row, u_idx) = -1;       % -u(i,k+1)
            A_eq(eq_row, nu_k_idx) = -1;    % -nu(i,k)
            A_eq(eq_row, nu_kp1_idx) = 1;   % +nu(i,k+1)
            b_eq(eq_row) = 0;               % = 0

            eq_row = eq_row + 1;
        end
    end

    %% =====================================================
    % 3. 线性不等式约束 (暂无, 所有非线性约束在二次约束中)
    % ======================================================

    % C_ineq, d_ineq: 一般线性不等式约束 d_ <= Cx <= d^ (对应 HPIPM C, d)
    % 本代码未启用 (设 ng=0), 仅保留字段以对齐 HPIPM API
    C_ineq = [];
    d_ineq = [];

    %% =====================================================
    % 4. 二次不等式约束 (QCQP)
    %    对应论文公式 (20a) 转向锥 + (20b) 轮速
    % ======================================================

    Hq_list = {};
    gq_list = {};
    uq_list = [];

    % 构造 û 对应的 ν̂ 序列 (用于转向锥凸化的 B 项, 论文 (16))
    % ν̂_k = current_nu + Σ_{j=1}^k û_j (论文 Appendix A: {ν̂_k} defined by ν̂_{k+1}=ν̂_k+û_{k+1})
    nu_hat = zeros(3, K);
    nu_hat(:, 1) = current_nu + u_hat(:, 1);
    for k = 1:K-1
        nu_hat(:, k+1) = nu_hat(:, k) + u_hat(:, k+1);
    end

    % -------------------------------------------------------
    % 4.1 轮速 SOC 约束 (论文公式 20b): ||H_n·ν_k|| ≤ vimax
    %     凸约束, 直接转二次形式: ν_k^T·M_n·ν_k ≤ vimax^2
    %     QCQP 形式: 0.5·x^T·Hq·x + gq^T·x ≤ uq
    %       Hq 在 nu(:,k) 块 = 2·M_n, gq = 0, uq = vimax^2
    %     共 N×K = 4×6 = 24 条约束
    % -------------------------------------------------------

    for k = 1:K
        for n = 1:num_wheels
            Hq_k = zeros(n_var, n_var);
            gq_k = zeros(n_var, 1);

            % nu(:,k) 的索引范围
            nu_k_start = nu_start + (k-1)*3;
            nu_k_end = nu_start + (k-1)*3 + 2;

            % ||H_n·ν_k||^2 = ν_k^T·(H_n^T·H_n)·ν_k = ν_k^T·M_n·ν_k
            % 在 0.5·x^T·Hq·x 形式: Hq = 2·M_n (在 nu_k 块)
            Hq_k(nu_k_start:nu_k_end, nu_k_start:nu_k_end) = 2 * M_cell{n};

            Hq_list{end+1} = Hq_k;
            gq_list{end+1} = gq_k;
            uq_list(end+1) = vimax^2;  % 右端 = z_max^2 (论文: 5^2=25)
        end
    end

    % -------------------------------------------------------
    % 4.2 第一组转向锥凸化约束 (论文公式 12, R_1 = R(π/2-δ_ω))
    % -------------------------------------------------------
    % 论文公式 (12): (R_i·z_n(t_k))^T·z_n(t_{k+1}) ≥ 0
    %   R_1 = R(π/2 - δ_ω), R_2 = R_1^T
    % 代入 z_n = H_n·ν_k (论文 (4)) 得论文公式 (13):
    %   ν_{k-1}^T·H_n^T·R_i^T·H_n·(ν_{k-1} + u_k) ≥ 0  (非凸双线性)
    %
    % 经 Proposition 1 凸化后 (论文公式 15-16):
    %   C^k_{i,n}(u, û) = A^k_{i,n}(u) - B^k_{i,n}(û) - L^k_{i,n}(u, û) ≤ 0
    %   A = 1/2·||H_n·ν_{k-1}||^2 + 1/2·||H_n·(ν_{k-1}+u_k)||^2  (凸, 关于 u)
    %   B = 1/2·||(I+R_i)·H_n·ν̂_{k-1} + R_i·H_n·û_k||^2          (常数, 在 û 处)
    %   L = ∇_u B|_û · (u - û)                                    (B 的一阶展开, 线性)

    % R_1 = R(π/2 - δ_ω) (论文 (12) 式)
    R1 = [
        sin(delta_theta), -cos(delta_theta);
        cos(delta_theta),  sin(delta_theta)
    ];

    % 构造 N×K = 24 条凸化转向锥约束 (R_1 组)
    result = add_steering_cone_constraints( ...
        Hq_list, gq_list, uq_list, ...
        R1, K, num_wheels, H_cell, M_cell, ...
        u_hat, nu_hat, current_nu, ...
        n_var, nu_start, delta_theta ...
    );

    Hq_list = result{1};
    gq_list = result{2};
    uq_list = result{3};

    % -------------------------------------------------------
    % 4.3 第二组转向锥凸化约束 (论文公式 12, R_2 = R_1^T)
    % -------------------------------------------------------
    % R_2 = R_1^T (论文 (12): "R_2 = R_1^T")
    R2 = [
        sin(delta_theta),  cos(delta_theta);
       -cos(delta_theta),  sin(delta_theta)
    ];

    % 构造 N×K = 24 条凸化转向锥约束 (R_2 组)
    result = add_steering_cone_constraints( ...
        Hq_list, gq_list, uq_list, ...
        R2, K, num_wheels, H_cell, M_cell, ...
        u_hat, nu_hat, current_nu, ...
        n_var, nu_start, delta_theta ...
    );

    Hq_list = result{1};
    gq_list = result{2};
    uq_list = result{3};

    %% =====================================================
    % 5. 返回 QCQP 问题结构 (供 HPIPM 求解)
    % ======================================================

    qp_problem.H = H_mat;       % 代价 Hessian
    qp_problem.g = g_vec;       % 代价线性项
    qp_problem.A = A_eq;        % 等式约束矩阵 (动力学)
    qp_problem.b = b_eq;        % 等式约束右端
    qp_problem.C = C_ineq;      % 一般线性不等式矩阵 (空, ng=0)
    qp_problem.d = d_ineq;      % 一般线性不等式右端 (空)
    qp_problem.lb = [];         % box 下界 (空, nb=0)
    qp_problem.ub = [];         % box 上界 (空, nb=0)

    % 二次约束 (轮速 24 + 转向锥 48 = 72 条)
    qp_problem.Hq = Hq_list;
    qp_problem.gq = gq_list;
    qp_problem.uq = uq_list;

    % 元数据
    qp_problem.n_var = n_var;           % 36
    qp_problem.K = K;                  % 6
    qp_problem.n_eq = n_eq;            % 18
    qp_problem.n_qcqp = length(uq_list);  % 72
    qp_problem.objective_constant = objective_constant;

end


%% =========================================================
% 构造转向锥凸化约束 (论文 Proposition 1, 公式 15-16)
% ==========================================================
% 论文公式 (15): C^k_{i,n}(u, û) = A - B - L ≤ 0
%
% 论文公式 (16):
%   A^k_{i,n}(u) = 1/2·||H_n·ν_{k-1}||^2 + 1/2·||H_n·(ν_{k-1}+u_k)||^2
%   B^k_{i,n}(û) = 1/2·||(I+R_i)·H_n·ν̂_{k-1} + R_i·H_n·û_k||^2
%   L^k_{i,n}(u, û) = tr((∇_u B|_û)^T·(u - û))  (B 的一阶展开)
%
% 论文 Appendix A 公式 (24-25): L 的分块计算
%   L = Σ_{l=1}^{k-1} (∇_{u_l} B|_û)^T·(u_l - û_l) + (∇_{u_k} B|_û)^T·(u_k - û_k)
%   其中 ∇_{u_l} b = (I+R_i)·H_n  (l < k),  ∇_{u_k} b = R_i·H_n  (l = k)
%   b = (I+R_i)·H_n·ν̂_{k-1} + R_i·H_n·û_k
%
% 展开后的二次形式 (0.5·x^T·Hq·x + gq^T·x ≤ uq):
%   k > 1 时:
%     A 项 (二次): ||H_n·ν_{k-1}||^2 + ||H_n·(ν_{k-1}+u_k)||^2
%       = 2·ν_{k-1}^T·M_n·ν_{k-1} + 2·ν_{k-1}^T·M_n·u_k + u_k^T·M_n·u_k
%       → Hq(ν_{k-1}, ν_{k-1}) = 4·M_n, Hq(ν_{k-1}, u_k) = 2·M_n, Hq(u_k, u_k) = 2·M_n
%     L 项 (线性): -Σ_{l<k} 2·b^T·(I+R_i)·H_n·(u_l - û_l) - 2·b^T·R_i·H_n·(u_k - û_k)
%     B 项 (常数): ||b||^2 (移到 uq)
% ==========================================================

function result = add_steering_cone_constraints( ...
    Hq_list, gq_list, uq_list, ...
    R, K, num_wheels, H_cell, M_cell, ...
    u_hat, nu_hat, current_nu, ...
    n_var, nu_start, delta_theta ...
)

    for k = 1:K
        for n = 1:num_wheels
            Hq_k = zeros(n_var, n_var);
            gq_k = zeros(n_var, 1);

            Hn = H_cell{n};  % 论文 (3): 2×3 特征矩阵
            Mn = M_cell{n};  % M_n = H_n^T·H_n

            if k > 1
                % =====================================================
                % k > 1: 标准情况 (ν_{k-1} 是决策变量)
                % =====================================================
                % b = (I+R)·H_n·ν̂_{k-1} + R·H_n·û_k (论文 (16) B 项的核)
                %   (ν̂_{k-1} 来自上一次迭代的 û 序列)
                lv = (eye(2) + R) * Hn * nu_hat(:, k-1) + R * Hn * u_hat(:, k);

                % ----- A 项 (二次部分, 论文 (16)) -----
                % ||H_n·ν_{k-1}||^2 + ||H_n·(ν_{k-1}+u_k)||^2
                % = ν_{k-1}^T·M_n·ν_{k-1} + (ν_{k-1}+u_k)^T·M_n·(ν_{k-1}+u_k)
                % 在 0.5·x^T·Hq·x 形式:
                %   Hq(ν_{k-1}, ν_{k-1}) = 4·M_n  (来自两项各贡献 2·M_n)
                %   Hq(ν_{k-1}, u_k) = 2·M_n      (交叉项)
                %   Hq(u_k, u_k) = 2·M_n

                nu_k1_start = nu_start + (k-2)*3;    % ν_{k-1} 块起始
                nu_k1_end = nu_start + (k-2)*3 + 2;
                u_k_start = (k-1)*3 + 1;             % u_k 块起始
                u_k_end = k*3;

                Hq_k(nu_k1_start:nu_k1_end, nu_k1_start:nu_k1_end) = 4 * Mn;
                Hq_k(nu_k1_start:nu_k1_end, u_k_start:u_k_end) = 2 * Mn;
                Hq_k(u_k_start:u_k_end, nu_k1_start:nu_k1_end) = 2 * Mn;
                Hq_k(u_k_start:u_k_end, u_k_start:u_k_end) = 2 * Mn;

                % ----- L 项 (线性部分, 论文 Appendix A (24)-(25)) -----
                % L = Σ_{l=1}^{k-1} b^T·(I+R)·H_n·(u_l - û_l) + b^T·R·H_n·(u_k - û_k)
                % 在约束 C = A - B - L ≤ 0 中, L 移到左边变 -L:
                %   对 u_l (l<k): gq += -2·b^T·(I+R)·H_n  (系数 -2 因 0.5 形式)
                %   对 u_k:       gq += -2·b^T·R·H_n

                for l = 1:k-1
                    % l < k: ∇_{u_l} b = (I+R)·H_n (论文 (25))
                    coeff = -2 * ((eye(2) + R)' * lv)' * Hn;  % 1×3
                    u_l_start = (l-1)*3 + 1;
                    gq_k(u_l_start:u_l_start+2) = gq_k(u_l_start:u_l_start+2) + coeff';
                end
                % l = k: ∇_{u_k} b = R·H_n (论文 (25))
                coeff_k = -2 * (R' * lv)' * Hn;  % 1×3
                gq_k(u_k_start:u_k_start+2) = gq_k(u_k_start:u_k_start+2) + coeff_k';

                % ----- B 项 + L 中的 û 常数 (移到右端 uq) -----
                % uq = ||b||^2 - Σ 2·b^T·(I+R)·H_n·û_l - 2·b^T·R·H_n·û_k
                u_hat_const = 0;
                for l = 1:k-1
                    u_hat_const = u_hat_const + 2 * lv' * (eye(2) + R) * Hn * u_hat(:, l);
                end
                u_hat_const = u_hat_const + 2 * lv' * R * Hn * u_hat(:, k);

                uq_val = lv' * lv - u_hat_const;  % ||b||^2 - L 中的 û 常数部分

            else
                % =====================================================
                % k = 1: 边界情况 (ν_0 = current_nu 是已知常数, 非决策变量)
                % =====================================================
                % b = (I+R)·H_n·ν_0 + R·H_n·û_1 (ν_0 替代 ν̂_0)
                lv = (eye(2) + R) * Hn * current_nu + R * Hn * u_hat(:, 1);

                % ----- A 项 (二次部分) -----
                % ||H_n·ν_0||^2 (常数, 移到 uq) + ||H_n·(ν_0+u_1)||^2
                % 后者 = ν_0^T·M_n·ν_0 + 2·ν_0^T·M_n·u_1 + u_1^T·M_n·u_1
                % 决策变量部分: u_1^T·M_n·u_1 → Hq(u_1, u_1) = 2·M_n

                u_1_start = 1;
                u_1_end = 3;

                Hq_k(u_1_start:u_1_end, u_1_start:u_1_end) = 2 * Mn;

                % ----- 线性部分 -----
                % 一次项: 2·ν_0^T·M_n·u_1 (来自 A) - 2·b^T·R·H_n·u_1 (来自 -L)
                gq_k(u_1_start:u_1_end) = 2 * Mn' * current_nu - 2 * Hn' * R' * lv;

                % ----- 常数项 (移到 uq) -----
                % uq = -2·ν_0^T·M_n·ν_0 (A 中常数) + ||b||^2 (B) - 2·b^T·R·H_n·û_1 (L 中 û 常数)
                uq_val = -2 * current_nu' * Mn * current_nu + lv' * lv - 2 * lv' * R * Hn * u_hat(:, 1);
            end

            % 确保 Hq 对称 (数值稳定性)
            Hq_k = 0.5 * (Hq_k + Hq_k');

            Hq_list{end+1} = Hq_k;
            gq_list{end+1} = gq_k;
            uq_list(end+1) = uq_val;
        end
    end

    result = {Hq_list, gq_list, uq_list};
end
