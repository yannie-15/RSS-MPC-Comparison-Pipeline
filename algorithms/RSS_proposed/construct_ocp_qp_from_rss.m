function ocp = construct_ocp_qp_from_rss(path, step, v0, state, u_hat, params)
% CONSTRUCT_OCP_QP_FROM_RSS
% 误差状态 OCP QP 构造 (HPIPM ocp_qp Python 接口, 线性化二次约束)
%
% 状态 x_k = [e_k; v_k], e_k = xi_k - xi_k^ref (跟踪误差)
% 控制 u_n = v_{n+1} - v_n (HPIPM 编号, = 论文 u_{n+1})
% 动力学 x_{n+1} = A x_n + B u_n + b_n (含时变 bias)
%
% HPIPM 字段名: A, B, b, Q, R, r, C, D, lg, ug (一般线性约束)
% HPIPM dim 字段: nx, nu, ng, nbx
% ×2 约定: 构造时统一×2 (Q_eff, R_eff, r 均已×2), set 时不再×2; const 不×2
%
% 注: 由于 HPIPM ocp_qcqp IPM solver 存在 bug (status=3 NAN_SOL),
%     改用 ocp_qp + 线性化二次约束 (SCP 外层迭代).
%     二次约束在当前线性化点 (v_hat, u_hat) 处做一阶泰勒展开, 转为一般线性约束.
%     外层 SCP 迭代 (control_RSS.m 中 m=1..3) 保证收敛到原二次约束的解.

    %% =====================================================
    % 参数提取 (论文 IV-A 实验设置, 与 dense 版一致)
    % =====================================================
    K = 6;                                     % 预测时域 (论文 IV-A: K=6)
    dt = params.dt;                            % tau: 离散化间隔 (论文 IV-A: 0.01s)
    phidotmax = params.phidotmax;              % omega_max: 最大转向角速率 (论文 (5): 5π rad/s)
    vimax = params.vimax;                      % z_max: 最大轮速 (论文 (20b): 5 m/s)
    wheel_pos = params.wheel_pos;              % d_n: 轮子位置 (论文 (3): [dx, dy])
    num_wheels = size(wheel_pos, 1);           % N: 轮数 (论文: 4)

    current_xy = [state(1); state(2)];         % 当前世界系位置 xi_w(t0) 的 xy 分量
    psi0 = state(3);                           % 当前航向 psi_w(t0)
    R_psi0 = [cos(psi0), -sin(psi0); sin(psi0), cos(psi0)];  % 论文 (1) 旋转矩阵 R(psi_w)

    % 代价函数权重 (论文公式 18: Q=diag(30,30,1), R=diag(0.3,0.3,0.3))
    w_pos = 30;      % Q 的位置分量
    w_psi = 1;       % Q 的姿态分量
    w_control = 0.3; % R 的对角元
    rho = 0.01;      % rho: RSS 强凸正则化参数 (论文 (17))

    % 论文公式 (3): 特征矩阵 H_n (2×3), H_n = [1, 0, -dy_n; 0, 1, dx_n]
    Hn = cell(1, num_wheels);
    for n = 1:num_wheels
        Hn{n} = [1, 0, -wheel_pos(n,2); 0, 1, wheel_pos(n,1)];
    end

    % 论文公式 (11): delta_theta = omega_max * tau (单步最大转向角变化)
    delta_theta = dt * phidotmax;

    %% =====================================================
    % 1. 动力学矩阵 (常数, LTI)
    %    A (6x6): 状态转移矩阵, 右上分块 = blkdiag(R(psi0)*dt, dt) (论文公式19)
    %    B (6x3): 控制矩阵 (HPIPM 动力学字段 B), = [0; I3]
    % =====================================================
    A  = [eye(3), [R_psi0*dt, zeros(2,1); zeros(1,2), dt]; zeros(3,3), eye(3)];  % 6×6
    B  = [zeros(3,3); eye(3)];                % 6×3 控制矩阵 (HPIPM 字段 B)

    % 状态选择矩阵 C = [I3, 0] (e_k = C * x_k), 用于构造 Q_eff = C'*Q*C
    C = [eye(3), zeros(3,3)];                 % 3×6

    %% =====================================================
    % 2. 代价权重 (论文符号, 系数1; ×2 在第7节构造时统一做)
    %    Q (3x3): 跟踪权重 = diag(30, 30, 1)
    %    R (3x3): 控制权重 = 0.3 * I3
    %    rho:     RSS 强凸参数 = 0.01
    % =====================================================
    Q = diag([w_pos, w_pos, w_psi]);          % 3×3 跟踪权重 (论文符号)
    R = w_control * eye(3);                   % 3×3 控制权重 (论文符号)

    %% =====================================================
    % 3. dim 设置 (HPIPM 字段 nx, nu, ng, nbx; 逐阶段指定)
    %    nx (n=0..K):  6  — 状态 x_k = [e_k; v_k] (n=0 通过 box bounds 固定)
    %    nu (n=0..K-1): 3  — 控制 u_n = v_{n+1} - v_n
    %    nu (n=K):     0  — 末阶段无控制
    %    nbx (n=0):    6  — 初始状态通过 box bounds 固定
    %    nbx (n=1..K): 0
    %    ng (n=0):     2*N (仅转向锥; x_0 固定故无轮速约束)
    %    ng (n=1..K-1): N + 2*N = 3*N (轮速 + 转向锥)
    %    ng (n=K-1):   2*N + 2*N = 4*N (轮速 + 转向锥 + 终端轮速)
    %    ng (n=K):     0 (终端 nu=0, 无约束)
    % =====================================================
    N_stages = K + 1;                         % n=0..K
    ng_wheels = num_wheels;                   % 轮速约束数 (每 stage)
    ng_cone = 2 * num_wheels;                 % 转向锥约束数 (R1+R2, 每 stage)

    %% =====================================================
    % 4. 时变量: b_list, r_list, const
    %    b_list: 动力学 bias 列表 (HPIPM 字段 b), n=0..K-1, 每个 6×1
    %            b_n = [xi_k^ref - xi_{k+1}^ref; 0] (参考轨迹差分, 仅前3维非零)
    %    r_list: 控制线性项列表 (HPIPM 字段 r), n=0..K-1, 每个 3×1
    %            r_n = -2*rho * u_hat(:, n+1)  (已×2; off-by-one: u_hat 用论文编号)
    %    const:  常数项 = rho * sum ||u_hat||^2 (不×2; 不影响 argmin, 仅用于 obj 比较)
    % =====================================================
    const = 0;                                % 累加器初始化 (最终非零, 不×2)
    b_list = cell(1, K);                      % n=0..K-1
    r_list = cell(1, K);
    for n = 0:K-1
        ref_idx_k   = min(size(path,2), step + n);       % xi_k^ref
        ref_idx_kp1 = min(size(path,2), step + n + 1);   % xi_{k+1}^ref
        ref_k   = path(:, ref_idx_k);
        ref_kp1 = path(:, ref_idx_kp1);
        b_xi = ref_k - ref_kp1;                % 3×1, 参考轨迹差分
        b_n = [b_xi; zeros(3,1)];              % 6×1 动力学 bias
        b_list{n+1} = b_n;

        r_n = -2 * rho * u_hat(:, n+1);        % 3×1, 已×2 (off-by-one: u_hat(:,n+1) 是论文编号)
        r_list{n+1} = r_n;

        const = const + rho * (u_hat(:,n+1)' * u_hat(:,n+1));  % rho*||u_hat||^2 (不×2)
    end

    %% =====================================================
    % 5. 初始状态 (通过 box bounds 固定, 不吸收进 b_0)
    %    x_0 = [e_0; v_0^cur], e_0 = xi_cur - xi_0^ref
    %    lbx(0) = ubx(0) = x_0 (HPIPM box bounds)
    % =====================================================
    ref_idx_0 = min(size(path,2), step);       % xi_0^ref
    ref_0 = path(:, ref_idx_0);
    e0 = [current_xy - ref_0(1:2); psi0 - ref_0(3)];   % 3×1
    x0 = [e0; v0];                             % 6×1 (已知, 通过 box bounds 固定)

    %% =====================================================
    % 6. 线性化约束 (ocp_qp) — HPIPM 字段 C, D, lg, ug
    %    所有二次约束在当前线性化点 (v_hat, u_hat) 处做一阶泰勒展开,
    %    转为一般线性约束 D*u + C*x + lg <= 0 <= ug (HPIPM ng 约束).
    % =====================================================
    % 预计算 û 对应的 v̂ 序列 (转向锥凸化 B 项用, 论文 Appendix A)
    %   v̂_k = v0 + Σ_{j=1}^k û_j  (论文编号; 代码 u_hat(:,j) 即论文 û_j)
    %   nu_hat(:,1) = v̂_0 = v0, nu_hat(:,k+1) = v̂_k
    nu_hat = zeros(3, K+1);                    % k=0..K (MATLAB 1..K+1)
    nu_hat(:, 1) = v0;                         % v̂_0 = v0
    for k = 1:K
        nu_hat(:, k+1) = nu_hat(:, k) + u_hat(:, k);
    end

    % 转向锥旋转矩阵 R1, R2 (论文 (12)): R1 = R(pi/2 - delta_theta), R2 = R1^T
    R1 = [sin(delta_theta), -cos(delta_theta); cos(delta_theta),  sin(delta_theta)];
    R2 = [sin(delta_theta),  cos(delta_theta); -cos(delta_theta), sin(delta_theta)];

    % 约束存储 (stage K-1 需要额外 N 条终端轮速约束, 故 ngc_max = 4N)
    ngc_max = 2*ng_wheels + ng_cone;
    Cmat = cell(N_stages, ngc_max);            % 状态系数 (ng×nx)
    Dmat = cell(N_stages, ngc_max);            % 控制系数 (ng×nu)
    lg = cell(N_stages, ngc_max);              % 下界 (ng×1)
    ug = cell(N_stages, ngc_max);              % 上界 (ng×1)
    ng_per_stage = zeros(1, N_stages);         % 每 stage 实际约束数

    % -------------------------------------------------------
    % 6.1 轮速 SOC 约束 (论文公式 20b): ||H_i * v_k||^2 <= vimax^2
    %     线性化: ||H v_hat||^2 + 2*v_hat'M*(v - v_hat) <= vimax^2
    %     注意 ||H v_hat||^2 = v_hat'M v_hat, 所以 RHS = vimax^2
    %     => (2*M*v_hat)' * v <= vimax^2  (线性, 作用在 x 的 v 分量)
    %     约束形式: C*x + D*u + lg <= 0 <= ug
    %       C = [0,0,0, 2*M*v_hat] (1×6), D = zeros(1,3), ug = vimax^2
    % -------------------------------------------------------
    for n = 1:K-1                              % 不含终端 stage K
        v_hat_n = nu_hat(:, n+1);              % v̂_n (MATLAB 1-indexed)
        for i = 1:num_wheels
            Mi = Hn{i}' * Hn{i};               % 3×3
            coeff_v = 2 * Mi * v_hat_n;        % 3×1, grad = 2*M*v_hat
            A_hat_i = v_hat_n' * Mi * v_hat_n; % A_hat = f(v_hat) = v_hat'M v_hat
            C_ni = zeros(1, 6);
            C_ni(4:6) = coeff_v';
            D_ni = zeros(1, 3);
            lg_ni = -1e8;                      % 无下界 (大负数, 避免 IPM 边界问题)
            ug_ni = vimax^2 + A_hat_i;         % 上界: vimax^2 + A_hat (一阶泰勒补偿)

            idx = ng_per_stage(n+1) + 1;
            Cmat{n+1, idx} = C_ni;
            Dmat{n+1, idx} = D_ni;
            lg{n+1, idx} = lg_ni;
            ug{n+1, idx} = ug_ni;
            ng_per_stage(n+1) = idx;
        end
    end

    % -------------------------------------------------------
    % 6.1b 终端轮速约束移到 stage K-1 (同二次约束版本结构)
    %   原终端约束: ||H_i * v_K||^2 <= vimax^2
    %   v_K = v_{K-1} + u_{K-1} (动力学)
    %   线性化在 (v_hat_{K-1}, u_hat_{K-1}):
    %     => 2*M*(v_hat+u_hat)' * (v+u) <= vimax^2
    %     => C = [0,0,0, 2*M*(v_hat+u_hat)], D = [2*M*(v_hat+u_hat)], ug = vimax^2
    % -------------------------------------------------------
    n_term = K - 1;                            % HPIPM stage K-1 (MATLAB 索引 K)
    v_hat_term = nu_hat(:, n_term+1);          % v̂_{K-1}
    u_hat_term = u_hat(:, n_term+1);           % û_{K-1} (论文编号 K)
    w_hat_term = v_hat_term + u_hat_term;      % ŵ = v̂ + û (终端速度线性化点)
    for i = 1:num_wheels
        Mi = Hn{i}' * Hn{i};
        A_hat_term = w_hat_term' * Mi * w_hat_term;  % A_hat = w_hat'M w_hat
        coeff_vu = 2 * Mi * w_hat_term;        % 3×1, grad = 2*M*w_hat
        C_ni = zeros(1, 6);
        C_ni(4:6) = coeff_vu';
        D_ni = coeff_vu';                      % 1×3
        lg_ni = -1e8;
        ug_ni = vimax^2 + A_hat_term;          % 上界: vimax^2 + A_hat (一阶泰勒补偿)

        idx = ng_per_stage(n_term+1) + 1;
        Cmat{n_term+1, idx} = C_ni;
        Dmat{n_term+1, idx} = D_ni;
        lg{n_term+1, idx} = lg_ni;
        ug{n_term+1, idx} = ug_ni;
        ng_per_stage(n_term+1) = idx;
    end

    % -------------------------------------------------------
    % 6.2 Steering cone convexified (paper Prop.1 / eq 15-16):
    %   f(x,u) = x'M(x+u) - 0.5*||Tx+Uu||^2 >= 0
    %   M = H'*Rg*H (NON-SYMMETRIC!), T=(I+Rg)H, U=Rg*H
    %   Linearise at (xh,uh): gx=(M+M')x+M'u-T'Tx-T'Uu, gu=M'x-U'Uu-U'Tx
    %   HPIPM: -gx'x -gu'u <= -(f_hat+gx'xh+gu'uh)
    % -------------------------------------------------------
    for n = 0:K-1
        k = n + 1;
        xh = nu_hat(:, k);        % x_hat = v_hat_{k-1}
        uh = u_hat(:, k);         % u_hat_k
        for i = 1:num_wheels
            Hi = Hn{i};           % 2x3
            for gg = 1:2
                if gg == 1;  Rg = R1;  else;  Rg = R2;  end

                M = Hi' * Rg * Hi;              % 3x3 NON-SYMMETRIC!
                T = (eye(2) + Rg) * Hi;         % 2x3
                U = Rg * Hi;                    % 2x3
                TtT = T' * T;   UtU = U' * U;
                TtU = T' * U;   UtT = U' * T;

                ell = T * xh + U * uh;
                B_const = 0.5 * (ell' * ell);
                f_hat = xh' * M * (xh + uh) - B_const;

                gx = (M + M') * xh + M' * uh  -  TtT * xh  -  TtU * uh;
                gu = M' * xh                  -  UtU * uh   -  UtT * xh;

                rhs = f_hat + gx' * xh + gu' * uh;
                coeff_v = -gx;
                coeff_u = -gu;
                ug_ni   = -rhs;
                lg_ni   = -1e8;

                C_ni = zeros(1, 6);   C_ni(4:6) = coeff_v';
                D_ni = coeff_u';
                idx = ng_per_stage(n+1) + 1;
                Cmat{n+1, idx} = C_ni;  Dmat{n+1, idx} = D_ni;
                lg{n+1, idx} = lg_ni;   ug{n+1, idx} = ug_ni;
                ng_per_stage(n+1) = idx;
            end
        end
    end

    %% =====================================================
    % 7. 返回 (字段名对齐 HPIPM ocp_qp API)
    %    ×2 约定: 构造时统一×2 (Q_eff, R_eff, r 均已×2)
    %    set 时不再×2; const 不×2 (HPIPM const 无 1/2 前缀)
    % =====================================================
    Q_eff = 2 * (C' * Q * C);                  % 6×6 blkdiag(Q, 0_3), 已×2
    R_eff = 2 * (R + rho * eye(3));            % 3×3 (R + rho*I), 已×2

    ocp.A  = A;      ocp.B  = B;                % HPIPM 字段 A, B
    ocp.Q_eff = Q_eff;  ocp.R_eff = R_eff;     % 已×2 (set 时不×2)
    ocp.b  = b_list;                            % HPIPM 字段 b (bias)
    ocp.r  = r_list;                            % HPIPM 字段 r (控制线性项, 已×2)
    ocp.const = const;                         % 不×2 (HPIPM const 无 1/2 前缀)
    ocp.K  = K;      ocp.N_stages = N_stages;
    % dim 设置 (逐 stage):
    ocp.nx = repmat(6, 1, N_stages);           % [6, 6, ..., 6] (所有 stage nx=6)
    ocp.nu = [repmat(3, 1, K), 0];            % [3, 3, ..., 3, 0]
    % ng 分布 (终端轮速约束移到 stage K-1, 终端 stage ng=0):
    %   stage 0:     ng_cone = 2N (仅转向锥; x_0 固定故无轮速约束)
    %   stage 1..K-2: ng_wheels + ng_cone = 3N (轮速 + 转向锥)
    %   stage K-1:   2*ng_wheels + ng_cone = 4N (轮速 + 转向锥 + 终端轮速)
    %   stage K:     0 (终端 nu=0, 无约束)
    if K >= 2
        ng_arr = [ng_cone, repmat(ng_wheels+ng_cone, 1, K-2), 2*ng_wheels+ng_cone, 0];
    else
        % K=1 边界情况: 只有 stage 0 和 stage 1
        ng_arr = [2*ng_wheels+ng_cone, 0];
    end
    ocp.ng = ng_arr;
    ocp.nbx = [6, repmat(0, 1, K)];           % [6, 0, ..., 0] (仅 stage 0 有 box bounds)
    ocp.ng_per_stage = ng_per_stage;
    ocp.Cmat = Cmat; ocp.Dmat = Dmat;         % HPIPM 字段 C, D (一般线性约束系数)
    ocp.lg = lg; ocp.ug = ug;                 % HPIPM 字段 lg, ug (下界/上界)
    ocp.x0 = x0;                               % 初始状态 (box bounds: lbx=ubx=x0)
    ocp.idxbx = [0:5];                         % stage 0 的 box bound 索引 (0-indexed)
end
