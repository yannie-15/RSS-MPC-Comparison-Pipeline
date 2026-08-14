function ocp = construct_ocp_qcqp_from_rss(path, step, v0, state, u_anchor, params)
% CONSTRUCT_OCP_QCQP_FROM_RSS
% 误差状态 OCP QCQP 构造 (HPIPM ocp_qcqp Python 接口, 精确凸二次约束)
%
% 与 construct_complete_qp_from_rss.m (Dense QCQP oracle) 数学等价,
% 但利用 OCP 块三对角结构, 不需要 dense 化.
%
% 状态 x_n = [e_n; v_n] ∈ R^6, n=0..K (K=6, 7 stages)
%   e_n = ξ_n - ξ_n^ref (跟踪误差, 3 dim)
%   v_n = ν_n (车体系速度, 3 dim)
% 控制 u_n ∈ R^3, n=0..K-1 (u_n = v_{n+1} - v_n, HPIPM 编号 = 论文 u_{n+1})
%
% 动力学: x_{n+1} = A x_n + B u_n + b_n (LTI, 时变 bias)
%
% HPIPM OCP QCQP 标准形式 (每 stage n=0..N):
%   动力学:  x_{n+1} = A_n x_n + B_n u_n + b_n      (n=0..K-1)
%   代价:    min Σ_{n=0}^{K-1} [0.5 x'Q_n x + x'S_n' u + 0.5 u'R_n u + q_n'x + r_n'u]
%            + 0.5 x_K'Q_K x_K + q_K'x_K
%   二次约束(每 stage nq 条): 0.5 x'Qq x + x'Sq' u + 0.5 u'Rq u + qq'x + rq'u <= uq
%   box 约束: lbx <= x[idxbx] <= ubx
%
% 约定 (与 construct_ocp_qp_from_rss.m 一致):
%   - Q, R: 构造时×2 (HPIPM 0.5 前缀使其还原)
%   - S, q, r: 不×2 (HPIPM 无 0.5 前缀)
%   - Qq: ×2 (HPIPM 0.5 前缀), Rq: ×2 (HPIPM 0.5 前缀)
%   - Sq, qq, rq: 不×2 (HPIPM 无 0.5 前缀)
%   - uq, const: 不×2
%
% 72 条二次约束的 stage 分布:
%   nq_per_stage = [8, 12, 12, 12, 12, 12, 4]  (stage 0..6)
%   stage 0:     8 steering (k=1, 4 wheels × 2 rotations)
%   stage 1..5:  4 wheel (k=stage, on x_stage) + 8 steering (k=stage+1) = 12
%   stage 6:     4 wheel (k=6, on x_6, terminal)
%   总计: 8 + 5*12 + 4 = 72

    %% =====================================================
    % 参数提取 (论文 IV-A 实验设置, 与 Dense QCQP 一致)
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

    % 转向锥旋转矩阵 R1, R2 (论文 (12)): R1 = R(pi/2 - delta_theta), R2 = R1^T
    R1 = [sin(delta_theta), -cos(delta_theta); cos(delta_theta),  sin(delta_theta)];
    R2 = [sin(delta_theta),  cos(delta_theta); -cos(delta_theta), sin(delta_theta)];

    %% =====================================================
    % 1. 动力学矩阵 (常数, LTI, 与 construct_ocp_qp_from_rss.m 一致)
    % =====================================================
    A  = [eye(3), [R_psi0*dt, zeros(2,1); zeros(1,2), dt]; zeros(3,3), eye(3)];  % 6×6
    B  = [zeros(3,3); eye(3)];                % 6×3 控制矩阵

    %% =====================================================
    % 2. dim 设置
    % =====================================================
    N_stages = K + 1;                         % n=0..K (7 stages)
    nx = 6;
    nu = 3;

    % nq_per_stage = [8, 12, 12, 12, 12, 12, 4]
    nq_per_stage = zeros(1, N_stages);
    nq_per_stage(1) = 2 * num_wheels;         % stage 0: 8 steering
    for s = 2:K                                % stage 1..K-1: 4 wheel + 8 steering = 12
        nq_per_stage(s) = num_wheels + 2 * num_wheels;
    end
    nq_per_stage(K+1) = num_wheels;           % stage K: 4 wheel (terminal)
    total_nq = sum(nq_per_stage);             % 72

    %% =====================================================
    % 3. 时变量: b_list, r_list, const
    %    (与 construct_ocp_qp_from_rss.m 完全一致)
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

        r_n = -2 * rho * u_anchor(:, n+1);    % 3×1, RSS 锚点 u_anchor (不×2, 自然系数)
        r_list{n+1} = r_n;

        const = const + rho * (u_anchor(:,n+1)' * u_anchor(:,n+1));  % rho*||u_anchor||^2 (不×2)
    end

    %% =====================================================
    % 4. 初始状态 (通过 box bounds 固定)
    %    x_0 = [e_0; v_0], e_0 = xi_cur - xi_0^ref
    % =====================================================
    ref_idx_0 = min(size(path,2), step);       % xi_0^ref
    ref_0 = path(:, ref_idx_0);
    e0 = [current_xy - ref_0(1:2); psi0 - ref_0(3)];   % 3×1
    x0 = [e0; v0];                             % 6×1 (已知, 通过 box bounds 固定)

    %% =====================================================
    % 5. 逐 stage 代价矩阵 (与 construct_ocp_qp_from_rss.m 一致)
    %    stage 0 (k=0): Q=0 (e_0 是已知常数, 不加代价)
    %    stage 1 (k=1): Q = 2*diag(0, 0, w_psi, 0, 0, 0) (仅姿态)
    %    stage 2..K:    Q = 2*diag(w_pos, w_pos, w_psi, 0, 0, 0) (位置+姿态)
    %    R = 2*(w_control*I + rho*I), r = -2*rho*u_anchor, S=0, q=0
    % =====================================================
    Q_pos_psi = 2 * diag([w_pos, w_pos, w_psi, 0, 0, 0]);  % stage 2..K (已×2)
    Q_psi_only = 2 * diag([0, 0, w_psi, 0, 0, 0]);         % stage 1 (已×2)
    Q_zero = zeros(6);                                      % stage 0

    Q_stack = zeros(nx, nx * N_stages);         % 2D 堆叠: 每 6 列为一个 stage 的 Q
    Q_stack(:, 1:6) = Q_zero;                   % stage 0
    Q_stack(:, 7:12) = Q_psi_only;              % stage 1
    for s = 3:N_stages
        Q_stack(:, (s-1)*nx+1:s*nx) = Q_pos_psi;   % stage 2..K
    end

    R_eff = 2 * (w_control * eye(nu) + rho * eye(nu));  % 3×3 (R + rho*I), 已×2
    S_eff = zeros(nu, nx);                               % HPIPM S (nu*nx), 无交叉项
    q_stack = zeros(nx, N_stages);                       % HPIPM q, 无线性项

    %% =====================================================
    % 6. 预计算 u_anchor 对应的 v̂ 序列 (转向锥凸化 B/L 项用)
    %    v̂_k = v0 + Σ_{j=1}^k û_j (论文 Appendix A)
    %    nu_hat_anchor(:, 1) = v̂_0 = v0
    %    nu_hat_anchor(:, k+1) = v̂_k = v̂_{k-1} + û_k
    % =====================================================
    nu_hat_anchor = zeros(3, K+1);             % k=0..K (MATLAB 1..K+1)
    nu_hat_anchor(:, 1) = v0;                  % v̂_0 = v0
    for k = 1:K
        nu_hat_anchor(:, k+1) = nu_hat_anchor(:, k) + u_anchor(:, k);
    end

    %% =====================================================
    % 7. 构造 72 条精确二次约束
    % =====================================================
    % 约束 stacks (3D/2D, 按 stage 顺序排列)
    Qq_stack = zeros(nx, nx, total_nq);         % (6, 6, 72)
    Sq_stack = zeros(nu, nx, total_nq);         % (3, 6, 72) HPIPM 原生方向 (nu, nx)
    Rq_stack = zeros(nu, nu, total_nq);         % (3, 3, 72)
    qq_stack = zeros(nx, total_nq);             % (6, 72)
    rq_stack = zeros(nu, total_nq);             % (3, 72)
    uq_stack = zeros(total_nq, 1);              % (72, 1)

    % 元数据 (用于 Dense/OCP 约束对应)
    metadata_kind = cell(total_nq, 1);          % 'wheel' / 'steering'
    metadata_k = zeros(total_nq, 1);            % paper step k (1..K)
    metadata_wheel = zeros(total_nq, 1);        % wheel index (1..4)
    metadata_rotation = cell(total_nq, 1);      % 'R1' / 'R2' / '' (wheel)
    metadata_stage = zeros(total_nq, 1);        % OCP stage (0..K)
    metadata_local_index = zeros(total_nq, 1);  % local index within stage

    idx = 0;  % 全局约束索引 (0-based, MATLAB 中用 idx+1)

    % -------------------------------------------------------
    % 7.1 Stage 0: 8 steering 约束 (k=1, involves x_0 and u_0)
    % -------------------------------------------------------
    s = 0;  % OCP stage 0
    k = 1;  % paper step k=1
    v_hat = nu_hat_anchor(:, k);         % v̂_0 = v0 (anchor velocity at stage 0)
    u_hat_k = u_anchor(:, k);            % û_1 (anchor control)
    for n = 1:num_wheels
        H_n = Hn{n};
        M = H_n' * H_n;                  % 3×3 对称
        for gg = 1:2
            if gg == 1; Rg = R1; rot_name = 'R1'; else; Rg = R2; rot_name = 'R2'; end
            T = (eye(2) + Rg) * H_n;     % 2×3
            U = Rg * H_n;                % 2×3
            ell = T * v_hat + U * u_hat_k;  % 2×1

            idx = idx + 1;
            % Qq: ×2 (0.5 前缀), Qq(4:6,4:6) = 2*M
            Qq = zeros(nx, nx);
            Qq(4:6, 4:6) = 2 * M;
            Qq_stack(:, :, idx) = Qq;
            % Sq: 不×2 (无 0.5 前缀), Sq(:,4:6) = M
            Sq = zeros(nu, nx);
            Sq(:, 4:6) = M;
            Sq_stack(:, :, idx) = Sq;
            % Rq: ×2 (0.5 前缀), A 项系数 0.5 → Rq = 2*0.5*M = M
            Rq = M;
            Rq_stack(:, :, idx) = Rq;
            % qq: 不×2, qq(4:6) = -T'*ell
            qq = zeros(nx, 1);
            qq(4:6) = -T' * ell;
            qq_stack(:, idx) = qq;
            % rq: 不×2, rq = -U'*ell
            rq = -U' * ell;
            rq_stack(:, idx) = rq;
            % uq: 不×2, uq = -0.5*ell'ell
            uq_stack(idx) = -0.5 * (ell' * ell);

            % metadata
            metadata_kind{idx} = 'steering';
            metadata_k(idx) = k;
            metadata_wheel(idx) = n;
            metadata_rotation{idx} = rot_name;
            metadata_stage(idx) = s;
            metadata_local_index(idx) = idx;
        end
    end

    % -------------------------------------------------------
    % 7.2 Stage 1..K-1: 4 wheel (k=stage) + 8 steering (k=stage+1)
    % -------------------------------------------------------
    for s = 1:K-1
        % --- 4 wheel 约束 (k=s, on x_s) ---
        k_wheel = s;  % paper step k = stage
        for n = 1:num_wheels
            H_n = Hn{n};
            M = H_n' * H_n;

            idx = idx + 1;
            % Qq: ×2, Qq(4:6,4:6) = 2*M
            Qq = zeros(nx, nx);
            Qq(4:6, 4:6) = 2 * M;
            Qq_stack(:, :, idx) = Qq;
            % Sq, Rq, rq: 全零 (轮速约束不依赖 u)
            % (Sq_stack, Rq_stack, rq_stack 已初始化为零)
            % qq: 全零
            % (qq_stack 已初始化为零)
            % uq: vimax^2
            uq_stack(idx) = vimax^2;

            % metadata
            metadata_kind{idx} = 'wheel';
            metadata_k(idx) = k_wheel;
            metadata_wheel(idx) = n;
            metadata_rotation{idx} = '';
            metadata_stage(idx) = s;
            metadata_local_index(idx) = idx - sum(nq_per_stage(1:s));
        end

        % --- 8 steering 约束 (k=s+1, involves x_s and u_s) ---
        k_steer = s + 1;  % paper step k = stage + 1
        v_hat = nu_hat_anchor(:, k_steer);     % v̂_s (anchor velocity at stage s)
        u_hat_k = u_anchor(:, k_steer);        % û_{s+1} (anchor control)
        for n = 1:num_wheels
            H_n = Hn{n};
            M = H_n' * H_n;
            for gg = 1:2
                if gg == 1; Rg = R1; rot_name = 'R1'; else; Rg = R2; rot_name = 'R2'; end
                T = (eye(2) + Rg) * H_n;
                U = Rg * H_n;
                ell = T * v_hat + U * u_hat_k;

                idx = idx + 1;
                Qq = zeros(nx, nx);
                Qq(4:6, 4:6) = 2 * M;
                Qq_stack(:, :, idx) = Qq;
                Sq = zeros(nu, nx);
                Sq(:, 4:6) = M;
                Sq_stack(:, :, idx) = Sq;
                Rq = M;
                Rq_stack(:, :, idx) = Rq;
                qq = zeros(nx, 1);
                qq(4:6) = -T' * ell;
                qq_stack(:, idx) = qq;
                rq = -U' * ell;
                rq_stack(:, idx) = rq;
                uq_stack(idx) = -0.5 * (ell' * ell);

                metadata_kind{idx} = 'steering';
                metadata_k(idx) = k_steer;
                metadata_wheel(idx) = n;
                metadata_rotation{idx} = rot_name;
                metadata_stage(idx) = s;
                metadata_local_index(idx) = idx - sum(nq_per_stage(1:s));
            end
        end
    end

    % -------------------------------------------------------
    % 7.3 Stage K: 4 wheel 约束 (k=K, on x_K, terminal)
    % -------------------------------------------------------
    s = K;  % terminal stage
    k_wheel = K;
    for n = 1:num_wheels
        H_n = Hn{n};
        M = H_n' * H_n;

        idx = idx + 1;
        Qq = zeros(nx, nx);
        Qq(4:6, 4:6) = 2 * M;
        Qq_stack(:, :, idx) = Qq;
        % 终端 stage nu=0: Sq/Rq/rq 不设置 (保持零, Python 端跳过)
        % qq: 全零
        % uq: vimax^2
        uq_stack(idx) = vimax^2;

        metadata_kind{idx} = 'wheel';
        metadata_k(idx) = k_wheel;
        metadata_wheel(idx) = n;
        metadata_rotation{idx} = '';
        metadata_stage(idx) = s;
        metadata_local_index(idx) = idx - sum(nq_per_stage(1:s));
    end

    % 校验总数
    assert(idx == total_nq, '约束总数不匹配: idx=%d, total_nq=%d', idx, total_nq);

    %% =====================================================
    % 8. 返回 OCP QCQP 结构体
    % =====================================================
    ocp.A  = A;      ocp.B  = B;
    ocp.Q_eff = Q_stack;  ocp.R_eff = R_eff;
    ocp.S_eff = S_eff;    ocp.q_stack = q_stack;
    ocp.b  = b_list;      ocp.r  = r_list;
    ocp.const = const;
    ocp.K  = K;      ocp.N_stages = N_stages;

    % dim 设置 (逐 stage)
    ocp.nx = repmat(nx, 1, N_stages);           % [6, 6, ..., 6]
    ocp.nu = [repmat(nu, 1, K), 0];             % [3, 3, ..., 3, 0]
    ocp.nbx = [nx, repmat(0, 1, K)];            % [6, 0, ..., 0]
    ocp.nq = nq_per_stage;                      % [8, 12, 12, 12, 12, 12, 4]
    ocp.nq_per_stage = nq_per_stage;

    % 二次约束 stacks
    ocp.Qq_stack = Qq_stack;
    ocp.Sq_stack = Sq_stack;
    ocp.Rq_stack = Rq_stack;
    ocp.qq_stack = qq_stack;
    ocp.rq_stack = rq_stack;
    ocp.uq_stack = uq_stack;

    % 初始状态
    ocp.x0 = x0;
    ocp.idxbx = [0:5];                          % stage 0 的 box bound 索引 (0-indexed)

    % 元数据 (用于 Dense/OCP 约束对应)
    ocp.metadata.kind = metadata_kind;
    ocp.metadata.k = metadata_k;
    ocp.metadata.wheel = metadata_wheel;
    ocp.metadata.rotation = metadata_rotation;
    ocp.metadata.stage = metadata_stage;
    ocp.metadata.local_index = metadata_local_index;
end
