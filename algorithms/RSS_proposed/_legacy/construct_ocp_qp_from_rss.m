function ocp = construct_ocp_qp_from_rss(path, step, v0, state, u_anchor, u_cut, params)
% 修改: 分离 u_anchor (RSS 凸化锚点 B/L/r/const) 和 u_cut (线性化切点)
% 模式 B (bounded multi-cut): u_cut 可以是 cell array {3×K, 3×K, ...}, 最多 3 组
% 向后兼容: 若只传 6 参数, u_cut = u_anchor (退化为旧行为)
if nargin < 7
    params = u_cut;
    u_cut = u_anchor;
end
% 统一 u_cut 为 cell array (模式 A: 1 组; 模式 B: 多组)
if ~iscell(u_cut)
    u_cut = {u_cut};
end
n_cuts = numel(u_cut);  % cut point 数量
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
%     三次 OCP QP 是对原凸二次子问题的有限切平面近似;
%     可在特定问题上接近 Dense QCQP, 但不保证一般性严格等价.

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

        r_n = -2 * rho * u_anchor(:, n+1);    % 3×1, 已×2 (RSS 锚点 u_anchor)
        r_list{n+1} = r_n;

        const = const + rho * (u_anchor(:,n+1)' * u_anchor(:,n+1));  % rho*||u_anchor||^2 (不×2)
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
    % u_anchor 的 nu_hat 序列 (用于 B, L 项: RSS 凸化锚点)
    nu_hat_anchor = zeros(3, K+1);             % k=0..K (MATLAB 1..K+1)
    nu_hat_anchor(:, 1) = v0;                  % v̂_0 = v0
    for k = 1:K
        nu_hat_anchor(:, k+1) = nu_hat_anchor(:, k) + u_anchor(:, k);
    end
    % u_cut 的 nu_hat 序列 (用于线性化切点: A 项切平面)
    % 模式 B: 对每个 cut point 计算 nu_hat_cut
    nu_hat_cut_list = cell(1, n_cuts);
    for c = 1:n_cuts
        nu_hat_cut_c = zeros(3, K+1);
        nu_hat_cut_c(:, 1) = v0;
        for k = 1:K
            nu_hat_cut_c(:, k+1) = nu_hat_cut_c(:, k) + u_cut{c}(:, k);
        end
        nu_hat_cut_list{c} = nu_hat_cut_c;
    end

    % 转向锥旋转矩阵 R1, R2 (论文 (12)): R1 = R(pi/2 - delta_theta), R2 = R1^T
    R1 = [sin(delta_theta), -cos(delta_theta); cos(delta_theta),  sin(delta_theta)];
    R2 = [sin(delta_theta),  cos(delta_theta); -cos(delta_theta), sin(delta_theta)];

    % 约束存储 (模式 B: 每 stage 约束数 = n_cuts * (2*ng_wheels + ng_cone))
    ngc_max = n_cuts * (2*ng_wheels + ng_cone);
    Cmat = cell(N_stages, ngc_max);            % 状态系数 (ng×nx)
    Dmat = cell(N_stages, ngc_max);            % 控制系数 (ng×nu)
    lg = cell(N_stages, ngc_max);              % 下界 (ng×1)
    ug = cell(N_stages, ngc_max);              % 上界 (ng×1)
    ng_per_stage = zeros(1, N_stages);         % 每 stage 实际约束数

    % -------------------------------------------------------
    % 6.1 轮速 SOC 约束 (论文公式 20b): ||H_i * v_k||^2 <= vimax^2
    %     原约束: v'M v <= vimax^2  (M = H_i' * H_i)
    %     一阶泰勒线性化在 v_hat 处:
    %       v_hat'M v_hat + 2*v_hat'M*(v - v_hat) <= vimax^2
    %     => 2*v_hat'M * v <= vimax^2 + v_hat'M*v_hat
    %     即 RHS = vimax^2 + A_hat, 其中 A_hat = v_hat'M*v_hat
    %     约束形式: C*x + D*u + lg <= 0 <= ug
    %       C = [0,0,0, 2*M*v_hat] (1×6), D = zeros(1,3), ug = vimax^2 + A_hat
    % -------------------------------------------------------
    for c = 1:n_cuts
        nu_hat_cut = nu_hat_cut_list{c};       % 当前 cut point 的 nu_hat 序列
        u_cut_c = u_cut{c};                    % 当前 cut point 的 u
    for n = 1:K-1                              % 不含终端 stage K
        v_hat_n = nu_hat_cut(:, n+1);          % v̂_n (切点, from u_cut)
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
    end % for c = 1:n_cuts

    % -------------------------------------------------------
    % 6.1b 终端轮速约束移到 stage K-1 (同二次约束版本结构)
    %   原终端约束: ||H_i * v_K||^2 <= vimax^2
    %   v_K = v_{K-1} + u_{K-1} (动力学)
    %   令 w = v + u (终端速度), 原约束: w'M w <= vimax^2
    %   一阶泰勒线性化在 w_hat = v_hat + u_hat 处:
    %     2*w_hat'M * w <= vimax^2 + w_hat'M*w_hat
    %   即 RHS = vimax^2 + A_hat_term, 其中 A_hat_term = w_hat'M*w_hat
    %     => C = [0,0,0, 2*M*w_hat], D = [2*M*w_hat], ug = vimax^2 + A_hat_term
    % -------------------------------------------------------
    for c = 1:n_cuts
        nu_hat_cut = nu_hat_cut_list{c};
        u_cut_c = u_cut{c};
    n_term = K - 1;                            % HPIPM stage K-1 (MATLAB 索引 K)
    v_hat_term = nu_hat_cut(:, n_term+1);      % v̂_{K-1} (切点, from u_cut)
    u_hat_term = u_cut_c(:, n_term+1);         % u_cut_{K-1} (切点)
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
    end % for c = 1:n_cuts

    % -------------------------------------------------------
    % 6.2 Steering cone convexified (paper Prop.1 / eq 15-16):
    %   C(x,u) = A(x,u) - B(u_hat) - L(u,u_hat) <= 0  (凸二次约束)
    %   A = 0.5*||H*x||^2 + 0.5*||H*(x+u)||^2  (凸, Mn = H*H 对称)
    %   B = 0.5*||(I+R)*H*xh + R*H*uh||^2      (常数, 在 u_hat 处)
    %   L = b*(u - uh),  b = grad_u B|_uh = U*(T*xh + U*uh)  (线性)
    %   T = (I+R)*H, U = R*H
    %
    %   线性化 C 在 (xh,uh) 处 (一阶泰勒, 用于 OCP QP 线性约束):
    %     注意: SCP 中线性化点 (xh,uh) = 上一轮 u_hat, 故 L(uh,uh) = 0
    %     C(xh,uh) = A(xh,uh) - B_const
    %     grad_x C = 2*Mn*xh + Mn*uh
    %     grad_u C = Mn*(xh+uh) - b
    %     rhs = grad_x_C*xh + grad_u_C*uh - C(xh,uh)
    %   HPIPM: coeff_v*v + coeff_u*u <= rhs
    % -------------------------------------------------------
    for c = 1:n_cuts
        nu_hat_cut = nu_hat_cut_list{c};
        u_cut_c = u_cut{c};
    for n = 0:K-1
        k = n + 1;
        x_anchor = nu_hat_anchor(:, k);  % 锚点 v̂_{k-1} (from u_anchor, 用于 B/L)
        u_anchor_k = u_anchor(:, k);     % 锚点 û_k (from u_anchor, 用于 B/L)
        xh = nu_hat_cut(:, k);           % 切点 v_cut_{k-1} (from u_cut, 用于 A 切平面)
        uh = u_cut_c(:, k);              % 切点 u_cut_k (from u_cut, 用于 A 切平面)
        for i = 1:num_wheels
            Hi = Hn{i};           % 2x3
            Mn = Hi' * Hi;        % 3x3 对称 (M_n = H_n*H_n, 与 dense QCQP 一致)
            for gg = 1:2
                if gg == 1;  Rg = R1;  else;  Rg = R2;  end

                T = (eye(2) + Rg) * Hi;         % 2x3
                U = Rg * Hi;                    % 2x3

                % b = T*xh + U*uh  (论文 b, 在 u_hat 处)
                ell_hat = T * x_anchor + U * u_anchor_k;  % (I+R)*H*x_anchor + R*H*u_anchor_k (锚点)
                B_const = 0.5 * (ell_hat' * ell_hat);  % B 项 (常数)

                % L 项 = b'*T*(x_v - xh) + b'*U*(u - uh)  (论文公式 16 L)
                %   利用 Sigma_{l<k}(u_l - u_hat_l) = nu_{k-1} - nu_hat_{k-1} = x_v - xh
                %   L 同时依赖 x 和 u!
                grad_xv_L = T' * ell_hat;       % grad_{x_v} L = T'*b, 3x1
                grad_u_L  = U' * ell_hat;       % grad_u L = U'*b, 3x1

                % A 项 (凸二次, 在切点 (xh,uh)=(v_cut,u_cut) 处求值)
                A_xh_uh = xh'*Mn*xh + xh'*Mn*uh + 0.5*uh'*Mn*uh;

                % L 项 (在切点处, 锚点固定为 u_anchor; 当 u_cut≠u_anchor 时 L≠0)
                %   L = ell_anchor' * (T*(x_cut - x_anchor) + U*(u_cut - u_anchor_k))
                L_cut = ell_hat' * (T*(xh - x_anchor) + U*(uh - u_anchor_k));

                % C(xh,uh) = A(xh,uh) - B_const - L_cut
                C_xh_uh = A_xh_uh - B_const - L_cut;

                % C 的梯度: grad C = grad A - grad L (B 是常数)
                grad_x_C = 2*Mn*xh + Mn*uh - grad_xv_L;  % grad_x A - grad_x L
                grad_u_C = Mn*(xh+uh) - grad_u_L;         % grad_u A - grad_u L

                % 线性化: grad_x_C'*x + grad_u_C'*u <= rhs
                rhs = grad_x_C'*xh + grad_u_C'*uh - C_xh_uh;
                coeff_v = grad_x_C;             % 3x1
                coeff_u = grad_u_C;             % 3x1
                ug_ni   = rhs;
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
    end % for c = 1:n_cuts

    %% =====================================================
    % 7. 返回 (字段名对齐 HPIPM ocp_qp API)
    %    ×2 约定: 构造时统一×2 (Q_eff, R_eff, r 均已×2)
    %    set 时不再×2; const 不×2 (HPIPM const 无 1/2 前缀)
    % =====================================================
    % Q_eff 逐 stage 不同 (与 Dense QCQP construct_complete_qp_from_rss.m 一致):
    %   Dense QCQP: 位置代价 k=2..K, 姿态代价 k=1..K
    %   OCP QP stage n 对应论文 k=n:
    %     stage 0 (k=0): Q=0 (e_0 是已知常数, 不加代价)
    %     stage 1 (k=1): Q = 2*diag(0, 0, w_psi, 0, 0, 0) (仅姿态)
    %     stage 2..K:    Q = 2*diag(w_pos, w_pos, w_psi, 0, 0, 0) (位置+姿态)
    % 用 2D 堆叠 (6, 6*N_stages) 避免 MATLAB 3D→numpy 3D 维度转置问题
    Q_pos_psi = 2 * diag([w_pos, w_pos, w_psi, 0, 0, 0]);  % stage 2..K (已×2)
    Q_psi_only = 2 * diag([0, 0, w_psi, 0, 0, 0]);         % stage 1 (已×2)
    Q_zero = zeros(6);                                      % stage 0

    Q_eff = zeros(6, 6*N_stages);              % 2D 堆叠: 每 6 列为一个 stage 的 Q
    Q_eff(:, 1:6) = Q_zero;                    % stage 0 (k=0)
    Q_eff(:, 7:12) = Q_psi_only;               % stage 1 (k=1)
    for s = 3:N_stages
        Q_eff(:, (s-1)*6+1:s*6) = Q_pos_psi;   % stage 2..K (k=2..K)
    end

    R_eff = 2 * (R + rho * eye(3));            % 3×3 (R + rho*I), 已×2
    S_eff = zeros(3, 6);                        % HPIPM S (nu*nx), no cross term, no x2
    q_stack = zeros(6, N_stages);               % HPIPM q (nx*N_stages), no linear term, no x2

    ocp.A  = A;      ocp.B  = B;                % HPIPM 字段 A, B
    ocp.Q_eff = Q_eff;  ocp.R_eff = R_eff;     % Q_eff 为 2D (6, 6*N_stages), 已×2
    ocp.S_eff = S_eff;                          % HPIPM S (nu*nx), no x2
    ocp.q_stack = q_stack;                      % HPIPM q (nx*N_stages), no x2
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
        ng_arr = n_cuts * [ng_cone, repmat(ng_wheels+ng_cone, 1, K-2), 2*ng_wheels+ng_cone, 0];
    else
        % K=1 边界情况: 只有 stage 0 和 stage 1
        ng_arr = n_cuts * [2*ng_wheels+ng_cone, 0];
    end
    ocp.ng = ng_arr;
    ocp.nbx = [6, repmat(0, 1, K)];           % [6, 0, ..., 0] (仅 stage 0 有 box bounds)
    ocp.ng_per_stage = ng_per_stage;
    ocp.Cmat = Cmat; ocp.Dmat = Dmat;         % HPIPM 字段 C, D (一般线性约束系数)
    ocp.lg = lg; ocp.ug = ug;                 % HPIPM 字段 lg, ug (下界/上界)
    ocp.x0 = x0;                               % 初始状态 (box bounds: lbx=ubx=x0)
    ocp.idxbx = [0:5];                         % stage 0 的 box bound 索引 (0-indexed)
end
