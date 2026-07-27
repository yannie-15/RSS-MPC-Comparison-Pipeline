function qp_problem = construct_complete_qp_from_rss(path, step, current_nu, state, u_hat, params)
% CONSTRUCT_COMPLETE_QP_FROM_RSS
% 从RSS问题显式构造完整的QCQP问题
%
% 问题形式 (QCQP)：
%   min  0.5*x'*H*x + g'*x
%   s.t. A*x = b                 (动力学等式约束)
%        C*x <= d                (线性不等式约束)
%        0.5*x'*Hq_i*x + gq_i'*x <= uq_i  (二次不等式约束)
%
% 二次不等式约束包括：
%   - 轮速约束: ||H{n}*nu(:,k)||^2 <= vimax^2
%   - 转向锥凸化约束 (两组)
%
% 变量排列: x = [u(:); nu(:)]
%   u(i,k) 的全局索引 = (k-1)*3 + i         (i=1,2,3; k=1,...,K)
%   nu(i,k) 的全局索引 = nu_start + (k-1)*3 + (i-1)
%   其中 nu_start = 3*K + 1
%
% 修复记录:
%   [Bug1] nu索引从 +1/+2/+3 改为 +0/+1/+2 (0-based偏移)
%   [Bug2] 加入跨预测阶段的Hessian交叉项
%   [Bug3] 位置代价常数部分加入 current_nu(1:2)*dt 位移
%   [Bug4] Hessian系数改为 2*w*dt^2 (0.5*x'Hx形式)
%   [Bug5] 等式约束递推段索引从 +i 改为 +(i-1)

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

    % 代价函数权重 (与 control_RSS.m / control_RSS_v2.m CVX分支完全一致)
    w_pos = 30;
    w_psi = 1;        % k1 in original
    w_control = 0.3;
    rho = 0.01;

    % 轮子特征矩阵 H{n} (2x3)
    H_cell = cell(1, num_wheels);
    for n = 1:num_wheels
        H_cell{n} = [
            1, 0, -wheel_pos(n, 2);
            0, 1,  wheel_pos(n, 1)
        ];
    end

    % M{n} = H{n}'*H{n} (3x3) 用于二次约束
    M_cell = cell(1, num_wheels);
    for n = 1:num_wheels
        M_cell{n} = H_cell{n}' * H_cell{n};
    end

    delta_theta = dt * phidotmax;

    %% =====================================================
    % 决策变量排列：x = [u(:); nu(:)]
    % u 的全局索引: (k-1)*3 + i        (i=1,2,3; k=1..K)
    % nu 的全局索引: nu_start + (k-1)*3 + (i-1)
    % ======================================================

    n_var = 6 * K;
    nu_start = 3*K + 1;  % nu(1,1) 的全局索引 = 19 (K=6时)

    %% =====================================================
    % 1. Hessian矩阵 H 和一次项 g (代价函数)
    % ======================================================
    % 目标函数: 0.5*x'*H*x + g'*x 应等价于 CVX 中的
    %   J_pos + J_psi + w_control*||u||^2 + rho*||u-u_hat||^2
    %
    % CVX 的 sum_square(v) = ||v||^2 = v'*v
    % 在 0.5*x'Hx 形式中: ||v||^2 = 0.5*x'*(2*A)'*...
    % 所以 Hessian 系数需要乘 2 才能匹配 0.5*x'Hx 形式

    H_mat = zeros(n_var, n_var);
    g_vec = zeros(n_var, 1);

    % 1.1 位置跟踪代价: sum_{k=2}^{K} w_pos * ||position_error_k||^2
    %
    % CVX 表达式:
    %   NU(:,1) = current_nu(1:2)*dt;                (常数，不含nu变量)
    %   NU(:,k) = current_nu(1:2)*dt + dt*sum_{j=1}^{k-1} nu(1:2,j)   (k>=2)
    %   position_error_k = current_xy - ref_xy_k + R_psi0*NU(:,k)
    %
    % 展开:
    %   position_error_k = c_k + R_psi0*dt*S_k
    %   其中 c_k = current_xy - ref_xy_k + R_psi0*current_nu(1:2)*dt  (常数)
    %         S_k = sum_{j=1}^{k-1} nu(1:2,j)                          (nu变量)
    %
    % 代价 = w_pos * ||c_k + R_psi0*dt*S_k||^2
    %       = w_pos*(||c_k||^2 + 2*c_k'*R_psi0*dt*S_k + dt^2*||S_k||^2)
    %
    % 由于 R_psi0 是旋转矩阵, R_psi0'*R_psi0 = I
    % ||R_psi0*dt*S_k||^2 = dt^2*||S_k||^2 = dt^2*(sum_x)^2 + dt^2*(sum_y)^2
    %
    % 在 0.5*x'Hx 形式:
    %   二次项: dt^2*(sum_x)^2 = dt^2*sum_{i,j} nu_x(i)*nu_x(j)
    %           → H(nu_x(i),nu_x(j)) = 2*w_pos*dt^2  (对所有 i,j 对)
    %           同理 nu_y 对: H(nu_y(i),nu_y(j)) = 2*w_pos*dt^2
    %           nu_x 与 nu_y 无交叉 (因为 ||S||^2 = sum_x^2 + sum_y^2)
    %
    %   一次项: 2*w_pos*dt*(R_psi0'*c_k)' * S_k
    %           → g(nu_x(j)) += 2*w_pos*dt*(R_psi0'*c_k)(1)  (j=1..k-1)
    %           → g(nu_y(j)) += 2*w_pos*dt*(R_psi0'*c_k)(2)  (j=1..k-1)

    for k = 2:K
        ref_idx = min(size(path, 2), step + k);
        ref_xy = path(1:2, ref_idx);

        % [Bug3修复] 常数部分必须包含 current_nu(1:2)*dt 产生的位移
        c_k = current_xy - ref_xy + R_psi0 * current_nu(1:2) * dt;

        % 计算 R_psi0' * c_k (梯度方向)
        grad_dir = R_psi0' * c_k;

        % NU(:,k) 依赖于 nu(1:2, 1:k-1)
        for i = 1:k-1
            for j = 1:k-1
                % [Bug2+Bug4修复] 加入所有跨阶段交叉项，系数 = 2*w_pos*dt^2
                nu_x_i = nu_start + (i-1)*3;       % [Bug1修复] 偏移0，不是1
                nu_x_j = nu_start + (j-1)*3;
                nu_y_i = nu_start + (i-1)*3 + 1;   % [Bug1修复] 偏移1，不是2
                nu_y_j = nu_start + (j-1)*3 + 1;

                H_mat(nu_x_i, nu_x_j) = H_mat(nu_x_i, nu_x_j) + 2 * w_pos * dt^2;
                H_mat(nu_y_i, nu_y_j) = H_mat(nu_y_i, nu_y_j) + 2 * w_pos * dt^2;
            end
        end

        % 一次项
        for j = 1:k-1
            nu_x_j = nu_start + (j-1)*3;
            nu_y_j = nu_start + (j-1)*3 + 1;

            g_vec(nu_x_j) = g_vec(nu_x_j) + 2 * w_pos * dt * grad_dir(1);
            g_vec(nu_y_j) = g_vec(nu_y_j) + 2 * w_pos * dt * grad_dir(2);
        end
    end

    % 1.2 姿态跟踪代价: sum_{k=1}^{K} w_psi * (psi_k - ref_psi_k)^2
    %
    % CVX 表达式:
    %   psi(1) = psi0 + current_nu(3)*dt;          (常数，不含nu变量)
    %   psi(k) = psi0 + current_nu(3)*dt + dt*sum_{j=1}^{k-1} nu(3,j)  (k>=2)
    %
    % 姿态误差: psi(k) - ref_psi_k = psi_c_k + dt*sum_{j=1}^{k-1} nu(3,j)
    %   其中 psi_c_k = psi0 + current_nu(3)*dt - ref_psi_k  (常数)
    %
    % 代价 = w_psi * (psi_c_k + dt*sum_{j=1}^{k-1} nu(3,j))^2
    %
    % 在 0.5*x'Hx 形式:
    %   二次项: w_psi*dt^2*(sum_{j=1}^{k-1} nu(3,j))^2
    %           = w_psi*dt^2*sum_{i,j=1}^{k-1} nu(3,i)*nu(3,j)
    %           → H(nu_psi(i), nu_psi(j)) = 2*w_psi*dt^2  (对所有 i,j 对)
    %
    %   一次项: 2*w_psi*dt*psi_c_k * sum_{j=1}^{k-1} nu(3,j)
    %           → g(nu_psi(j)) += 2*w_psi*dt*psi_c_k  (j=1..k-1)
    %
    % 注意: k=1 时没有nu变量贡献(psi(1)是常数)，但仍产生常数代价

    for k = 1:K
        ref_idx = min(size(path, 2), step + k);
        ref_psi = path(3, ref_idx);

        psi_c_k = psi0 + current_nu(3)*dt - ref_psi;

        for i = 1:k-1
            for j = 1:k-1
                % [Bug2+Bug4修复] 加入所有跨阶段交叉项
                nu_psi_i = nu_start + (i-1)*3 + 2;  % [Bug1修复] 偏移2，不是3
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

    % 1.3 控制正则化: w_control * ||u||^2
    % ||u||^2 在 0.5*x'Hx 形式: H 对角 = 2*w_control
    for k = 1:K
        for i = 1:3
            u_idx = (k-1)*3 + i;
            H_mat(u_idx, u_idx) = H_mat(u_idx, u_idx) + 2 * w_control;
        end
    end

    % 1.4 RSS约束松弛: rho * ||u - u_hat||^2
    % = rho*(||u||^2 - 2*u'*u_hat + ||u_hat||^2)
    % 在 0.5*x'Hx+g'x 形式: H += 2*rho*I_u, g -= 2*rho*u_hat(:)
    for k = 1:K
        for i = 1:3
            u_idx = (k-1)*3 + i;
            H_mat(u_idx, u_idx) = H_mat(u_idx, u_idx) + 2 * rho;
            g_vec(u_idx) = g_vec(u_idx) - 2 * rho * u_hat(i, k);
        end
    end

    % 确保 H 对称
    H_mat = 0.5 * (H_mat + H_mat');

    %% =====================================================
    % 1.5 目标函数常数项 objective_constant
    % ======================================================
    % 完整 CVX 目标: J = 0.5*x'*H*x + g'*x + c
    % 其中 c 与决策变量无关，不影响最优解，但影响 objective value 的数值比较。
    %
    % c 包含三部分:
    %   (a) 位置跟踪常数: w_pos * ||c_k||^2  (k=2..K)
    %   (b) 姿态跟踪常数: w_psi * psi_c_k^2   (k=1..K)
    %   (c) RSS 正则常数: rho * ||u_hat||^2

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

    % (c) RSS 正则常数: rho * ||u_hat(:)||^2
    objective_constant = objective_constant + rho * sum(u_hat(:).^2);

    %% =====================================================
    % 2. 等式约束 (动力学递推)
    % ======================================================
    % nu(:,1) = current_nu + u(:,1)  =>  -u(:,1) + nu(:,1) = current_nu
    % nu(:,k+1) = nu(:,k) + u(:,k+1) => -u(:,k+1) - nu(:,k) + nu(:,k+1) = 0

    n_eq = 3*K;  % 3 for initial + 3*(K-1) for recursion
    A_eq = zeros(n_eq, n_var);
    b_eq = zeros(n_eq, 1);

    eq_row = 1;

    % 2.1 初始条件: nu(:,1) - u(:,1) = current_nu
    for i = 1:3
        u_idx = i;
        % [已确认正确] nu(i,1) = nu_start + (i-1)
        nu_idx = nu_start + i - 1;

        A_eq(eq_row, u_idx) = -1;
        A_eq(eq_row, nu_idx) = 1;
        b_eq(eq_row) = current_nu(i);

        eq_row = eq_row + 1;
    end

    % 2.2 递推约束: nu(:,k+1) - nu(:,k) - u(:,k+1) = 0
    for k = 1:K-1
        for i = 1:3
            u_idx = k*3 + i;
            % [Bug5修复] nu(i,k) = nu_start + (k-1)*3 + (i-1), 不是 +(i)
            nu_k_idx = nu_start + (k-1)*3 + i - 1;
            nu_kp1_idx = nu_start + k*3 + i - 1;

            A_eq(eq_row, u_idx) = -1;
            A_eq(eq_row, nu_k_idx) = -1;
            A_eq(eq_row, nu_kp1_idx) = 1;
            b_eq(eq_row) = 0;

            eq_row = eq_row + 1;
        end
    end

    %% =====================================================
    % 3. 线性不等式约束
    % ======================================================
    % (暂无纯线性不等式约束，所有非线性约束都在二次约束中)

    C_ineq = [];
    d_ineq = [];

    %% =====================================================
    % 4. 二次不等式约束 (QCQP)
    % ======================================================

    Hq_list = {};
    gq_list = {};
    uq_list = [];

    % 构造 nu_hat 序列 (用于转向锥线性化)
    nu_hat = zeros(3, K);
    nu_hat(:, 1) = current_nu + u_hat(:, 1);
    for k = 1:K-1
        nu_hat(:, k+1) = nu_hat(:, k) + u_hat(:, k+1);
    end

    % -------------------------------------------------------
    % 4.1 轮速约束: ||H{n}*nu(:,k)||^2 <= vimax^2
    % 等价于: nu(:,k)' * M{n} * nu(:,k) <= vimax^2
    % QCQP形式: 0.5*x'*Hq*x + gq'*x <= uq
    %   Hq 在 nu(:,k) 块有 2*M{n}, 其余为0
    %   gq = 0
    %   uq = vimax^2
    % -------------------------------------------------------

    for k = 1:K
        for n = 1:num_wheels
            Hq_k = zeros(n_var, n_var);
            gq_k = zeros(n_var, 1);

            % nu(:,k) 的索引范围 (已确认正确: block start = nu_start+(k-1)*3)
            nu_k_start = nu_start + (k-1)*3;
            nu_k_end = nu_start + (k-1)*3 + 2;

            Hq_k(nu_k_start:nu_k_end, nu_k_start:nu_k_end) = 2 * M_cell{n};

            Hq_list{end+1} = Hq_k;
            gq_list{end+1} = gq_k;
            uq_list(end+1) = vimax^2;
        end
    end

    % -------------------------------------------------------
    % 4.2 第一组转向锥凸化约束
    % R = [sin(dtheta), -cos(dtheta); cos(dtheta), sin(dtheta)]
    % -------------------------------------------------------

    R1 = [
        sin(delta_theta), -cos(delta_theta);
        cos(delta_theta),  sin(delta_theta)
    ];

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
    % 4.3 第二组转向锥凸化约束
    % R = [sin(dtheta), cos(dtheta); -cos(dtheta), sin(dtheta)]
    % -------------------------------------------------------

    R2 = [
        sin(delta_theta),  cos(delta_theta);
       -cos(delta_theta),  sin(delta_theta)
    ];

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
    % 5. 返回QCQP问题结构
    % ======================================================

    qp_problem.H = H_mat;
    qp_problem.g = g_vec;
    qp_problem.A = A_eq;
    qp_problem.b = b_eq;
    qp_problem.C = C_ineq;
    qp_problem.d = d_ineq;
    qp_problem.lb = [];
    qp_problem.ub = [];

    % QCQP约束
    qp_problem.Hq = Hq_list;
    qp_problem.gq = gq_list;
    qp_problem.uq = uq_list;

    % 元数据
    qp_problem.n_var = n_var;
    qp_problem.K = K;
    qp_problem.n_eq = n_eq;
    qp_problem.n_qcqp = length(uq_list);
    qp_problem.objective_constant = objective_constant;

end


%% =========================================================
% 构造转向锥凸化约束
% ==========================================================

function result = add_steering_cone_constraints( ...
    Hq_list, gq_list, uq_list, ...
    R, K, num_wheels, H_cell, M_cell, ...
    u_hat, nu_hat, current_nu, ...
    n_var, nu_start, delta_theta ...
)
% 对每组R，构造K*num_wheels个二次不等式约束
%
% 约束形式 (k>1):
%   ||H{n}*nu(:,k-1)||^2 + ||H{n}*(nu(:,k-1)+u(:,k))||^2
%     - ||lv||^2 - LH <= 0
%
% 其中:
%   lv = (I+R)*H{n}*nu_hat(:,k-1) + R*H{n}*u_hat(:,k)  [常数]
%   LH = sum_{l=1}^{k-1} 2*lv'*(I+R)*H{n}*(u(:,l)-u_hat(:,l))
%         + 2*lv'*R*H{n}*(u(:,k)-u_hat(:,k))  [u的线性函数]
%
% k=1时:
%   ||H{n}*current_nu||^2 + ||H{n}*(current_nu+u(:,1))||^2
%     - ||lv||^2 - 2*lv'*R*H{n}*(u(:,1)-u_hat(:,1)) <= 0

    for k = 1:K
        for n = 1:num_wheels
            Hq_k = zeros(n_var, n_var);
            gq_k = zeros(n_var, 1);

            Hn = H_cell{n};
            Mn = M_cell{n};

            if k > 1
                % lv = (I+R)*H{n}*nu_hat(:,k-1) + R*H{n}*u_hat(:,k)
                lv = (eye(2) + R) * Hn * nu_hat(:, k-1) + R * Hn * u_hat(:, k);

                % 二次项:
                % ||H{n}*nu(:,k-1)||^2 + ||H{n}*(nu(:,k-1)+u(:,k))||^2
                % = nu_k1'*Mn*nu_k1 + (nu_k1+u_k)'*Mn*(nu_k1+u_k)
                % = 2*nu_k1'*Mn*nu_k1 + 2*nu_k1'*Mn*u_k + u_k'*Mn*u_k
                %
                % 在 0.5*x'*Hq*x 形式中:
                % Hq(nu_k1,nu_k1) = 4*Mn
                % Hq(nu_k1,u_k) = 2*Mn
                % Hq(u_k,nu_k1) = 2*Mn
                % Hq(u_k,u_k) = 2*Mn

                % nu(:,k-1) block: start = nu_start+(k-2)*3, end = nu_start+(k-2)*3+2
                nu_k1_start = nu_start + (k-2)*3;
                nu_k1_end = nu_start + (k-2)*3 + 2;
                u_k_start = (k-1)*3 + 1;
                u_k_end = k*3;

                Hq_k(nu_k1_start:nu_k1_end, nu_k1_start:nu_k1_end) = 4 * Mn;
                Hq_k(nu_k1_start:nu_k1_end, u_k_start:u_k_end) = 2 * Mn;
                Hq_k(u_k_start:u_k_end, nu_k1_start:nu_k1_end) = 2 * Mn;
                Hq_k(u_k_start:u_k_end, u_k_start:u_k_end) = 2 * Mn;

                % 一次项 (来自 -LH):
                % u的线性系数:
                % 对 u(:,l), l=1..k-1: -2*lv'*(I+R)*H{n}
                % 对 u(:,k): -2*lv'*R*H{n}

                for l = 1:k-1
                    coeff = -2 * ((eye(2) + R)' * lv)' * Hn;  % 1x3
                    u_l_start = (l-1)*3 + 1;
                    gq_k(u_l_start:u_l_start+2) = gq_k(u_l_start:u_l_start+2) + coeff';
                end
                coeff_k = -2 * (R' * lv)' * Hn;  % 1x3
                gq_k(u_k_start:u_k_start+2) = gq_k(u_k_start:u_k_start+2) + coeff_k';

                % 常数项 (RHS):
                % uq = ||lv||^2 - [来自-LH的u_hat常数项]
                u_hat_const = 0;
                for l = 1:k-1
                    u_hat_const = u_hat_const + 2 * lv' * (eye(2) + R) * Hn * u_hat(:, l);
                end
                u_hat_const = u_hat_const + 2 * lv' * R * Hn * u_hat(:, k);

                uq_val = lv' * lv - u_hat_const;

            else
                % k = 1
                % lv = (I+R)*H{n}*current_nu + R*H{n}*u_hat(:,1)
                lv = (eye(2) + R) * Hn * current_nu + R * Hn * u_hat(:, 1);

                % 二次项:
                % ||H{n}*(current_nu+u(:,1))||^2 中决策变量部分: u_1'*Mn*u_1
                % 0.5*x'*Hq*x 形式: Hq(u_1,u_1) = 2*Mn

                u_1_start = 1;
                u_1_end = 3;

                Hq_k(u_1_start:u_1_end, u_1_start:u_1_end) = 2 * Mn;

                % 一次项:
                % 2*current_nu'*Mn*u_1 - 2*lv'*R*H{n}*u_1

                gq_k(u_1_start:u_1_end) = 2 * Mn' * current_nu - 2 * Hn' * R' * lv;

                % 常数项:
                % 约束: [0.5*x'Hq*x + gq'x <= uq]
                % 即所有常数项移到 uq

                uq_val = -2 * current_nu' * Mn * current_nu + lv' * lv - 2 * lv' * R * Hn * u_hat(:, 1);
            end

            % 确保 Hq 对称
            Hq_k = 0.5 * (Hq_k + Hq_k');

            Hq_list{end+1} = Hq_k;
            gq_list{end+1} = gq_k;
            uq_list(end+1) = uq_val;
        end
    end

    result = {Hq_list, gq_list, uq_list};
end
