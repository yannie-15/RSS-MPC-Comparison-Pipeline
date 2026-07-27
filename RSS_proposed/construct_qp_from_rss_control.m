function qp_problem = construct_qp_from_rss_control(path, step, current_nu, state, u_hat, params)
% CONSTRUCT_QP_FROM_RSS_CONTROL
% 从RSS控制问题显式构造QP问题矩阵
%
% 问题形式（标准QP）：
%   min  0.5*[u; nu]'*H*[u; nu] + g'*[u; nu]
%   s.t. A_eq*[u; nu] = b_eq       (动力学约束)
%        A_ineq*[u; nu] <= b_ineq (轮速约束)
%        
% 输入：
%   path       : 3×N 参考路径
%   step       : 当前MPC步
%   current_nu : 3×1 当前车体坐标系速度
%   state      : 1×3 当前世界坐标系状态 [x, y, psi]
%   u_hat      : 3×K 上一次迭代的控制序列
%   params     : 配置参数
%
% 输出：
%   qp_problem : 结构体，包含H, g, A_eq, b_eq, A_ineq, b_ineq

    %% =====================================================
    % 提取参数
    % ======================================================
    
    K = 6;  % 预测步数
    dt = params.dt;
    phidotmax = params.phidotmax;
    vimax = params.vimax;
    wheel_pos = params.wheel_pos;
    num_wheels = size(wheel_pos, 1);
    
    % 当前状态
    current_xy = [state(1); state(2)];
    psi0 = state(3);
    R_psi0 = [cos(psi0), -sin(psi0); sin(psi0), cos(psi0)];
    
    % 代价函数权重
    w_pos = 30;      % 位置跟踪权重
    w_psi = 0.15;    % 姿态跟踪权重
    w_control = 0.3; % 控制正则化权重
    rho = 0.01;      % RSS约束松弛权重
    
    %% =====================================================
    % 决策变量编排
    % ======================================================
    % 变量顺序：[u(:); nu(:)]
    % u的索引：1:3K
    % nu的索引：3K+1:6K
    
    n_var = 6 * K;
    u_start = 1;
    u_end = 3 * K;
    nu_start = 3 * K + 1;
    nu_end = 6 * K;
    
    %% =====================================================
    % 1. 构造Hessian矩阵 H 和一次项 g
    % ======================================================
    
    H = zeros(n_var, n_var);
    g = zeros(n_var, 1);
    
    % 由于代价函数涉及复杂的非线性（位置和姿态），
    % 这里我们只能构造二次部分
    % 完整的二次规划化需要在RSS迭代的框架内进行
    
    % 控制正则化项：0.3*sum(u^2)
    for k = 1:K
        for i = 1:3
            u_idx = (k-1)*3 + i;
            H(u_idx, u_idx) = H(u_idx, u_idx) + w_control;
        end
    end
    
    % RSS约束松弛项：rho*sum((u - u_hat)^2)
    % = rho*sum(u^2) - 2*rho*u'*u_hat + rho*sum(u_hat^2)
    for k = 1:K
        for i = 1:3
            u_idx = (k-1)*3 + i;
            H(u_idx, u_idx) = H(u_idx, u_idx) + rho;
            g(u_idx) = g(u_idx) - 2*rho*u_hat(i, k);
        end
    end
    
    %% =====================================================
    % 2. 构造等式约束 (动力学递推)
    % ======================================================
    
    % nu(:, 1) = current_nu + u(:, 1)
    % => u(:, 1) - nu(:, 1) = -current_nu
    % 或者: u(:, 1) + (-1)*nu(:, 1) = -current_nu
    
    % nu(:, k+1) = nu(:, k) + u(:, k+1), for k=1:K-1
    % => u(:, k+1) - nu(:, k) + nu(:, k+1) = 0
    % => -nu(:, k) + u(:, k+1) + nu(:, k+1) = 0
    
    n_eq_constraints = 3 + 3*(K-1);  % 第一个等式 + (K-1)个动力学递推
    A_eq = zeros(n_eq_constraints, n_var);
    b_eq = zeros(n_eq_constraints, 1);
    
    eq_idx = 1;
    
    % 第一个等式：nu(:, 1) = current_nu + u(:, 1)
    for i = 1:3
        u_idx = (0)*3 + i;  % u(:, 1)的索引
        nu_idx = (0)*3 + i; % nu(:, 1)的索引
        A_eq(eq_idx, nu_start - 1 + nu_idx) = 1;
        A_eq(eq_idx, u_start - 1 + u_idx) = -1;
        b_eq(eq_idx) = -current_nu(i);
        eq_idx = eq_idx + 1;
    end
    
    % 递推约束：nu(:, k+1) = nu(:, k) + u(:, k+1)
    for k = 1:K-1
        for i = 1:3
            u_idx = k*3 + i;              % u(:, k+1)的索引
            nu_k_idx = (k-1)*3 + i;       % nu(:, k)的索引
            nu_kp1_idx = k*3 + i;         % nu(:, k+1)的索引
            
            A_eq(eq_idx, nu_start - 1 + nu_k_idx) = -1;
            A_eq(eq_idx, u_start - 1 + u_idx) = -1;
            A_eq(eq_idx, nu_start - 1 + nu_kp1_idx) = 1;
            b_eq(eq_idx) = 0;
            eq_idx = eq_idx + 1;
        end
    end
    
    %% =====================================================
    % 3. 构造不等式约束 (轮速限制)
    % ======================================================
    
    % 轮速约束：||H{n}*nu(:, k)|| <= vimax
    % 这是二次约束，需要特殊处理
    
    % 计算轮矩阵H
    H_wheel = cell(1, num_wheels);
    for n = 1:num_wheels
        H_wheel{n} = [
            1, 0, -wheel_pos(n, 2);
            0, 1,  wheel_pos(n, 1)
        ];
    end
    
    % 对于每个轮子和每个时步，生成约束：
    % (H{n}*nu)^2 <= vimax^2
    % => (H{n}*nu)_x^2 + (H{n}*nu)_y^2 <= vimax^2
    
    n_ineq_wheel_constraints = K * num_wheels * 2;
    A_ineq_wheel = zeros(n_ineq_wheel_constraints, n_var);
    b_ineq_wheel = zeros(n_ineq_wheel_constraints, 1);
    
    % 线性近似（对每个轮子的每个速度分量）
    ineq_idx = 1;
    for k = 1:K
        for n = 1:num_wheels
            % Hn*nu_k的两个分量
            for comp = 1:2
                % (H{n}(comp, :) * nu(:, k))^2 <= vimax^2
                % 线性近似：2*(ref_val)*(变量val) - (ref_val)^2 <= vimax^2
                
                nu_k_indices = nu_start - 1 + (k-1)*3 + (1:3);
                
                A_ineq_wheel(ineq_idx, nu_k_indices) = 2 * H_wheel{n}(comp, :);
                
                % 右侧常数项需要在循环外处理（因为涉及参考值）
                % 这里先设置为vimax^2
                b_ineq_wheel(ineq_idx) = vimax^2;
                
                ineq_idx = ineq_idx + 1;
            end
        end
    end
    
    % 合并不等式约束
    A_ineq = A_ineq_wheel;
    b_ineq = b_ineq_wheel;
    
    %% =====================================================
    % 4. 构造输出结构体
    % ======================================================
    
    qp_problem.H = H;
    qp_problem.g = g;
    qp_problem.A = A_eq;
    qp_problem.b = b_eq;
    qp_problem.C = A_ineq;
    qp_problem.d = b_ineq;
    qp_problem.lb = [];
    qp_problem.ub = [];
    
    % 额外信息
    qp_problem.n_var = n_var;
    qp_problem.K = K;
    qp_problem.n_u = 3*K;
    qp_problem.n_nu = 3*K;
    qp_problem.current_nu = current_nu;
    qp_problem.state = state;
    qp_problem.u_hat = u_hat;
    
end
