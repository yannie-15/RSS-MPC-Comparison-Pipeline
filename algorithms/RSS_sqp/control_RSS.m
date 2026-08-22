function [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state)
    
    params = config();
   
    % ================= Param Setup =================
    if isfield(params, 'K') && ~isempty(params.K)
        K = params.K;
    else
        K = 6;
    end
    % 正则化权重 rho: 优先取 pipeline 注入 (--rho 逐步序列经 config.rho),
    % 未注入时保持原默认 0.01
    if isfield(params, 'rho') && ~isempty(params.rho)
        rho = params.rho;
    else
        rho = 0.01;
    end
    k1 = 1; epsilon = 0;
    solve_time = 0;
    current_xy = [state(1), state(2)]';
    psi0 = state(3); current_nu = state_dot;
    R_psi0 = [cos(psi0), -sin(psi0); sin(psi0), cos(psi0)];

    % ================= H =================
    H = cell(1, 4);
    for n = 1:4
        H{n} = [1, 0, -params.wheel_pos(n, 2); 0, 1, params.wheel_pos(n, 1)];
    end

    % ================= 迭代 Setup =================
    max_iter = 3; u_hat = zeros(3,K);
    x0 = reshape(u_hat, 3*K, 1);

    % ================= 配置fmincon 有效集法 =================
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...      % 切换为有效集法
        'Display', 'off', ...               % 关闭打印
        'ConstraintTolerance', 1e-6, ...    % 有效集法适配的约束精度
        'OptimalityTolerance', 1e-6, ...    % 最优性精度
        'MaxFunctionEvaluations', 1e5, ...  % 最大函数评估次数
        'MaxIterations', 1);             % 最大迭代次数

    A = []; b = [];
    Aeq = []; beq = [];
    lb = []; ub = [];

    % ================= 构建目标 =================
    obj_fun = @(u) nmpc_objective(u, K, current_xy, path, step, psi0, R_psi0, ...
        current_nu, k1, rho, u_hat, params);
    nonlcon = @(u) nmpc_constraints(u, K, H, current_nu, u_hat, params);

    % ================= 迭代求解 =================
    t_start = tic;
    [u_opt, fval, exitflag, output] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    iter_time = toc(t_start);
    solve_time = solve_time + iter_time;
    iter_num = output.iterations;
    u_opt = reshape(u_opt, 3, K);
    
    % 打印收敛信息 (三算法统一口径):
    %   status =  1 : 收敛 (fmincon exitflag > 0; 1-5 为不同收敛判据)
    %   status =  0 : 达到 MaxIterations 停止 (fmincon exitflag = 0;
    %                 本算法为 "1 iteration edition" 设计, MaxIterations=1,
    %                 每步单次 SQP 迭代即返回, status=0 属预期而非失败)
    %   status = -1 : 求解失败 (fmincon exitflag < 0, 如无可行点)
    status = sign(exitflag);
    fprintf('exitflag=%d, status=%d\n', exitflag, status);
        
    u_hat = u_opt;
    x0 = reshape(u_hat, 3*K, 1);

    new_state_dot = [cos(state(3)), -sin(state(3)), 0;
                     sin(state(3)),  cos(state(3)), 0;
                                 0,              0, 1] * (state_dot + u_hat(:, 1));
    velocity = current_nu + u_hat(:, 1);
end

% ================= 目标函数 =================
function J = nmpc_objective(u, K, current_xy, path, step, psi0, R_psi0, ...
    current_nu, k1, rho, u_hat, params)
    u = reshape(u, 3, K);
    nu = zeros(3, K);
    NU = zeros(2, K);
    psi = zeros(1, K);

    NU(:, 1) = current_nu(1:2) * params.dt;
    psi(1) = psi0 + current_nu(3) * params.dt;
    nu(:, 1) = current_nu + u(:, 1);

    for k = 1:K-1
        nu(:, k+1) = nu(:, k) + u(:, k+1);
        NU(:, k+1) = NU(:, k) + nu(1:2, k) * params.dt;
        psi(k+1) = psi(k) + nu(3, k) * params.dt;
    end

% ================= 代价 =================
    J = 0;
    for k = 2:K
        idx = min(size(path, 2), step + k);
        pos_err = current_xy - path(1:2, idx) + R_psi0 * NU(:, k);
        J = J + 30 * sum(pos_err.^2);
    end

    for k = 1:K
        idx = min(size(path, 2), step + k);
        J = J + k1 * (psi(k) - path(3, idx))^2;
    end
    
    J = J + 0.3 * sum(u(:).^2);
end

% ================= 非线性约束 =================
function [c, ceq] = nmpc_constraints(u, K, H, current_nu, u_hat, params)
    u = reshape(u, 3, K);
    nu = zeros(3, K);
    c = [];
    ceq = [];
    delta_theta = params.dt * params.phidotmax;

    % ================= 动力学等式约束 =================
    nu(:, 1) = current_nu + u(:, 1);
    ceq = [ceq; nu(:, 1) - (current_nu + u(:, 1))];
    for k = 1:K-1
        nu(:, k+1) = nu(:, k) + u(:, k+1);
        ceq = [ceq; nu(:, k+1) - (nu(:, k) + u(:, k+1))];
    end

    % ================= 车轮速度凸约束 =================
    for k = 1:K
        for n = 1:4
            wheel_vel = H{n} * nu(:, k);
            c = [c; norm(wheel_vel, 2) - params.vimax];
        end
    end


    % 最原始约束
    % ================= v H R H (ν + u) ≥ 0 =================
    R1 = [sin(delta_theta), -cos(delta_theta);
          cos(delta_theta),  sin(delta_theta)];
    for k = 1:K 
        for n = 1:4
            if k == 1
                nu_km1 = current_nu;
            else
                nu_km1 = nu(:, k - 1);
            end
            term = nu_km1' * H{n}' * R1 * H{n} * (nu_km1 + u(:,k));
            c = [c; -term];
        end
    end

    
    for k = 1:K
        for n = 1:4
            if k == 1
                nu_km1 = current_nu;
            else
                nu_km1 = nu(:, k-1);
            end
            term = nu_km1' * H{n}' * R1' * H{n} * (nu_km1 + u(:,k));
            c = [c; -term];
        end
    end
end