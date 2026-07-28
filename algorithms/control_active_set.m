function [u, new_state_dot, velocity, solver_info] = control_active_set( ...
    path, step, state_dot, state, config ...
)
% CONTROL_ACTIVE_SET
% 论文对比算法: fmincon + Algorithm='active-set'
%
% 使用与 interior-point 相同的 P^K 问题定义 (对齐 RSS_fmincon):
%   - 决策变量仅 u (3K=18), nu 在目标/约束内部递推
%   - 等式约束用 ceq (非线性), 不用 Aeq/beq
%   - rho=0 (无正则项)
%   - 锥约束: 原始双线性形式
%     delta_theta = 0.5*pi - dt*phidotmax
%   - fmincon: Algorithm='active-set'
%
% 论文 Section IV.C: 此算法通常会失败（Infeasible），
% 因为 P_K 是非凸的且 active-set 需要可行初始点。

    if nargin < 5 || isempty(config)
        config = defaultConfig();
    end

    K = 6;
    rho = 0;      % 与 interior-point 一致: rho=0
    k1 = 1;

    current_xy = [state(1); state(2)];
    psi0 = state(3);
    current_nu = state_dot(:);
    R_psi0 = [cos(psi0), -sin(psi0); sin(psi0), cos(psi0)];

    % 车轮雅可比矩阵
    H = cell(1, 4);
    for n = 1:4
        H{n} = [1, 0, -config.wheel_pos(n, 2);
                0, 1,  config.wheel_pos(n, 1)];
    end

    % 初始点: u=0
    u_hat = zeros(3, K);
    x0 = u_hat(:);

    % fmincon 选项 (active-set)
    options = optimoptions('fmincon', ...
        'Algorithm', 'active-set', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 1e5, ...
        'MaxIterations', 2000, ...
        'OptimalityTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-6);

    % 构建目标和约束
    obj_fun = @(x) nmpc_objective(x, K, current_xy, path, step, psi0, R_psi0, ...
        current_nu, k1, rho, u_hat, config);
    nonlcon = @(x) nmpc_constraints(x, K, H, current_nu, u_hat, config);

    % 求解
    global solver_time_array;
    if isempty(solver_time_array), solver_time_array = []; end

    solver_tic = tic;
    solver_info = struct('exitflag', 0, 'iterations', 0, 'finite', true);

    try
        wstate = warning('off', 'optim:fmincon:RankDeficientMQ');
        wstate2 = warning('off', 'MATLAB:nearlySingularMatrix');
        cleanup = onCleanup(@() [warning(wstate); warning(wstate2)]);

        [x_sol, ~, exitflag, output] = fmincon( ...
            obj_fun, x0, [], [], [], [], [], [], nonlcon, options);

        solve_time = toc(solver_tic);
        solver_time_array(step) = solve_time;

        solver_info.exitflag = exitflag;
        solver_info.iterations = output.iterations;

        if ~all(isfinite(x_sol))
            solver_info.finite = false;
            u = zeros(3, K);
        else
            u = reshape(x_sol, 3, K);
        end

    catch ME
        solve_time = toc(solver_tic);
        solver_time_array(step) = solve_time;
        solver_info.exitflag = -999;
        solver_info.error = ME.message;
        solver_info.finite = false;
        u = zeros(3, K);
    end

    % 输出
    velocity = current_nu + u(:, 1);
    rotation_world_from_body = [cos(psi0), -sin(psi0), 0;
                                sin(psi0),  cos(psi0), 0;
                                0,          0,          1];
    new_state_dot = rotation_world_from_body * velocity;
end


%% =========================================================
% 目标函数 (对齐原始 RSS_fmincon, 与 interior-point 相同)
% ==========================================================
function J = nmpc_objective(x, K, current_xy, path, step, psi0, R_psi0, ...
    current_nu, k1, rho, u_hat, params)

    u = reshape(x, 3, K);
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

    % 代价
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


%% =========================================================
% 非线性约束 (对齐原始 RSS_fmincon, 与 interior-point 相同)
% ==========================================================
function [c, ceq] = nmpc_constraints(x, K, H, current_nu, u_hat, params)

    u = reshape(x, 3, K);
    nu = zeros(3, K);
    c = [];
    ceq = [];

    % 锥约束参数 (对齐原始 RSS_fmincon: delta_theta = 0.5*pi - dt*phidotmax)
    delta_theta = 0.5*pi - params.dt * params.phidotmax;

    % 动力学等式约束 (用 ceq)
    nu(:, 1) = current_nu + u(:, 1);
    ceq = [ceq; nu(:, 1) - (current_nu + u(:, 1))];
    for k = 1:K-1
        nu(:, k+1) = nu(:, k) + u(:, k+1);
        ceq = [ceq; nu(:, k+1) - (nu(:, k) + u(:, k+1))];
    end

    % 车轮速度约束: ||H_n * nu^k|| <= vimax
    for k = 1:K
        for n = 1:4
            wheel_vel = H{n} * nu(:, k);
            c = [c; norm(wheel_vel, 2) - params.vimax];
        end
    end

    % 转向锥约束 (原始 RSS_fmincon 双线性形式)
    % nu^{k-1}' * H' * R' * H * (nu^{k-1} + u^k) >= 0
    R1 = [cos(delta_theta), -sin(delta_theta);
          sin(delta_theta),  cos(delta_theta)];
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

    R2 = [cos(delta_theta),  sin(delta_theta);
         -sin(delta_theta),  cos(delta_theta)];
    for k = 1:K
        for n = 1:4
            if k == 1
                nu_km1 = current_nu;
            else
                nu_km1 = nu(:, k-1);
            end
            term = nu_km1' * H{n}' * R2' * H{n} * (nu_km1 + u(:,k));
            c = [c; -term];
        end
    end
end
