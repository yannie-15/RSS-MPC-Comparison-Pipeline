function [u, new_state_dot, velocity, solver_info] = control_eLMPC( ...
    path, step, state_dot, state, config ...
)
% CONTROL_ELMPC - aligned with RSS_sqp
%   Decision var: u (3K=18), ceq constraints, SQP MaxIter=1
%   delta_theta = dt*phidotmax, R1 = [sin,-cos;cos,sin]

    if nargin < 5 || isempty(config)
        config = defaultConfig();
    end

    K = 6; rho = 0.01; k1 = 1;

    current_xy = [state(1); state(2)];
    psi0 = state(3);
    current_nu = state_dot(:);
    R_psi0 = [cos(psi0), -sin(psi0); sin(psi0), cos(psi0)];

    H = cell(1, 4);
    for n = 1:4
        H{n} = [1, 0, -config.wheel_pos(n, 2);
                0, 1,  config.wheel_pos(n, 1)];
    end

    u_hat = zeros(3, K);
    x0 = reshape(u_hat, 3*K, 1);

    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'off', ...
        'ConstraintTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6, ...
        'MaxFunctionEvaluations', 1e5, ...
        'MaxIterations', 1);

    obj_fun = @(x) nmpc_objective(x, K, current_xy, path, step, psi0, R_psi0, ...
        current_nu, k1, rho, u_hat, config);
    nonlcon = @(x) nmpc_constraints(x, K, H, current_nu, u_hat, config);

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

    velocity = current_nu + u(:, 1);
    rot = [cos(psi0), -sin(psi0), 0;
           sin(psi0),  cos(psi0), 0;
           0,          0,          1];
    new_state_dot = rot * velocity;
end

function J = nmpc_objective(x, K, current_xy, path, step, psi0, R_psi0, ...
    current_nu, k1, rho, u_hat, params)

    u = reshape(x, 3, K);
    nu = zeros(3, K); NU = zeros(2, K); psi = zeros(1, K);
    NU(:,1) = current_nu(1:2) * params.dt;
    psi(1) = psi0 + current_nu(3) * params.dt;
    nu(:,1) = current_nu + u(:,1);
    for k = 1:K-1
        nu(:,k+1) = nu(:,k) + u(:,k+1);
        NU(:,k+1) = NU(:,k) + nu(1:2,k) * params.dt;
        psi(k+1) = psi(k) + nu(3,k) * params.dt;
    end
    J = 0;
    for k = 2:K
        idx = min(size(path,2), step+k);
        pos_err = current_xy - path(1:2,idx) + R_psi0 * NU(:,k);
        J = J + 30 * sum(pos_err.^2);
    end
    for k = 1:K
        idx = min(size(path,2), step+k);
        J = J + k1 * (psi(k) - path(3,idx))^2;
    end
    J = J + 0.3 * sum(u(:).^2);
end

function [c, ceq] = nmpc_constraints(x, K, H, current_nu, u_hat, params)

    u = reshape(x, 3, K);
    nu = zeros(3, K);
    c = []; ceq = [];
    delta_theta = params.dt * params.phidotmax;

    nu(:,1) = current_nu + u(:,1);
    ceq = [ceq; nu(:,1) - (current_nu + u(:,1))];
    for k = 1:K-1
        nu(:,k+1) = nu(:,k) + u(:,k+1);
        ceq = [ceq; nu(:,k+1) - (nu(:,k) + u(:,k+1))];
    end

    for k = 1:K
        for n = 1:4
            wheel_vel = H{n} * nu(:,k);
            c = [c; norm(wheel_vel,2) - params.vimax];
        end
    end

    R1 = [sin(delta_theta), -cos(delta_theta);
          cos(delta_theta),  sin(delta_theta)];
    for k = 1:K
        for n = 1:4
            if k == 1, nu_km1 = current_nu; else, nu_km1 = nu(:,k-1); end
            term = nu_km1' * H{n}' * R1 * H{n} * (nu_km1 + u(:,k));
            c = [c; -term];
        end
    end
    for k = 1:K
        for n = 1:4
            if k == 1, nu_km1 = current_nu; else, nu_km1 = nu(:,k-1); end
            term = nu_km1' * H{n}' * R1' * H{n} * (nu_km1 + u(:,k));
            c = [c; -term];
        end
    end
end
