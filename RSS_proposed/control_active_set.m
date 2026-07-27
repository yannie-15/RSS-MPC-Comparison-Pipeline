function [u, new_state_dot, velocity, solver_info] = control_active_set( ...
    path, step, state_dot, state, config ...
)
% CONTROL_ACTIVE_SET
% 论文对比算法: fmincon + Algorithm='active-set'
%
% 直接求解原始非凸优化问题 P_K (论文式14)
% 论文 Section IV.C: 此算法通常会失败（Infeasible），
% 因为 P_K 是非凸的且 active-set 需要可行初始点。
%
% 使用 build_baseline_pk 构造论文标准问题

    if nargin < 5 || isempty(config)
        config = defaultConfig();
    end

    % 构造论文标准 P^K
    pk = build_baseline_pk(path, step, state_dot, state, config);

    % fmincon 选项: active-set
    fmincon_opts = optimset( ...
        'Algorithm', 'active-set', ...
        'Display', 'off', ...
        'MaxIter', 200, ...
        'MaxFunEvals', 10000);

    % 求解
    global solver_time_array;
    if isempty(solver_time_array), solver_time_array = []; end

    solver_tic = tic;
    solver_info = struct('exitflag', 0, 'iterations', 0, 'finite', true);

    try
        % 抑制 active-set QP 子问题秩亏警告 (rank-deficient 是非凸问题正常现象)
        wstate = warning('off', 'optim:fmincon:RankDeficientMQ');
        wstate2 = warning('off', 'MATLAB:nearlySingularMatrix');
        cleanup = onCleanup(@() [warning(wstate); warning(wstate2)]);

        [x_sol, ~, exitflag, output] = fmincon( ...
            pk.objective, pk.x0_static, [], [], pk.Aeq, pk.beq, ...
            [], [], pk.nonlcon, fmincon_opts);

        solve_time = toc(solver_tic);
        solver_time_array(step) = solve_time;

        solver_info.exitflag = exitflag;
        solver_info.iterations = output.iterations;

        % 检查解是否有限
        if ~all(isfinite(x_sol))
            solver_info.finite = false;
            u = reshape(pk.x0_static(1:3*pk.K), 3, pk.K);
        else
            u = reshape(x_sol(1:3*pk.K), 3, pk.K);
        end

    catch ME
        solve_time = toc(solver_tic);
        solver_time_array(step) = solve_time;
        solver_info.exitflag = -999;
        solver_info.error = ME.message;
        solver_info.finite = false;
        u = zeros(3, pk.K);
    end

    % 输出
    velocity = pk.current_nu + u(:, 1);
    rotation_world_from_body = [cos(pk.psi0), -sin(pk.psi0), 0;
                                sin(pk.psi0),  cos(pk.psi0), 0;
                                0,             0,             1];
    new_state_dot = rotation_world_from_body * velocity;
end
