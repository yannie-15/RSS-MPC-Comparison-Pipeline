function summary = run_one_case(config, scenario)
% RUN_ONE_CASE 运行单次闭环仿真并计算完整指标
%
% 算法调用架构:
%   - proposed-3iter  → third_party/RSS_proposed/control_RSS.m  (git submodule)
%   - e-lmpc          → third_party/RSS_sqp/control_RSS.m       (git submodule)
%   - interior-point  → third_party/RSS_fmincon/control_RSS.m   (git submodule)
%   - active-set      → algorithms/control_active_set.m          (本地实现)
%
% 三个 submodule 的接口各不相同, 本函数负责适配:
%   RSS_proposed:  [new_state_dot] = control_RSS(path, k, state_dot, state)
%   RSS_sqp:       [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state)
%   RSS_fmincon:   [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state, params)
%   本地 active-set: [u, worldVelocity, bodyVelocity] = control_active_set(path, k, state_dot, state, config)
%
% 注: RSS_proposed 只输出 new_state_dot, 缺失的 u/solve_time 记为 NaN (用户已确认接受).

    if nargin < 1 || isempty(config)
        config = defaultConfig();
    end
    if nargin < 2 || isempty(scenario)
        scenario = struct();
        scenario.name = 'paper_fixed';
    end

    if ~isfield(scenario, 'name'), scenario.name = 'unnamed'; end
    if ~isfield(scenario, 'id'), scenario.id = 0; end

    if ~isfield(config, 'algorithm') || isempty(config.algorithm)
        config.algorithm = 'proposed-3iter';
    end

    algorithm = lower(config.algorithm);

    % 定位 submodule 路径
    script_dir = fileparts(mfilename('fullpath'));
    workspace_root = fileparts(script_dir);
    submodule_paths = struct( ...
        'proposed-3iter', fullfile(workspace_root, 'third_party', 'RSS_proposed'), ...
        'e-lmpc',         fullfile(workspace_root, 'third_party', 'RSS_sqp'), ...
        'interior-point', fullfile(workspace_root, 'third_party', 'RSS_fmincon') ...
    );

    %% =====================================================
    % 从 config / scenario 提取仿真参数
    % ======================================================
    num_steps = config.num_steps;
    num_wheels = size(config.wheel_pos, 1);

    % 生成参考轨迹
    path = generateReference(config, config.num_path_pts);

    % 初始状态
    if isfield(scenario, 'initialState') && ~isempty(scenario.initialState)
        state = scenario.initialState(:);
    else
        state = [0.05; 0.1; 0.2];
    end

    % 初始车体速度
    if isfield(scenario, 'initialVelocity') && ~isempty(scenario.initialVelocity)
        lastBodyVelocity = scenario.initialVelocity(:);
    else
        lastBodyVelocity = [0.01; 0.01; 0.01];
    end

    %% =====================================================
    % 预分配数组
    % ======================================================
    states = zeros(3, num_steps + 1);
    worldVelocities = zeros(3, num_steps);
    bodyVelocities = zeros(3, num_steps);
    wheelSpeeds = zeros(num_wheels, num_steps);
    wheelAngles = zeros(num_wheels, num_steps);
    executedU = zeros(3, num_steps);
    solverTimes_per_mpc_step = zeros(1, num_steps);

    success = true;
    failureReason = '';
    solved_count = 0;

    states(:, 1) = state;

    %% =====================================================
    % 闭环仿真主循环
    % ======================================================
    for k = 1:num_steps
        try
            % 按算法名分发调用
            switch algorithm
                case 'proposed-3iter'
                    % RSS_proposed: [new_state_dot] = control_RSS(path, k, state_dot, state)
                    addpath(submodule_paths.('proposed-3iter'));
                    worldVelocity = control_RSS(path, k, lastBodyVelocity, state');
                    rmpath(submodule_paths.('proposed-3iter'));
                    % 反推车体速度 (有 0.98 衰减, 近似)
                    R_bw = [cos(state(3)), sin(state(3)), 0;
                           -sin(state(3)), cos(state(3)), 0;
                            0,             0,             1];
                    bodyVelocity = R_bw * worldVelocity;
                    % u 无法精确获取 (0.98 衰减), 记为 NaN
                    u = NaN(3, 1);
                    solve_time = NaN;

                case 'e-lmpc'
                    % RSS_sqp: [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state)
                    addpath(submodule_paths.('e-lmpc'));
                    [worldVelocity, bodyVelocity, solve_time, ~] = ...
                        control_RSS(path, k, lastBodyVelocity, state');
                    rmpath(submodule_paths.('e-lmpc'));
                    % 反推 u(:,1) = bodyVelocity - lastBodyVelocity
                    u = bodyVelocity - lastBodyVelocity;

                case 'interior-point'
                    % RSS_fmincon: [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state, params)
                    addpath(submodule_paths.('interior-point'));
                    [worldVelocity, bodyVelocity, solve_time, ~] = ...
                        control_RSS(path, k, lastBodyVelocity, state', config);
                    rmpath(submodule_paths.('interior-point'));
                    u = bodyVelocity - lastBodyVelocity;

                case 'active-set'
                    % 本地: [u, worldVelocity, bodyVelocity] = control_active_set(path, k, state_dot, state, config)
                    [u_full, worldVelocity, bodyVelocity] = ...
                        control_active_set(path, k, lastBodyVelocity, state', config);
                    u = u_full(:, 1);
                    solve_time = NaN;  % 本地版本不输出 solve_time

                otherwise
                    error('run_one_case:UnknownAlgorithm', ...
                        'Unknown algorithm: %s', config.algorithm);
            end

            [wheelSpeed, wheelAngle] = computeWheelOutputs(bodyVelocity, config);

            executedU(:, k) = u;
            states(:, k) = state;
            worldVelocities(:, k) = worldVelocity;
            bodyVelocities(:, k) = bodyVelocity;
            wheelSpeeds(:, k) = wheelSpeed;
            wheelAngles(:, k) = wheelAngle;
            if ~isnan(solve_time)
                solverTimes_per_mpc_step(k) = solve_time;
            end
            solved_count = solved_count + 1;

            state = propagateState(state, worldVelocity, config);
            lastBodyVelocity = bodyVelocity;
            states(:, k+1) = state;

        catch ME
            success = false;
            failureReason = getReport(ME);
            break;
        end
    end

    %% =====================================================
    % 截断到实际完成的步数
    % ======================================================
    if solved_count == 0
        states = states(:, 1);
        worldVelocities = zeros(3, 0);
        bodyVelocities = zeros(3, 0);
        executedU = zeros(3, 0);
        wheelSpeeds = zeros(num_wheels, 0);
        wheelAngles = zeros(num_wheels, 0);
        solverTimes_per_mpc_step = [];
    else
        states = states(:, 1:solved_count+1);
        worldVelocities = worldVelocities(:, 1:solved_count);
        bodyVelocities = bodyVelocities(:, 1:solved_count);
        executedU = executedU(:, 1:solved_count);
        wheelSpeeds = wheelSpeeds(:, 1:solved_count);
        wheelAngles = wheelAngles(:, 1:solved_count);
        solverTimes_per_mpc_step = solverTimes_per_mpc_step(1:solved_count);
    end

    %% =====================================================
    % 计算指标
    % ======================================================
    metrics = computeMetrics(path, states, config, ...
        solverTimes_per_mpc_step, wheelSpeeds, wheelAngles, executedU);
    metrics.successRate = solved_count / num_steps;

    % 求解时间统计 (proposed 全为 NaN 时, 这些也是 NaN)
    valid_times = solverTimes_per_mpc_step(~isnan(solverTimes_per_mpc_step));
    if ~isempty(valid_times)
        metrics.meanSolveTime = mean(valid_times);
        metrics.maxSolveTime = max(valid_times);
        metrics.totalSolveTime = sum(valid_times);
    else
        metrics.meanSolveTime = NaN;
        metrics.maxSolveTime = NaN;
        metrics.totalSolveTime = NaN;
    end

    %% =====================================================
    % 组装输出
    % ======================================================
    summary = struct();
    summary.success = success;
    summary.failureReason = failureReason;
    summary.scenario = scenario;
    summary.config = config;
    summary.algorithm = algorithm;
    summary.states = states;
    summary.worldVelocities = worldVelocities;
    summary.bodyVelocities = bodyVelocities;
    summary.executedU = executedU;
    summary.wheelSpeeds = wheelSpeeds;
    summary.wheelAngles = wheelAngles;
    summary.solverTimes = solverTimes_per_mpc_step;
    summary.path = path;
    summary.metrics = metrics;
    summary.numWheels = num_wheels;

    if success
        fprintf('[Case %d: %s] RMSE=%.6f, meanSolveTime=%.4fs, J_total=%.4f, successRate=%.0f%%\n', ...
            scenario.id, scenario.name, metrics.rmse, metrics.meanSolveTime, ...
            metrics.trajectoryCost, metrics.successRate*100);
    else
        fprintf('[Case %d: %s] FAILED at step %d/%d: %s\n', ...
            scenario.id, scenario.name, solved_count, num_steps, failureReason);
    end
end
