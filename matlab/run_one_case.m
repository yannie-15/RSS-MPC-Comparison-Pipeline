function summary = run_one_case(config, scenario)
% RUN_ONE_CASE 运行单次闭环仿真并计算完整指标
%
% 算法调用架构:
%   - proposed-3iter  → algorithms/RSS_proposed/control_RSS.m  (git submodule)
%   - e-lmpc          → algorithms/RSS_sqp/control_RSS.m       (git submodule)
%   - interior-point  → algorithms/RSS_fmincon/control_RSS.m   (git submodule)
%   - active-set      → algorithms/control_active_set.m        (本地实现)
%
% 三个 submodule 的接口各不相同, 本函数负责适配:
%   RSS_proposed:  [new_state_dot] = control_RSS(path, k, state_dot, state)
%   RSS_sqp:       [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state)
%   RSS_fmincon:   [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state, params)
%   本地 active-set: [u, worldVelocity, bodyVelocity] = control_active_set(path, k, state_dot, state, config)
%
% Config 覆盖机制:
%   proposed / e-lmpc 的 control_RSS 不接收外部 config, 内部调 config()。
%   本函数通过 setup_config_override 在临时目录写 config.m, 返回 main 传入的
%   随机 config (含 seed 场景参数), 通过 addpath 覆盖 submodule 的 config.m。
%   interior-point / active-set 接收外部 config, 直接传入。
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
    algorithms_dir = fullfile(workspace_root, 'algorithms');
    submodule_dirs = containers.Map( ...
        {'proposed-3iter', 'e-lmpc', 'interior-point'}, ...
        {fullfile(algorithms_dir, 'RSS_proposed'), ...
         fullfile(algorithms_dir, 'RSS_sqp'), ...
         fullfile(algorithms_dir, 'RSS_fmincon')} ...
    );

    %% =====================================================
    % 设置 config 覆盖: 让 submodule 内部的 config() 返回 main 的随机 config
    % (proposed / e-lmpc 的 control_RSS 不接收外部 config, 内部调 config())
    % ======================================================
    override_dir = setup_config_override(config);
    cleanup_obj = onCleanup(@() rmpath(override_dir));

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
                    addpath(submodule_dirs('proposed-3iter'));
                    addpath(override_dir);  % 确保 main 的 config 覆盖 submodule 的 config
                    worldVelocity = control_RSS(path, k, lastBodyVelocity, state');
                    rmpath(submodule_dirs('proposed-3iter'));
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
                    addpath(submodule_dirs('e-lmpc'));
                    addpath(override_dir);  % 确保 main 的 config 覆盖 submodule 的 config
                    [worldVelocity, bodyVelocity, solve_time, ~] = ...
                        control_RSS(path, k, lastBodyVelocity, state');
                    rmpath(submodule_dirs('e-lmpc'));
                    % 反推 u(:,1) = bodyVelocity - lastBodyVelocity
                    u = bodyVelocity - lastBodyVelocity;

                case 'interior-point'
                    % RSS_fmincon: [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state, params)
                    addpath(submodule_dirs('interior-point'));
                    [worldVelocity, bodyVelocity, solve_time, ~] = ...
                        control_RSS(path, k, lastBodyVelocity, state', config);
                    rmpath(submodule_dirs('interior-point'));
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


%% =========================================================
% 辅助函数: 创建临时 config.m 覆盖 submodule 的 config
%% ==========================================================

function override_dir = setup_config_override(cfg)
%SETUP_CONFIG_OVERRIDE 创建临时 config.m, 用 main 的 config 覆盖 submodule 的 config
%
% proposed / e-lmpc 的 control_RSS.m 内部调 config() 取参数, 不接收外部 config。
% 本函数在临时目录写一个 config.m, 返回 main 传入的 cfg, 通过 addpath 覆盖。
%
% 调用前: addpath(override_dir) 确保 override 优先于 submodule 目录

    override_dir = fullfile(tempdir, 'rss_config_override');
    if ~exist(override_dir, 'dir'), mkdir(override_dir); end

    % 保存 cfg 到 mat (每次调用更新, 因为不同 seed 的 config 不同)
    cfg_file = fullfile(override_dir, 'config_data.mat');
    save(cfg_file, 'cfg', '-v7');

    % 写 config.m (只写一次, 复用)
    config_m = fullfile(override_dir, 'config.m');
    if ~exist(config_m, 'file')
        fid = fopen(config_m, 'w');
        fprintf(fid, 'function params = config()\n');
        fprintf(fid, '    f = fullfile(fileparts(mfilename(''fullpath'')), ''config_data.mat'');\n');
        fprintf(fid, '    loaded = load(f, ''cfg'');\n');
        fprintf(fid, '    params = loaded.cfg;\n');
        fprintf(fid, 'end\n');
        fclose(fid);
    end
end
