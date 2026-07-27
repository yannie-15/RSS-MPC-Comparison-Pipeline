function summary = run_one_case(config, scenario)
% RUN_ONE_CASE 运行单次闭环仿真并计算完整指标
%
% 输入：
%   config   : 配置结构体（由 defaultConfig 或 scenario_generator 生成）
%   scenario : 场景描述结构体（可选，用于标记）
%
% 输出：
%   summary : 结构体，包含仿真结果和评价指标
%
% 修复记录:
%   - [P1] 保存 executedU 以计算式(21)轨迹代价
%   - [P2] 求解时间正确聚合：区分 per-convex-subproblem vs per-MPC-step
%   - [P2] 记录最终状态（控制执行后的状态也保存）
%   - [P2] 动态分配轮速/轮角数组（不再写死4轮）
%   - [P2] 记录每个MPC步的真实求解器backend

    if nargin < 1 || isempty(config)
        config = defaultConfig();
    end
    if nargin < 2 || isempty(scenario)
        scenario = struct();
        scenario.name = 'paper_fixed';
    end

    % 确保 scenario 有必要字段
    if ~isfield(scenario, 'name'), scenario.name = 'unnamed'; end
    if ~isfield(scenario, 'id'), scenario.id = 0; end

    % 根据 config.algorithm 选择控制器函数句柄
    % 支持: 'proposed-3iter', 'proposed-1iter', 'e-lmpc',
    %        'active-set', 'interior-point'
    % 默认: 'proposed-3iter' (即 control_RSS_v2 的标准行为)

    if ~isfield(config, 'algorithm') || isempty(config.algorithm)
        config.algorithm = 'proposed-3iter';
    end

    algorithm = lower(config.algorithm);

    switch algorithm
        case 'proposed-3iter'
            controller_fn = @control_RSS_v2;
            max_outer_iter = 3;    % 3次RSS迭代
        case 'proposed-1iter'
            controller_fn = @control_RSS_v2;
            max_outer_iter = 1;    % 1次RSS迭代
            config.max_iter_override = 1;  % 传递给控制器
        case 'e-lmpc'
            controller_fn = @control_eLMPC;
            max_outer_iter = 1;    % e-LMPC: 1次SQP迭代
        case 'active-set'
            controller_fn = @control_active_set;
            max_outer_iter = 1;    % fmincon 一次求解占3个时间slot
        case 'interior-point'
            controller_fn = @control_interior_point;
            max_outer_iter = 1;    % fmincon 一次求解占3个时间slot
        otherwise
            error('run_one_case:UnknownAlgorithm', ...
                'Unknown algorithm: %s', config.algorithm);
    end

    %% =====================================================
    % 从 config / scenario 提取仿真参数
    % ======================================================
    num_steps = config.num_steps;
    num_wheels = size(config.wheel_pos, 1);

    % 生成参考轨迹
    path = generateReference(config, config.num_path_pts);

    % 初始状态: 优先从 scenario 取, 否则用默认值
    if isfield(scenario, 'initialState') && ~isempty(scenario.initialState)
        state = scenario.initialState(:);   % [x; y; psi]
    else
        state = [0.05; 0.1; 0.2];
    end

    % 初始车体速度: 优先从 scenario 取, 否则用默认值
    if isfield(scenario, 'initialVelocity') && ~isempty(scenario.initialVelocity)
        lastBodyVelocity = scenario.initialVelocity(:);  % [vx; vy; omega]
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

    % 求解时间: 区分三种时间度量 [P2修复]
    solverTimes_per_subproblem = zeros(1, max_outer_iter * num_steps);
    solverTimes_per_mpc_step = zeros(1, num_steps);
    solverTimes_per_iteration = zeros(1, num_steps);

    solverBackendPerStep = cell(1, num_steps);

    success = true;
    failureReason = '';
    solved_count = 0;

    % 初始化全局求解时间变量
    global solver_time_array;
    solver_time_array = [];

    % 记录初始状态
    states(:, 1) = state;

    for k = 1:num_steps
        try
            [u, worldVelocity, bodyVelocity] = controller_fn( ...
                path, k, lastBodyVelocity, state', config ...
            );

            [wheelSpeed, wheelAngle] = computeWheelOutputs(bodyVelocity, config);

            % 记录执行的控制增量 [P1新增]
            executedU(:, k) = u(:, 1);

            % 记录控制前的状态（索引 k，初始状态在索引 1）
            % [P2修复] 状态与时间对齐: state(:,k) 是第k步控制前的状态
            states(:, k) = state;
            worldVelocities(:, k) = worldVelocity;
            bodyVelocities(:, k) = bodyVelocity;
            wheelSpeeds(:, k) = wheelSpeed;
            wheelAngles(:, k) = wheelAngle;
            solved_count = solved_count + 1;

            % 推进状态
            state = propagateState(state, worldVelocity, config);
            lastBodyVelocity = bodyVelocity;

            % 记录控制后的状态 [P2修复: 保存最终状态]
            states(:, k+1) = state;

        catch ME
            success = false;
            failureReason = getReport(ME);
            break;
        end
    end

    % [P2修复] 正确聚合求解时间
    % solver_time_array 中有 max_outer_iter * solved_count 个元素
    % time_index = max_outer_iter * (step - 1) + m
    if ~isempty(solver_time_array)
        n_total_times = length(solver_time_array);

        % 保存所有 subproblem 时间
        n_save = min(n_total_times, max_outer_iter * solved_count);
        solverTimes_per_subproblem(1:n_save) = solver_time_array(1:n_save);

        % 按 MPC 步聚合: 每步有 max_outer_iter=3 次迭代
        for k_step = 1:solved_count
            idx_start = max_outer_iter * (k_step - 1) + 1;
            idx_end = max_outer_iter * k_step;
            if idx_end <= n_save
                % 每步总时间 = 3次迭代之和
                solverTimes_per_mpc_step(k_step) = sum(solverTimes_per_subproblem(idx_start:idx_end));
                % 每步第一次迭代时间 (对应论文 Fig.4 的 "单次迭代时间")
                solverTimes_per_iteration(k_step) = solverTimes_per_subproblem(idx_start);
            end
        end
    end

    % 截断到实际完成的步数
    if solved_count == 0
        % 控制器第一步就失败: 保留初始状态
        states = states(:, 1);
        worldVelocities = zeros(3, 0);
        bodyVelocities = zeros(3, 0);
        executedU = zeros(3, 0);
        wheelSpeeds = zeros(num_wheels, 0);
        wheelAngles = zeros(num_wheels, 0);
        solverTimes_per_mpc_step = [];
        solverTimes_per_iteration = [];
    else
        states = states(:, 1:solved_count+1);  % 包含最终状态
        worldVelocities = worldVelocities(:, 1:solved_count);
        bodyVelocities = bodyVelocities(:, 1:solved_count);
        executedU = executedU(:, 1:solved_count);
        wheelSpeeds = wheelSpeeds(:, 1:solved_count);
        wheelAngles = wheelAngles(:, 1:solved_count);
        solverTimes_per_mpc_step = solverTimes_per_mpc_step(1:solved_count);
        solverTimes_per_iteration = solverTimes_per_iteration(1:solved_count);
    end

    % 计算指标 (使用 per-MPC-step 时间作为主度量)
    metrics = computeMetrics(path, states, config, ...
        solverTimes_per_mpc_step, wheelSpeeds, wheelAngles, executedU);
    metrics.successRate = solved_count / num_steps;

    % [P2新增] 保存多种时间度量
    metrics.meanSolveTimePerStep = mean(solverTimes_per_mpc_step);
    metrics.maxSolveTimePerStep = max(solverTimes_per_mpc_step);
    metrics.totalSolveTimePerStep = sum(solverTimes_per_mpc_step);
    metrics.meanSolveTimePerIteration = mean(solverTimes_per_iteration);
    metrics.maxSolveTimePerIteration = max(solverTimes_per_iteration);
    metrics.solverTimesPerSubproblem = solverTimes_per_subproblem(1:min(length(solverTimes_per_subproblem), max_outer_iter*solved_count));

    % 组装输出
    summary = struct();
    summary.success = success;
    summary.failureReason = failureReason;
    summary.scenario = scenario;
    summary.config = config;
    summary.algorithm = algorithm;              % 记录算法标签
    summary.states = states;
    summary.worldVelocities = worldVelocities;
    summary.bodyVelocities = bodyVelocities;
    summary.executedU = executedU;           % [P1新增]
    summary.wheelSpeeds = wheelSpeeds;
    summary.wheelAngles = wheelAngles;
    summary.solverTimes = solverTimes_per_mpc_step;
    summary.solverTimesPerIteration = solverTimes_per_iteration;
    summary.path = path;
    summary.metrics = metrics;
    summary.numWheels = num_wheels;          % [P2新增] 记录轮数

    if success
        fprintf('[Case %d: %s] RMSE=%.6f, meanSolveTime=%.4fs, J_total=%.4f, successRate=%.0f%%\n', ...
            scenario.id, scenario.name, metrics.rmse, metrics.meanSolveTime, ...
            metrics.trajectoryCost, metrics.successRate*100);
    else
        fprintf('[Case %d: %s] FAILED at step %d/%d: %s\n', ...
            scenario.id, scenario.name, solved_count, num_steps, failureReason);
    end
end
