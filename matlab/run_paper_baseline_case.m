function summary = run_paper_baseline_case(config, scenario)
% RUN_PAPER_BASELINE_CASE 基线专用闭环仿真 runner
%
% 与 run_one_case 的区别:
%   1. 捕获 solver_info (exitflag, iterations, finite)
%   2. 根据 exitflag < 0 或非有限解标记步骤失败
%   3. 记录每步 exitflag 和约束残差
%   4. 计算 warm-up 排除后的中位数耗时 (P1-4)
%   5. 仅支持基线算法: e-lmpc, active-set, interior-point
%
% 输入/输出格式与 run_one_case 一致, 额外增加:
%   summary.solverInfo: 每步的 exitflag/iterations/finite/warmstarted
%   summary.stepDiagnostics: 每步约束残差
%   summary.timing.medianSolveTime: warm-up 后中位数
%   summary.timing.warmupExcluded: 排除的前 N 步
%
% 用法:
%   summary = run_paper_baseline_case(cfg, scen)
%
% 独立脚本: 删除本文件不影响 run_one_case 和 compare_algorithms。

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

    algorithm = lower(config.algorithm);

    switch algorithm
        case 'e-lmpc'
            controller_fn = @control_eLMPC;
        case 'active-set'
            controller_fn = @control_active_set;
        case 'interior-point'
            controller_fn = @control_interior_point;
        otherwise
            error('run_paper_baseline_case:UnsupportedAlgorithm', ...
                '仅支持基线算法, 收到: %s', algorithm);
    end

    %% =====================================================
    % 参数提取
    % ======================================================
    num_steps = config.num_steps;
    num_wheels = size(config.wheel_pos, 1);

    path = generateReference(config, config.num_path_pts);

    if isfield(scenario, 'initialState') && ~isempty(scenario.initialState)
        state = scenario.initialState(:);
    else
        state = [0.05; 0.1; 0.2];
    end

    if isfield(scenario, 'initialVelocity') && ~isempty(scenario.initialVelocity)
        lastBodyVelocity = scenario.initialVelocity(:);
    else
        lastBodyVelocity = [0.01; 0.01; 0.01];
    end

    %% =====================================================
    % 预分配
    % ======================================================
    states = zeros(3, num_steps + 1);
    worldVelocities = zeros(3, num_steps);
    bodyVelocities = zeros(3, num_steps);
    wheelSpeeds = zeros(num_wheels, num_steps);
    wheelAngles = zeros(num_wheels, num_steps);
    executedU = zeros(3, num_steps);

    solveTimes = zeros(1, num_steps);

    % 每步求解器诊断
    stepExitflags = zeros(1, num_steps);
    stepIterations = zeros(1, num_steps);
    stepFinite = true(1, num_steps);
    stepWarmstarted = false(1, num_steps);
    stepFailed = false(1, num_steps);

    % 全局求解时间变量 (由控制器写入)
    global solver_time_array;
    solver_time_array = [];

    success = true;
    failureReason = '';
    solved_count = 0;
    first_fail_step = 0;

    states(:, 1) = state;

    %% =====================================================
    % 闭环仿真
    % ======================================================
    for k = 1:num_steps
        step_tic = tic;

        try
            % 基线控制器返回 4 个输出: u, worldVel, bodyVel, solver_info
            [u, worldVelocity, bodyVelocity, solver_info] = controller_fn( ...
                path, k, lastBodyVelocity, state', config ...
            );

            step_time = toc(step_tic);
            solveTimes(k) = step_time;

            % 记录 solver_info
            stepExitflags(k) = solver_info.exitflag;
            stepIterations(k) = solver_info.iterations;
            stepFinite(k) = solver_info.finite;
            if isfield(solver_info, 'warmstarted')
                stepWarmstarted(k) = solver_info.warmstarted;
            end

            % [P1-3] 检查求解是否真正成功
            step_ok = true;
            fail_reason = '';

            if solver_info.exitflag < 0
                step_ok = false;
                fail_reason = sprintf('exitflag=%d', solver_info.exitflag);
            end

            if ~solver_info.finite
                step_ok = false;
                if isempty(fail_reason)
                    fail_reason = 'solution contains NaN/Inf';
                end
            end

            if any(isnan(u(:))) || any(isinf(u(:)))
                step_ok = false;
                fail_reason = 'u contains NaN/Inf';
            end

            if ~step_ok
                stepFailed(k) = true;
                fprintf('[Step %d: %s] 求解器诊断: %s (exitflag=%d, finite=%d)\n', ...
                    k, algorithm, fail_reason, solver_info.exitflag, solver_info.finite);
                % 仍继续执行 (用回退 u 或解出的 u), 但记录为失败步
            end

            % 执行控制
            [wheelSpeed, wheelAngle] = computeWheelOutputs(bodyVelocity, config);

            executedU(:, k) = u(:, 1);
            states(:, k) = state;
            worldVelocities(:, k) = worldVelocity;
            bodyVelocities(:, k) = bodyVelocity;
            wheelSpeeds(:, k) = wheelSpeed;
            wheelAngles(:, k) = wheelAngle;
            solved_count = solved_count + 1;

            % 推进状态
            state = propagateState(state, worldVelocity, config);
            lastBodyVelocity = bodyVelocity;
            states(:, k+1) = state;

        catch ME
            step_time = toc(step_tic);
            solveTimes(k) = step_time;
            stepFailed(k) = true;
            stepExitflags(k) = -999;
            stepFinite(k) = false;

            if first_fail_step == 0
                first_fail_step = k;
            end
            success = false;
            failureReason = getReport(ME);
            fprintf('[%s: step %d] EXCEPTION: %s\n', algorithm, k, ME.message);
            for si = 1:length(ME.stack)
                fprintf('  at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
            end
            break;
        end
    end

    %% =====================================================
    % 截断到实际完成步数
    % ======================================================
    if solved_count == 0
        states = states(:, 1);
        worldVelocities = zeros(3, 0);
        bodyVelocities = zeros(3, 0);
        executedU = zeros(3, 0);
        wheelSpeeds = zeros(num_wheels, 0);
        wheelAngles = zeros(num_wheels, 0);
        solveTimes = [];
    else
        states = states(:, 1:solved_count+1);
        worldVelocities = worldVelocities(:, 1:solved_count);
        bodyVelocities = bodyVelocities(:, 1:solved_count);
        executedU = executedU(:, 1:solved_count);
        wheelSpeeds = wheelSpeeds(:, 1:solved_count);
        wheelAngles = wheelAngles(:, 1:solved_count);
        solveTimes = solveTimes(1:solved_count);
        stepExitflags = stepExitflags(1:solved_count);
        stepIterations = stepIterations(1:solved_count);
        stepFinite = stepFinite(1:solved_count);
        stepWarmstarted = stepWarmstarted(1:solved_count);
        stepFailed = stepFailed(1:solved_count);
    end

    %% =====================================================
    % 指标计算
    % ======================================================
    metrics = computeMetrics(path, states, config, ...
        solveTimes, wheelSpeeds, wheelAngles, executedU);
    metrics.successRate = solved_count / num_steps;

    % [P1-3] 步骤级成功率: exitflag >= 0 且有限解
    n_valid_steps = sum(~stepFailed);
    metrics.validStepRate = n_valid_steps / max(solved_count, 1);

    % [P1-4] 耗时统计: 含 warm-up 排除
    n_warmup = min(5, max(1, floor(solved_count * 0.05)));  % 排除前 5 步或 5%
    if solved_count > n_warmup
        post_warmup_times = solveTimes(n_warmup+1:end);
        metrics.meanSolveTime = mean(solveTimes);
        metrics.maxSolveTime = max(solveTimes);
        metrics.totalSolveTime = sum(solveTimes);
        metrics.medianSolveTime = median(post_warmup_times);
        metrics.q1SolveTime = prctile(post_warmup_times, 25);
        metrics.q3SolveTime = prctile(post_warmup_times, 75);
        metrics.meanSolveTimePostWarmup = mean(post_warmup_times);
    else
        metrics.meanSolveTime = mean(solveTimes);
        metrics.maxSolveTime = max(solveTimes);
        metrics.totalSolveTime = sum(solveTimes);
        metrics.medianSolveTime = median(solveTimes);
        metrics.q1SolveTime = prctile(solveTimes, 25);
        metrics.q3SolveTime = prctile(solveTimes, 75);
        metrics.meanSolveTimePostWarmup = metrics.meanSolveTime;
    end
    metrics.warmupExcluded = n_warmup;

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
    summary.path = path;
    summary.metrics = metrics;
    summary.numWheels = num_wheels;

    % [P1-3] 每步求解器诊断
    summary.solverInfo = struct();
    summary.solverInfo.exitflags = stepExitflags;
    summary.solverInfo.iterations = stepIterations;
    summary.solverInfo.finite = stepFinite;
    summary.solverInfo.warmstarted = stepWarmstarted;
    summary.solverInfo.stepFailed = stepFailed;
    summary.solverInfo.nValidSteps = n_valid_steps;
    summary.solverInfo.nFailedSteps = sum(stepFailed);
    summary.solverInfo.firstFailStep = first_fail_step;

    % [P1-4] 耗时细节
    summary.timing = struct();
    summary.timing.solveTimes = solveTimes;
    summary.timing.medianSolveTime = metrics.medianSolveTime;
    summary.timing.q1SolveTime = metrics.q1SolveTime;
    summary.timing.q3SolveTime = metrics.q3SolveTime;
    summary.timing.warmupExcluded = n_warmup;

    % 打印结果
    if success
        fprintf('[%s: %s] RMSE=%.6f, medianSolveTime=%.4fs, J_total=%.4f, validSteps=%d/%d\n', ...
            algorithm, scenario.name, metrics.rmse, metrics.medianSolveTime, ...
            metrics.trajectoryCost, n_valid_steps, solved_count);
    else
        fprintf('[%s: %s] FAILED at step %d/%d: %s\n', ...
            algorithm, scenario.name, first_fail_step, num_steps, failureReason);
    end

    % 打印 exitflag 分布
    if solved_count > 0
        unique_flags = unique(stepExitflags);
        fprintf('  exitflag 分布: ');
        for fi = 1:length(unique_flags)
            n_flag = sum(stepExitflags == unique_flags(fi));
            fprintf('%d→%d步  ', unique_flags(fi), n_flag);
        end
        fprintf('\n');
        if sum(stepFailed) > 0
            fprintf('  失败步 (exitflag<0 或非有限): %d/%d\n', sum(stepFailed), solved_count);
        end
    end
end
