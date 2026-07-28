function summary = run_paper_baseline_case(config, scenario)
% RUN_PAPER_BASELINE_CASE 闭环仿真 runner (论文 Section IV 复现)
%
% 与 run_one_case 的区别:
%   1. 计算 warm-up 排除后的中位数耗时 (P1-4 复现要求)
%   2. 记录每步 iter_num (submodule 返回的求解器迭代数)
%   3. 检查解的有限性 (NaN/Inf) 标记失败步
%   4. 支持全部 4 种算法: proposed-3iter, e-lmpc, active-set, interior-point
%
% 算法调用架构 (与 run_one_case 一致):
%   - proposed-3iter  → algorithms/RSS_proposed/control_RSS.m  (git submodule)
%   - e-lmpc          → algorithms/RSS_sqp/control_RSS.m       (git submodule)
%   - interior-point  → algorithms/RSS_fmincon/control_RSS.m   (git submodule)
%   - active-set      → algorithms/control_active_set.m        (本地实现)
%
% 注: 每个算法使用各自 submodule 的 config.m 参数 (不使用外部传入的 config 覆盖)。
%     submodule 接口不返回 exitflag/warmstarted, 相关字段记为 NaN/false。
%     iter_num 从 submodule 的第 4 个输出获取 (e-lmpc/interior-point)。
%     RSS_proposed 只输出 new_state_dot, 缺失的 u/solve_time/iter_num 记为 NaN。
%
% 输出额外字段:
%   summary.solverInfo: 每步的 iterations/finite/stepFailed
%   summary.timing: warm-up 排除后的中位数/分位数耗时
%
% 用法:
%   summary = run_paper_baseline_case(cfg, scen)

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
        config.algorithm = 'e-lmpc';
    end

    algorithm = lower(config.algorithm);

    % 定位 submodule 路径 (与 run_one_case 一致)
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
    % 加载算法各自的 config.m (每个算法用各自的参数)
    % ======================================================
    switch algorithm
        case 'proposed-3iter'
            addpath(submodule_dirs('proposed-3iter'));
            alg_params = feval('config');
            rmpath(submodule_dirs('proposed-3iter'));
        case 'e-lmpc'
            addpath(submodule_dirs('e-lmpc'));
            alg_params = feval('config');
            rmpath(submodule_dirs('e-lmpc'));
        case 'interior-point'
            addpath(submodule_dirs('interior-point'));
            alg_params = feval('config');
            rmpath(submodule_dirs('interior-point'));
        case 'active-set'
            alg_params = defaultConfig();
        otherwise
            error('run_paper_baseline_case:UnsupportedAlgorithm', ...
                '不支持算法: %s', algorithm);
    end
    alg_params.algorithm = algorithm;
    config = alg_params;

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

    % 每步求解器诊断 (submodule 接口限制, 部分字段记为 NaN)
    stepExitflags = NaN(1, num_steps);       % submodule 不返回 exitflag
    stepIterations = zeros(1, num_steps);    % 从 submodule 第 4 输出获取
    stepFinite = true(1, num_steps);         % 自己检查 NaN/Inf
    stepWarmstarted = false(1, num_steps);   % submodule 不返回
    stepFailed = false(1, num_steps);

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
            % 按算法名分发 (与 run_one_case 一致, 但不使用外层 tic/toc 计时,
            % 而是使用 submodule 返回的 solve_time, 更准确)
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
                    solve_time = NaN;  % submodule 不输出 solve_time
                    iter_num = NaN;    % submodule 不输出 iter_num

                case 'e-lmpc'
                    % RSS_sqp: [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state)
                    addpath(submodule_paths.('e-lmpc'));
                    [worldVelocity, bodyVelocity, solve_time, iter_num] = ...
                        control_RSS(path, k, lastBodyVelocity, state');
                    rmpath(submodule_paths.('e-lmpc'));
                    u = bodyVelocity - lastBodyVelocity;

                case 'interior-point'
                    % RSS_fmincon: [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state, params)
                    addpath(submodule_paths.('interior-point'));
                    [worldVelocity, bodyVelocity, solve_time, iter_num] = ...
                        control_RSS(path, k, lastBodyVelocity, state', config);
                    rmpath(submodule_paths.('interior-point'));
                    u = bodyVelocity - lastBodyVelocity;

                case 'active-set'
                    % 本地: [u, worldVelocity, bodyVelocity] = control_active_set(path, k, state_dot, state, config)
                    [u_full, worldVelocity, bodyVelocity] = ...
                        control_active_set(path, k, lastBodyVelocity, state', config);
                    u = u_full(:, 1);
                    solve_time = toc(step_tic);  % 本地版本不输出 solve_time, 用外层计时
                    iter_num = NaN;

                otherwise
                    error('run_paper_baseline_case:UnsupportedAlgorithm', ...
                        '不支持算法: %s', algorithm);
            end

            solveTimes(k) = solve_time;
            stepIterations(k) = iter_num;

            % [P1-3] 检查解的有限性
            step_ok = true;
            fail_reason = '';

            if any(isnan(worldVelocity(:))) || any(isinf(worldVelocity(:)))
                step_ok = false;
                fail_reason = 'worldVelocity contains NaN/Inf';
            end

            if any(isnan(u(:))) || any(isinf(u(:)))
                step_ok = false;
                if isempty(fail_reason)
                    fail_reason = 'u contains NaN/Inf';
                end
            end

            stepFinite(k) = step_ok;

            if ~step_ok
                stepFailed(k) = true;
                fprintf('[Step %d: %s] 求解器诊断: %s\n', k, algorithm, fail_reason);
                % 仍继续执行, 但记录为失败步
            end

            % 执行控制
            [wheelSpeed, wheelAngle] = computeWheelOutputs(bodyVelocity, config);

            executedU(:, k) = u;
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
            solveTimes(k) = toc(step_tic);
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

    % [P1-3] 步骤级成功率: 解有限且未标记失败
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

    % [P1-3] 每步求解器诊断 (submodule 接口限制, exitflag/warmstarted 为 NaN/false)
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

    % 打印 iterations 分布 (替代原 exitflag 分布)
    if solved_count > 0
        valid_iters = stepIterations(~isnan(stepIterations));
        if ~isempty(valid_iters)
            unique_iters = unique(valid_iters);
            fprintf('  iterations 分布: ');
            for fi = 1:length(unique_iters)
                n_iter = sum(stepIterations == unique_iters(fi));
                fprintf('%d→%d步  ', unique_iters(fi), n_iter);
            end
            fprintf('\n');
        end
        if sum(stepFailed) > 0
            fprintf('  失败步 (解非有限): %d/%d\n', sum(stepFailed), solved_count);
        end
    end
end
