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
        {'proposed-3iter', 'e-lmpc', 'interior-point', 'active-set'}, ...
        {fullfile(algorithms_dir, 'RSS_proposed'), ...
         fullfile(algorithms_dir, 'RSS_sqp'), ...
         fullfile(algorithms_dir, 'RSS_fmincon'), ...
         fullfile(algorithms_dir, 'RSS_active_set')} ...
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
            addpath(submodule_dirs('active-set'));
            alg_params = feval('config');
            rmpath(submodule_dirs('active-set'));
        otherwise
            error('run_paper_baseline_case:UnsupportedAlgorithm', ...
                '不支持算法: %s', algorithm);
    end
    alg_params.algorithm = algorithm;
    config = alg_params;

    %% =====================================================
    % 循环外一次性 addpath 当前算法的 submodule
    % 避免循环内 addpath/rmpath 切换导致 control_RSS 函数缓存混乱
    % (此前 clear functions 不足以解决多版本同名函数解析冲突)
    % ======================================================
    switch algorithm
        case {'proposed-3iter', 'e-lmpc', 'interior-point', 'active-set'}
            addpath(submodule_dirs(algorithm));
            % onCleanup: 函数退出时自动 rmpath, 避免路径残留影响后续算法
            sub_to_remove = submodule_dirs(algorithm);
            cleanup_obj = onCleanup(@() rmpath(sub_to_remove));
            clear functions  % 清除残留的 control_RSS 解析, 确保用到当前 submodule 版本
            % proposed-3iter: HPIPM 求解器, 强制 reload Python 模块确保使用最新代码
            if strcmp(algorithm, 'proposed-3iter')
                try
                    % hpipm_qp_solver.py 现位于 algorithms/RSS_proposed/ (与 control_RSS.m 同目录)
                    rss_proposed_dir = submodule_dirs('proposed-3iter');
                    if exist(rss_proposed_dir, 'dir')
                        sys_mod = py.importlib.import_module('sys');
                        py.getattr(sys_mod, 'path').append(rss_proposed_dir);
                    end
                    % 强制重新加载 Python 模块 (清除 MATLAB Engine 的 Python 模块缓存)
                    solver_mod = py.importlib.import_module('hpipm_qp_solver');
                    py.importlib.reload(solver_mod);
                    fprintf('[proposed-3iter] Python HPIPM 求解器模块已重新加载\n');
                catch err
                    fprintf('[proposed-3iter] 警告: Python 模块 reload 失败: %s\n', err.message);
                end
            end
    end

    %% =====================================================
    % 参数提取
    % ======================================================
    num_steps = config.num_steps;
    num_wheels = size(config.wheel_pos, 1);

    % 路径生成 (所有算法参数已统一为 m 单位, 路径点数用各自 config.num_path_pts)
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
            % 默认 solver_info (非 fmincon 算法)
            step_solver_info = struct('exitflag', NaN, 'success', true, ...
                'max_ineq_violation', 0, 'max_eq_violation', 0);

            % 按算法名分发 (与 run_one_case 一致, 但不使用外层 tic/toc 计时,
            % 而是使用 submodule 返回的 solve_time, 更准确)
            switch algorithm
                case 'proposed-3iter'
                    % RSS_proposed 0121: [u, new_state_dot, velocity] = control_RSS(path, step, state_dot, state)
                    % 0121 版 cvx_solver ECOS 已在 cvx_begin 之后, 无需外部 CVX status 重置
                    clear K H R xInit solver_time_array
                    % 首步诊断: 确认 MATLAB 实际加载的 control_RSS.m 版本
                    if k == 1
                        crss_path = which('control_RSS');
                        fprintf('[proposed-3iter 诊断] control_RSS.m 路径: %s\n', crss_path);
                        try
                            fid = fopen(crss_path, 'r');
                            if fid ~= -1
                                line34 = '';
                                for li = 1:34
                                    line34 = fgetl(fid);
                                end
                                fclose(fid);
                                fprintf('[proposed-3iter 诊断] control_RSS.m 第34行: %s\n', line34);
                            end
                        catch
                        end
                        try
                            [sel_solver, avail_solvers] = cvx_solver;
                            fprintf('[proposed-3iter 诊断] CVX 当前 solver: %s\n', sel_solver);
                            fprintf('[proposed-3iter 诊断] CVX 可用 solvers: %s\n', strjoin(avail_solvers, ', '));
                        catch
                        end
                    end
                    [u_full, worldVelocity, bodyVelocity] = ...
                        control_RSS(path, k, lastBodyVelocity, state');
                    u = u_full(:, 1);
                    solve_time = NaN;  % submodule 不输出 solve_time (用 diary 捕获, 此处不解析)
                    iter_num = NaN;    % submodule 不输出 iter_num

                case 'e-lmpc'
                    % RSS_sqp: [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state)
                    % submodule 路径已在循环外一次性 addpath
                    [worldVelocity, bodyVelocity, solve_time, iter_num] = ...
                        control_RSS(path, k, lastBodyVelocity, state');
                    u = bodyVelocity - lastBodyVelocity;

                case 'interior-point'
                    % RSS_fmincon: [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, step, state_dot, state, params)
                    % submodule 路径已在循环外一次性 addpath
                    [worldVelocity, bodyVelocity, solve_time, iter_num] = ...
                        control_RSS(path, k, lastBodyVelocity, state', config);
                    u = bodyVelocity - lastBodyVelocity;

                case 'active-set'
                    % RSS_active_set: [new_state_dot, velocity, solve_time, iter_num, solver_info] = control_RSS(...)
                    % 第5输出 solver_info 含 exitflag/约束违反量, 用于判断求解是否真正成功
                    [worldVelocity, bodyVelocity, solve_time, iter_num, step_solver_info] = ...
                        control_RSS(path, k, lastBodyVelocity, state', config);
                    u = bodyVelocity - lastBodyVelocity;

                otherwise
                    error('run_paper_baseline_case:UnsupportedAlgorithm', ...
                        '不支持算法: %s', algorithm);
            end

            solveTimes(k) = solve_time;
            stepIterations(k) = iter_num;
            stepExitflags(k) = step_solver_info.exitflag;  % active-set 记录 exitflag

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

            % active-set: 基于 exitflag 和约束违反量判断成功
            % (exitflag<=0 或约束违反量>1e-6 视为求解失败, 不应记为有效步)
            if strcmp(algorithm, 'active-set') && ~step_solver_info.success
                step_ok = false;
                if isempty(fail_reason)
                    fail_reason = sprintf('fmincon failed: exitflag=%d, ineq_viol=%.2e, eq_viol=%.2e', ...
                        step_solver_info.exitflag, step_solver_info.max_ineq_violation, ...
                        step_solver_info.max_eq_violation);
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
            % 所有算法的 control_RSS 期望传入 body frame 速度 (state_dot):
            %   - RSS_proposed 0121: last_vel = velocity (第3输出, body frame)
            %   - RSS_sqp/RSS_fmincon: state_dot = bodyVelocity
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
