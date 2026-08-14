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
    % 保存用户传入的 solver 模式 (config 会被 alg_params 覆盖)
    user_solver = '';
    if isfield(config, 'solver') && ~isempty(config.solver)
        user_solver = config.solver;
    end

    alg_params.algorithm = algorithm;
    config = alg_params;

    % 恢复 solver 模式 ('ocpqp' 默认 / 'denseqcqp' 精确二次约束)
    if ~isempty(user_solver)
        config.solver = user_solver;
    end

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
                    % 强制重新加载 Python 模块 (用 sys.modules.pop + import, 比 reload 可靠)
                    % reload 不会重新解析文件路径, 必须先 pop 再 import
                    reload_py = fullfile(rss_proposed_dir, '_reload_solver_tmp.py');
                    fid = fopen(reload_py, 'w');
                    fprintf(fid, 'import sys\n');
                    fprintf(fid, 'sys.modules.pop("hpipm_qp_solver", None)\n');
                    fprintf(fid, 'import hpipm_qp_solver\n');
                    fclose(fid);
                    py.runpy.run_path(reload_py);
                    delete(reload_py);
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
    stepApproximate = false(1, num_steps);  % 使用 approximate incumbent (非严格可行)
    step_solver_call_count = zeros(1, num_steps);  % 真实 HPIPM solver.solve() 调用数 (proposed-3iter)
    % 约束违反量 per-step (由 run_paper_baseline_case 维护, 不依赖 persistent)
    step_wheel_viol = zeros(1, num_steps);  % 原始轮速约束违反
    step_cone_viol = zeros(1, num_steps);   % 原始转向锥约束违反
    step_incumbent_type = cell(1, num_steps); % 'strict'/'approximate'/'none'

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
                    % RSS_proposed: [u, new_state_dot, velocity, diagnostics] = control_RSS(path, step, state_dot, state)
                    % 求解器通过 cfg.solver 切换 (三模式 switch-case):
                    %   'ocpqp'    (默认) → control_RSS           (HPIPM OCP QP + 线性化约束, legacy)
                    %   'ocpqcqp'         → control_RSS_ocpqcqp   (HPIPM OCP QCQP, 精确二次约束)
                    %   'denseqcqp'       → control_RSS_denseqcqp (HPIPM Dense QCQP, Golden oracle)
                    global RSS_WARMSTART_UHAT RSS_SOLVER_MODE;
                    if k == 1
                        RSS_WARMSTART_UHAT = [];
                        if ~isfield(config, 'solver') || isempty(config.solver)
                            config.solver = 'ocpqcqp';  % 默认
                        end
                        solver_mode = lower(config.solver);
                        switch solver_mode
                            case 'ocpqp'
                                fprintf('[proposed-3iter] 求解器: control_RSS (OCP QP, 线性化约束, legacy)\n');
                            case 'ocpqcqp'
                                fprintf('[proposed-3iter] 求解器: control_RSS_ocpqcqp (OCP QCQP, 精确二次约束)\n');
                            case 'denseqcqp'
                                fprintf('[proposed-3iter] 求解器: control_RSS_denseqcqp (Dense QCQP, Golden oracle)\n');
                            otherwise
                                warning('未知 solver mode: %s, 使用默认 ocpqp', config.solver);
                                solver_mode = 'ocpqp';
                        end
                        RSS_SOLVER_MODE = solver_mode;
                    end
                    solver_mode = RSS_SOLVER_MODE;
                    clear K H R xInit control_RSS control_RSS_denseqcqp control_RSS_ocpqcqp
                    switch solver_mode
                        case 'ocpqp'
                            [u_full, worldVelocity, bodyVelocity, diagnostics] = ...
                                control_RSS(path, k, lastBodyVelocity, state');
                        case 'ocpqcqp'
                            [u_full, worldVelocity, bodyVelocity, diagnostics] = ...
                                control_RSS_ocpqcqp(path, k, lastBodyVelocity, state');
                        case 'denseqcqp'
                            [u_full, worldVelocity, bodyVelocity, diagnostics] = ...
                                control_RSS_denseqcqp(path, k, lastBodyVelocity, state');
                    end
                    u = u_full(:, 1);
                    solve_time = diagnostics.total_solve_time;  % SCP 迭代总耗时
                    % 使用真实 solver_call_count (不是 max_iter, 后者始终为 3 但不代表实际调用次数)
                    iter_num = diagnostics.solver_call_count;
                    step_solver_call_count(k) = diagnostics.solver_call_count;

                    % ========== 记录 per-step 约束违反量 (不用 persistent) ==========
                    if isfield(diagnostics, 'orig_wheel_viol_final')
                        step_wheel_viol(k) = diagnostics.orig_wheel_viol_final;
                        step_cone_viol(k) = diagnostics.orig_cone_viol_final;
                    end
                    if isfield(diagnostics, 'incumbent_type')
                        step_incumbent_type{k} = diagnostics.incumbent_type;
                    end
                    if isfield(diagnostics, 'step_approximate') && diagnostics.step_approximate
                        stepApproximate(k) = true;
                    end

                    % ========== 检查 step_failed: 不得忽略 ==========
                    % control_RSS.m 在三次 outer 后无任何 incumbent 时设置 step_failed=true
                    % 此时 u = u_seed (未验证), 不得用于推进闭环状态
                    if isfield(diagnostics, 'step_failed') && diagnostics.step_failed
                        step_ok = false;
                        fail_reason = 'diagnostics.step_failed=true (三次 outer 后无 incumbent)';
                        stepFailed(k) = true;
                        stepFinite(k) = false;
                        if first_fail_step == 0
                            first_fail_step = k;
                        end
                        success = false;
                        failureReason = fail_reason;
                        fprintf('[%s: step %d] STEP FAILED: %s\n', algorithm, k, fail_reason);
                        break;
                    end

                    % ========== 检查 solver_call_count == 3 ==========
                    if diagnostics.solver_call_count ~= 3
                        fprintf('[%s: step %d] 警告: solver_call_count=%d (预期 3)\n', ...
                            algorithm, k, diagnostics.solver_call_count);
                    end

                    % ========== Warm start: 由 control_RSS.m 内部管理 RSS_WARMSTART_UHAT ==========
                    % control_RSS.m 已保存 RSS_WARMSTART_UHAT = u (3×K 矩阵), 此处不再覆盖

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
        stepExitflags = stepExitflags(1:0);
        stepIterations = stepIterations(1:0);
        stepFinite = stepFinite(1:0);
        stepWarmstarted = stepWarmstarted(1:0);
        stepFailed = stepFailed(1:0);
        stepApproximate = stepApproximate(1:0);
        step_solver_call_count = step_solver_call_count(1:0);
        step_wheel_viol = step_wheel_viol(1:0);
        step_cone_viol = step_cone_viol(1:0);
        step_incumbent_type = step_incumbent_type(1:0);
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
        stepApproximate = stepApproximate(1:solved_count);
        step_solver_call_count = step_solver_call_count(1:solved_count);
        step_wheel_viol = step_wheel_viol(1:solved_count);
        step_cone_viol = step_cone_viol(1:solved_count);
        step_incumbent_type = step_incumbent_type(1:solved_count);
    end

    %% =====================================================
    % 指标计算
    % ======================================================
    metrics = computeMetrics(path, states, config, ...
        solveTimes, wheelSpeeds, wheelAngles, executedU);
    metrics.successRate = solved_count / num_steps;

    % [P1-3] 步骤级成功率: 解有限且未标记失败
    n_valid_steps = sum(~stepFailed);
    if solved_count == 0
        n_valid_steps = 0;
    end
    metrics.validStepRate = n_valid_steps / max(solved_count, 1);

    % J 分解: J_position, J_heading, J_control, J_total
    % J_total = sum_k [ 30*(ex^2+ey^2) + 1*e_psi^2 + 0.3*||u_k||^2 ]
    Q_w = [30, 30, 1];      % [w_pos, w_pos, w_psi]
    R_w = [0.3, 0.3, 0.3];
    J_position = 0; J_heading = 0; J_control = 0;
    num_ref = size(path, 2);
    for kk = 1:solved_count
        ref_idx = min(kk, num_ref);
        e = states(:, kk) - path(:, ref_idx);
        e(3) = mod(e(3) + pi, 2*pi) - pi;  % wrap angle
        J_position = J_position + Q_w(1) * (e(1)^2 + e(2)^2);
        J_heading  = J_heading  + Q_w(3) * e(3)^2;
        J_control  = J_control  + R_w(1) * sum(executedU(:, kk).^2);
    end
    metrics.J_position = J_position;
    metrics.J_heading = J_heading;
    metrics.J_control = J_control;
    metrics.J_total = J_position + J_heading + J_control;

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
    % 约束违反量汇总 (由 run_paper_baseline_case 维护, 不依赖 persistent)
    % ======================================================
    if solved_count > 0
        max_wheel_viol_reported = max(step_wheel_viol);
        max_cone_viol_reported = max(step_cone_viol);
        [max_wheel_viol_reported, step_wheel_max] = max(step_wheel_viol);
        [max_cone_viol_reported, step_cone_max] = max(step_cone_viol);
        cnt_wheel_viol = sum(step_wheel_viol > 1e-6);
        cnt_cone_viol = sum(step_cone_viol > 1e-6);
        cnt_strict = sum(strcmp(step_incumbent_type, 'strict'));
        cnt_approximate = sum(strcmp(step_incumbent_type, 'approximate'));
        cnt_none = sum(strcmp(step_incumbent_type, 'none'));
    else
        max_wheel_viol_reported = 0; max_cone_viol_reported = 0;
        step_wheel_max = 0; step_cone_max = 0;
        cnt_wheel_viol = 0; cnt_cone_viol = 0;
        cnt_strict = 0; cnt_approximate = 0; cnt_none = 0;
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
    summary.solverInfo.stepApproximate = stepApproximate;
    summary.solverInfo.nValidSteps = n_valid_steps;
    summary.solverInfo.nFailedSteps = sum(stepFailed);
    summary.solverInfo.nApproximateSteps = sum(stepApproximate);
    summary.solverInfo.firstFailStep = first_fail_step;
    summary.solverInfo.solverCallCounts = step_solver_call_count;  % 真实 HPIPM solver.solve() 调用数

    % 约束违反量汇总 (per-step 数组 + 汇总)
    summary.constraintInfo = struct();
    summary.constraintInfo.wheelViolPerStep = step_wheel_viol;
    summary.constraintInfo.coneViolPerStep = step_cone_viol;
    summary.constraintInfo.incumbentTypePerStep = step_incumbent_type;
    summary.constraintInfo.maxWheelViol = max_wheel_viol_reported;
    summary.constraintInfo.maxConeViol = max_cone_viol_reported;
    summary.constraintInfo.stepWheelMax = step_wheel_max;
    summary.constraintInfo.stepConeMax = step_cone_max;
    summary.constraintInfo.cntWheelViol = cnt_wheel_viol;
    summary.constraintInfo.cntConeViol = cnt_cone_viol;
    summary.constraintInfo.cntStrictIncumbent = cnt_strict;
    summary.constraintInfo.cntApproximateIncumbent = cnt_approximate;
    summary.constraintInfo.cntNoneIncumbent = cnt_none;

    % [P1-4] 耗时细节
    summary.timing = struct();
    summary.timing.solveTimes = solveTimes;
    summary.timing.medianSolveTime = metrics.medianSolveTime;
    summary.timing.q1SolveTime = metrics.q1SolveTime;
    summary.timing.q3SolveTime = metrics.q3SolveTime;
    summary.timing.warmupExcluded = n_warmup;

    % 打印结果 (仓库 README 风格: 简洁汇总行)
    if success
        fprintf('[%s: %s] RMSE=%.6f, medianSolveTime=%.4fs, J_total=%.4f, validSteps=%d/%d\n', ...
            algorithm, scenario.name, metrics.rmse, metrics.medianSolveTime, ...
            metrics.J_total, n_valid_steps, solved_count);
    else
        fprintf('[%s: %s] FAILED at step %d/%d: %s\n', ...
            algorithm, scenario.name, first_fail_step, num_steps, failureReason);
    end

    % 打印 iterations 分布
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

        % 约束违反量汇总 (proposed-3iter 专用, 简洁格式)
        if strcmp(algorithm, 'proposed-3iter')
            fprintf('  solver_calls/step: %d | incumbent: strict=%d, approx=%d, none=%d\n', ...
                step_solver_call_count(1), cnt_strict, cnt_approximate, cnt_none);
            fprintf('  约束违反: maxWheel=%.2e (step %d), maxCone=%.2e (step %d)\n', ...
                max_wheel_viol_reported, step_wheel_max, max_cone_viol_reported, step_cone_max);
        end
    end
end
