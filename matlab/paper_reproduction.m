function comparison = paper_reproduction(algorithms)
% PAPER_REPRODUCTION  论文 Section IV 固定轨迹实验复现
%
% 使用论文 Section IV 的固定初始状态对四种算法进行仿真:
%   - proposed-3iter
%   - e-lmpc
%   - active-set
%   - interior-point
%
% 每个算法使用各自 submodule 的 config.m 参数 (不在本脚本中覆盖):
%   - proposed-3iter → algorithms/RSS_proposed/config.m
%   - e-lmpc         → algorithms/RSS_sqp/config.m
%   - interior-point → algorithms/RSS_fmincon/config.m
%   - active-set     → defaultConfig.m (本地, 对齐 RSS_fmincon)
%
% 固定初始状态 (论文无随机性):
%   - initialState  = [0.05; 0.1; 0.2]
%   - initialVelocity = [0.01; 0.01; 0.01]
%
% 输出格式与 main 一致:
%   - results/paper_reproduction/per_algorithm/{算法名}/{算法名}_summary.png
%   - results/paper_reproduction/per_algorithm/{算法名}/seed_0001.png
%   - results/paper_reproduction/{算法名}_results.csv
%   - results/paper_reproduction/paper_reproduction.mat
%
% 增量复现: 默认只运行指定算法, 保留未运行算法的已有结果 (从 paper_reproduction.mat 加载)。
% 若需全量重跑, 先删除 results/paper_reproduction/paper_reproduction.mat 或传入全部算法列表。
%
% 用法:
%   paper_reproduction                              % 默认: 重跑全部 4 种算法
%   paper_reproduction({'proposed-3iter'})          % 只重跑 proposed-3iter, 保留其余
%   paper_reproduction({'e-lmpc','active-set'})     % 只重跑指定算法, 保留其余
%
% 完全独立脚本: 删除本文件不影响其他功能。

    %% =====================================================
    % 1. 路径设置
    %% =====================================================
    script_dir = fileparts(mfilename('fullpath'));
    workspace_root = fileparts(script_dir);
    addpath(fullfile(workspace_root, 'matlab'));
    addpath(fullfile(workspace_root, 'algorithms'));

    %% =====================================================
    % 2. 算法列表与场景设置
    % ======================================================
    all_algorithms = {'proposed-3iter', 'e-lmpc', 'active-set', 'interior-point'};
    if nargin < 1 || isempty(algorithms)
        algorithms = all_algorithms;
    end
    % 校验算法名合法性
    for i = 1:length(algorithms)
        if ~any(strcmp(all_algorithms, algorithms{i}))
            error('paper_reproduction:UnknownAlgorithm', ...
                '未知算法: %s。可选: %s', algorithms{i}, strjoin(all_algorithms, ', '));
        end
    end
    num_algorithms = length(algorithms);

    fprintf('========================================\n');
    fprintf('论文 Section IV 固定轨迹实验复现\n');
    fprintf('本次运行算法: %s\n', strjoin(algorithms, ', '));
    fprintf('(每个算法使用各自 submodule 的 config.m 参数)\n');
    fprintf('========================================\n');

    % cfg 仅传递算法名, 实际参数由 run_paper_baseline_case 加载各自 config.m
    cfg = struct();

    % 固定初始状态 (论文无随机性)
    scen = struct();
    scen.name = 'paper_fixed';
    scen.id = 1;
    scen.initialState = [0.05; 0.1; 0.2];
    scen.initialVelocity = [0.01; 0.01; 0.01];

    %% =====================================================
    % 3. 输出目录 (增量模式: 不清空整个目录, 只清理本次运行算法的输出)
    % ======================================================
    results_dir = fullfile(workspace_root, 'results', 'paper_reproduction');
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end

    % 清理本次将重跑算法的旧输出文件 (CSV + 汇总图 + 逐 seed 图目录)
    for i = 1:num_algorithms
        alg = algorithms{i};
        csv_file = fullfile(results_dir, sprintf('%s_results.csv', alg));
        if exist(csv_file, 'file'), delete(csv_file); end
        summary_png = fullfile(results_dir, 'per_algorithm', sprintf('%s_summary.png', alg));
        if exist(summary_png, 'file'), delete(summary_png); end
        per_alg_dir = fullfile(results_dir, 'per_algorithm', alg);
        if exist(per_alg_dir, 'dir')
            try rmdir(per_alg_dir, 's'); catch, end
        end
    end

    %% =====================================================
    % 4. 初始化 comparison 结构体
    % 增量模式: 加载已有 paper_reproduction.mat, 保留未重跑算法的旧结果
    % ======================================================
    final_file = fullfile(results_dir, 'paper_reproduction.mat');
    prev_comparison = [];
    if exist(final_file, 'file')
        try
            loaded = load(final_file, 'comparison');
            prev_comparison = loaded.comparison;
            fprintf('[增量模式] 已加载旧结果: %s\n', final_file);
        catch
            fprintf('[增量模式] 旧 .mat 加载失败, 将全新开始。\n');
        end
    end

    comparison = struct();
    comparison.seeds = 1;
    comparison.numSeeds = 1;
    comparison.totalElapsed = 0;

    if ~isempty(prev_comparison) && isfield(prev_comparison, 'results') && ~isempty(prev_comparison.results)
        % 合并: 保留旧结果中不在本次运行列表里的算法
        kept_results = {};
        kept_pairs = cell(0, 2);
        for i = 1:length(prev_comparison.results)
            r = prev_comparison.results{i};
            if isfield(r, 'algorithm') && ~any(strcmp(algorithms, r.algorithm))
                kept_results{end+1} = r;
                kept_pairs(end+1, :) = {1, r.algorithm};
            end
        end
        comparison.results = kept_results;
        comparison.completedPairs = kept_pairs;
        if ~isempty(kept_results)
            fprintf('[增量模式] 保留已有算法结果: %s\n', ...
                strjoin({kept_pairs{:, 2}}, ', '));
        end
    else
        comparison.results = {};
        comparison.completedPairs = cell(0, 2);
    end

    % comparison.algorithms 始终记录全部 4 种算法 (用于 Table II 汇总打印完整)
    all_algorithms_full = {'proposed-3iter', 'e-lmpc', 'active-set', 'interior-point'};
    comparison.algorithms = all_algorithms_full;
    comparison.numAlgorithms = length(all_algorithms_full);

    %% =====================================================
    % 5. 逐算法运行 (仅本次指定算法)
    % ======================================================
    total_tic = tic;

    for alg_idx = 1:num_algorithms
        algorithm = algorithms{alg_idx};
        cfg.algorithm = algorithm;

        % 切换算法时清除 e-LMPC warm start
        clear functions;

        fprintf('\n--- [%s] 论文固定场景 ---\n', algorithm);

        case_tic = tic;
        try
            summary = run_paper_baseline_case(cfg, scen);
            summary.elapsedTime = toc(case_tic);
        catch ME
            summary = struct();
            summary.success = false;
            summary.failureReason = getReport(ME);
            summary.scenario = scen;
            summary.config = cfg;
            summary.algorithm = algorithm;
            summary.seed = 1;
            summary.elapsedTime = toc(case_tic);
            summary.metrics = struct( ...
                'rmse', NaN, 'trajectoryCost', NaN, 'successRate', 0);
            fprintf('[%s] EXCEPTION: %s\n', algorithm, ME.message);
            for si = 1:length(ME.stack)
                fprintf('  at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
            end
        end

        summary.seed = 1;
        comparison.results{end+1} = summary;
        comparison.completedPairs(end+1, :) = {1, algorithm};

        % 生成图像 + CSV (按算法名在 all_algorithms_full 中的索引定位子图)
        alg_full_idx = find(strcmp(all_algorithms_full, algorithm));
        plot_one_algorithm_paper(comparison, alg_full_idx, results_dir);
        save_algorithm_csv_paper(comparison, alg_full_idx, results_dir);
    end

    comparison.totalElapsed = toc(total_tic);

    %% =====================================================
    % 6. 汇总统计 (覆盖全部 4 种算法, 含保留的旧结果)
    % ======================================================
    fprintf('\n========================================\n');
    for alg_idx = 1:length(all_algorithms_full)
        algorithm = all_algorithms_full{alg_idx};
        alg_results = get_alg_results_paper(comparison, alg_idx);
        if isempty(alg_results)
            fprintf('[算法 %s] 无结果\n', algorithm);
        else
            n_success = sum([alg_results.success]);
            fprintf('[算法 %s] 成功: %d/1, 失败: %d/1\n', ...
                algorithm, n_success, 1 - n_success);
        end
    end

    %% =====================================================
    % 7. 保存完整结果
    % ======================================================
    final_file = fullfile(results_dir, 'paper_reproduction.mat');
    save(final_file, 'comparison', '-v7.3');
    fprintf('\n完整结果已保存到: %s\n', final_file);

    % Table II 格式汇总
    print_table_ii_paper(comparison);
end


%% =========================================================
% 单算法图像生成 (与 compare_algorithms 风格一致)
% ==========================================================

function plot_one_algorithm_paper(comparison, alg_idx, results_dir)

    algorithm = comparison.algorithms{alg_idx};
    alg_results = get_alg_results_paper(comparison, alg_idx);

    n_success = sum([alg_results.success]);
    if n_success == 0
        fprintf('[算法 %s] 无成功案例，跳过绘图\n', algorithm);
        return;
    end

    fig_dir = fullfile(results_dir, 'per_algorithm');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    rmse_vals = arrayfun(@(r) safe_metric(r, 'rmse'), alg_results);
    rmse_vals(~[alg_results.success]) = NaN;

    traj_costs = arrayfun(@(r) safe_metric(r, 'trajectoryCost'), alg_results);
    traj_costs(~[alg_results.success]) = NaN;

    fig = figure('Position', [50, 50, 1200, 500], 'Visible', 'off');

    subplot(1, 3, 1);
    valid_rmse = rmse_vals([alg_results.success]);
    bar(valid_rmse, 'FaceColor', [0.2, 0.6, 0.8]);
    xlabel('Scenario'); ylabel('RMSE (m)');
    title(sprintf('%s: RMSE', algorithm));
    grid on;

    subplot(1, 3, 2);
    valid_costs = traj_costs([alg_results.success]);
    bar(valid_costs, 'FaceColor', [0.8, 0.2, 0.2]);
    xlabel('Scenario'); ylabel('J_total');
    title(sprintf('%s: Trajectory cost', algorithm));
    grid on;

    subplot(1, 3, 3);
    solve_times = arrayfun(@(r) safe_metric(r, 'meanSolveTime'), alg_results) * 1000;
    solve_times(~[alg_results.success]) = NaN;
    valid_times = solve_times([alg_results.success]);
    bar(valid_times, 'FaceColor', [0.2, 0.8, 0.2]);
    xlabel('Scenario'); ylabel('Mean solve time (ms)');
    title(sprintf('%s: Solve time', algorithm));
    grid on;

    fig_file = fullfile(fig_dir, sprintf('%s_summary.png', algorithm));
    saveas(fig, fig_file);
    close(fig);

    fprintf('[算法 %s] 汇总图已保存: %s\n', algorithm, fig_file);

    % ===== 逐 seed 轨迹图 (风格与 main.m 一致) =====
    per_case_dir = fullfile(fig_dir, algorithm);
    if ~exist(per_case_dir, 'dir'), mkdir(per_case_dir); end

    success_mask = [alg_results.success];
    success_results = alg_results(success_mask);

    fn = 'Times New Roman'; fs = 11;

    for i = 1:length(success_results)
        r = success_results(i);
        if ~isfield(r, 'states') || isempty(r.states), continue; end
        if ~isfield(r, 'path') || isempty(r.path), continue; end

        num_w = r.numWheels;
        dt = r.config.dt;
        num_steps = size(r.states, 2) - 1;
        t_vec = (0:num_steps-1) * dt;

        wheel_colors = lines(max(num_w, 5));
        ref_color  = wheel_colors(4, :);
        act_color  = wheel_colors(2, :);
        cons_color = wheel_colors(5, :);

        fig = figure('Position', [100, 100, 1500, 300], 'Visible', 'off');

        % 子图1: 轨迹跟踪性能
        subplot(1, 3, 1);
        hold on;
        plot(r.path(1,:), r.path(2,:), '-', ...
            'Color', ref_color, 'LineWidth', 2, 'DisplayName', 'Reference Trajectory');
        plot(r.states(1,:), r.states(2,:), 'x', ...
            'Color', act_color, 'LineWidth', 1.5, 'DisplayName', 'Simulation Results');
        xlabel('x_w (m)', 'FontName', fn, 'FontSize', fs);
        ylabel('y_w (m)', 'FontName', fn, 'FontSize', fs);
        title('Trajectory Tracking Performance', 'FontName', fn, 'FontSize', fs);
        legend('FontSize', 6, 'FontName', fn);
        grid on; axis equal;
        set(gca, 'FontName', fn, 'FontSize', fs);
        hold off;

        % 子图2: 各轮转向速率
        subplot(1, 3, 2);
        hold on;
        yline(r.config.phidotmax, '--', 'Color', cons_color, ...
            'LineWidth', 2, 'DisplayName', 'Constraints');
        yline(-r.config.phidotmax, '--', 'Color', cons_color, ...
            'LineWidth', 2, 'HandleVisibility', 'off');
        if isfield(r, 'wheelAngles') && size(r.wheelAngles, 2) > 1
            d_ang = diff(r.wheelAngles, 1, 2);
            d_ang = mod(d_ang + pi, 2*pi) - pi;
            phidot = d_ang / dt;
            for w = 1:num_w
                plot(t_vec(2:end), phidot(w,:), '-', ...
                    'Color', wheel_colors(w,:), 'LineWidth', 1.7, ...
                    'DisplayName', sprintf('Wheel %d', w));
            end
        end
        xlabel('Time(s)', 'FontName', fn, 'FontSize', fs);
        ylabel('Steering rate (rad/s)', 'FontName', fn, 'FontSize', fs);
        title('Steering Rate of Each Wheel', 'FontName', fn, 'FontSize', fs);
        legend('FontSize', 6, 'FontName', fn);
        grid on;
        set(gca, 'FontName', fn, 'FontSize', fs);
        % interior-point: 数据范围较小, 强制对齐原论文纵坐标范围 -20~20
        if strcmp(algorithm, 'interior-point')
            ylim([-20, 20]);
        end
        hold off;

        % 子图3: 各轮输出速度
        subplot(1, 3, 3);
        hold on;
        yline(r.config.vimax, '--', 'Color', cons_color, ...
            'LineWidth', 2, 'DisplayName', 'Constraints');
        if isfield(r, 'wheelSpeeds') && ~isempty(r.wheelSpeeds)
            for w = 1:num_w
                plot(t_vec, r.wheelSpeeds(w,:), '-', ...
                    'Color', wheel_colors(w,:), 'LineWidth', 1.7, ...
                    'DisplayName', sprintf('Wheel %d', w));
            end
        end
        xlabel('Time (s)', 'FontName', fn, 'FontSize', fs);
        ylabel('Output Velocity (m/s)', 'FontName', fn, 'FontSize', fs);
        title('Output Velocity', 'FontName', fn, 'FontSize', fs);
        legend('FontSize', 6, 'FontName', fn);
        grid on;
        set(gca, 'FontName', fn, 'FontSize', fs);
        hold off;

        case_fig = fullfile(per_case_dir, sprintf('seed_%04d.png', r.seed));
        saveas(fig, case_fig);
        close(fig);
    end

    fprintf('[算法 %s] 逐seed轨迹图已保存: %s (%d个)\n', ...
        algorithm, per_case_dir, length(success_results));
end


%% =========================================================
% 单算法 CSV 导出 (与 compare_algorithms 风格一致)
% ==========================================================

function save_algorithm_csv_paper(comparison, alg_idx, results_dir)

    algorithm = comparison.algorithms{alg_idx};
    alg_results = get_alg_results_paper(comparison, alg_idx);
    n = length(alg_results);

    col_names = { ...
        'seed', 'algorithm', 'success', ...
        'rmse', 'trajectoryCost', ...
        'meanSolveTime', 'maxSolveTime', 'totalSolveTime', ...
        'medianSolveTime', 'q1SolveTime', 'q3SolveTime', ...
        'wheelSpeedViolationRatio', 'maxWheelSpeedViolation', ...
        'steeringRateViolationRatio', 'maxSteeringRateViolation', ...
        'successRate', 'validStepRate', 'elapsedTime' ...
    };

    T = table();
    for c = 1:length(col_names)
        if strcmp(col_names{c}, 'algorithm')
            T.(col_names{c}) = cell(n, 1);
        else
            T.(col_names{c}) = zeros(n, 1);
        end
    end

    for i = 1:n
        r = alg_results(i);
        m = r.metrics;

        T.seed(i) = r.seed;
        T.algorithm{i} = algorithm;
        T.success(i) = r.success;

        fnames = {'rmse', 'trajectoryCost', ...
                  'meanSolveTime', 'maxSolveTime', 'totalSolveTime', ...
                  'medianSolveTime', 'q1SolveTime', 'q3SolveTime', ...
                  'wheelSpeedViolationRatio', 'maxWheelSpeedViolation', ...
                  'steeringRateViolationRatio', 'maxSteeringRateViolation', ...
                  'successRate', 'validStepRate'};
        for f = 1:length(fnames)
            if isfield(m, fnames{f}) && ~isempty(m.(fnames{f})) && isnumeric(m.(fnames{f}))
                T.(fnames{f})(i) = m.(fnames{f});
            else
                T.(fnames{f})(i) = NaN;
            end
        end

        T.elapsedTime(i) = r.elapsedTime;
    end

    csv_file = fullfile(results_dir, sprintf('%s_results.csv', algorithm));

    try
        writetable(T, csv_file);
        fprintf('[算法 %s] CSV 已保存: %s\n', algorithm, csv_file);
    catch
        csv_tmp = fullfile(results_dir, sprintf('%s_results_new.csv', algorithm));
        writetable(T, csv_tmp);
        fprintf('[算法 %s] 原CSV被占用, 已保存到: %s\n', algorithm, csv_tmp);
    end
end


%% =========================================================
% Table II 格式汇总
% ==========================================================

function print_table_ii_paper(comparison)

    algorithms = comparison.algorithms;
    num_algorithms = length(algorithms);

    fprintf('\n');
    fprintf('====================================================================\n');
    fprintf('       Table II: Paper Reproduction Summary (P1-3/P1-4 enhanced)\n');
    fprintf('====================================================================\n');
    fprintf('| %-15s | %5s | %8s | %8s | %8s | %5s | %11s |\n', ...
        'Algorithm', 'Valid', 'J_total', 'Med(ms)', 'Mean(ms)', 'VRate', 'Constraints');
    fprintf('|-----------------|-------|----------|----------|----------|-------|-------------|\n');

    for alg_idx = 1:num_algorithms
        alg = algorithms{alg_idx};
        alg_results = get_alg_results_paper(comparison, alg_idx);
        n_success = sum([alg_results.success]);

        if n_success > 0
            success_mask = [alg_results.success];
            success_results = alg_results(success_mask);

            traj_costs = arrayfun(@(r) r.metrics.trajectoryCost, success_results);
            mean_J = nanmean(traj_costs);

            % P1-4: median solve time (warm-up excluded)
            if isfield(success_results(1).metrics, 'medianSolveTime')
                med_times = arrayfun(@(r) r.metrics.medianSolveTime, success_results) * 1000;
                med_time = nanmean(med_times);
            else
                med_time = NaN;
            end

            mean_times = arrayfun(@(r) safe_metric(r, 'meanSolveTime'), success_results) * 1000;
            mean_time = nanmean(mean_times);

            ws_viol = nanmean(arrayfun(@(r) safe_metric(r, 'wheelSpeedViolationRatio'), success_results));
            sr_viol = nanmean(arrayfun(@(r) safe_metric(r, 'steeringRateViolationRatio'), success_results));

            % P1-3: valid step rate
            if isfield(success_results(1), 'solverInfo') && isfield(success_results(1).solverInfo, 'nValidSteps')
                vrate = nanmean(arrayfun(@(r) r.metrics.validStepRate, success_results));
                feas_str = sprintf('%d/1', n_success);
            else
                vrate = NaN;
                feas_str = sprintf('%d/1', n_success);
            end

            J_str = sprintf('%.2f', mean_J);

            if isnan(med_time)
                med_str = 'N/A';
            else
                med_str = sprintf('%.2f', med_time);
            end

            time_str = sprintf('%.1f', mean_time);

            if isnan(vrate)
                vrate_str = 'N/A';
            else
                vrate_str = sprintf('%.0f%%', vrate * 100);
            end

            if ws_viol > 0.01 || sr_viol > 0.01
                constr_str = 'Violated';
            else
                constr_str = 'Satisfied';
            end
        else
            feas_str = '0/1';
            J_str = 'N/A';
            med_str = 'N/A';
            time_str = 'N/A';
            vrate_str = 'N/A';
            constr_str = 'N/A';
        end

        fprintf('| %-15s | %5s | %8s | %8s | %8s | %5s | %11s |\n', ...
            alg, feas_str, J_str, med_str, time_str, vrate_str, constr_str);
    end

    fprintf('====================================================================\n');
    fprintf('Valid = 成功完成; Med = warm-up后中位数耗时; VRate = 有效步比例\n');
    fprintf('有效步 = exitflag>=0 且解有限\n');
end


%% =========================================================
% 辅助函数
% ==========================================================

function alg_results = get_alg_results_paper(comparison, alg_idx)
    % 从 cell 数组中安全提取某算法的所有 results 为 struct 数组
    algorithm = comparison.algorithms{alg_idx};
    alg_mask = strcmp(comparison.completedPairs(:, 2), algorithm);
    n_match = sum(alg_mask);
    if n_match == 0
        alg_results = struct([]);
    elseif n_match == 1
        alg_results = comparison.results{find(alg_mask, 1)};
    else
        alg_results = [comparison.results{alg_mask}];
    end
end


function val = safe_metric(r, field_name)
    if isfield(r, 'metrics') && isfield(r.metrics, field_name)
        val = r.metrics.(field_name);
        if ~isnumeric(val), val = NaN; end
    else
        val = NaN;
    end
end
