function comparison = compare_algorithms(seed_range, algorithms, force_regen)
% COMPARE_ALGORITHMS  四算法批量对比 + per-(seed,algorithm) 增量续跑
%
% 用法:
%   comparison = compare_algorithms(1:10, {'proposed-3iter', 'e-lmpc', 'active-set', 'interior-point'})
%   comparison = compare_algorithms(1:100, {'proposed-3iter', 'e-lmpc', 'active-set', 'interior-point'})
%   comparison = compare_algorithms(1:50, {'proposed-3iter'})  % 只跑一种算法
%   comparison = compare_algorithms(1:200, {'proposed-3iter', 'e-lmpc', ...}, true)  % 强制重新生成场景
%
% 核心设计:
%   1. 每个 seed = 1 个独立可复现场景
%   2. 所有算法面对同一个 seed 时，场景完全相同
%   3. 每完成 1 个 (seed, algorithm) 组合 → 立即保存 checkpoint
%   4. 中断后重新运行 → 自动跳过已完成的组合，从断点续跑
%   5. 扩大 seed_range（如 1:50 → 1:200）→ 只跑新增的 seed
%
% 输入:
%   seed_range   : 种子范围，如 1:10, 1:100
%   algorithms   : 算法列表（默认 4 种，不含 proposed-1iter）
%   force_regen  : 是否强制重新生成场景文件（默认 false）

    if nargin < 1 || isempty(seed_range), seed_range = 1:10; end
    if nargin < 2 || isempty(algorithms)
        algorithms = {'proposed-3iter', 'e-lmpc', 'active-set', 'interior-point'};
    end
    if nargin < 3 || isempty(force_regen), force_regen = false; end

    % 自动确保项目路径已添加 (避免用户忘记运行 setup_paths)
    script_dir = fileparts(mfilename('fullpath'));
    workspace_root = fileparts(script_dir);
    addpath(fullfile(workspace_root, 'matlab'));
    addpath(fullfile(workspace_root, 'RSS_proposed'));

    seeds = seed_range(:)';  % 确保为行向量
    num_seeds = length(seeds);
    num_algorithms = length(algorithms);

    fprintf('========================================\n');
    fprintf('四算法自动对比批量仿真\n');
    fprintf('种子范围: %d-%d (共%d个)\n', seeds(1), seeds(end), num_seeds);
    fprintf('算法: %s\n', strjoin(algorithms, ', '));
    fprintf('总组合数: %d\n', num_seeds * num_algorithms);
    fprintf('========================================\n');

    %% =====================================================
    % 1. 确定保存路径和断点续跑
    % ======================================================

    results_dir = fullfile(workspace_root, 'results');

    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    checkpoint_file = fullfile(results_dir, 'comparison_checkpoint.mat');

    %% =====================================================
    % 2. 检查断点续跑
    % ======================================================

    if exist(checkpoint_file, 'file') && ~force_regen
        fprintf('发现已有结果文件: %s\n', checkpoint_file);
        loaded = load(checkpoint_file, 'comparison');
        comparison = loaded.comparison;

        % 兼容旧格式: 旧版 completedPairs 是数值矩阵 [seed, alg_idx]
        % 需要转换为 cell 数组 {seed, algorithm_name}
        if ~iscell(comparison.completedPairs)
            old_pairs = comparison.completedPairs;
            new_pairs = cell(size(old_pairs, 1), 2);
            for r = 1:size(old_pairs, 1)
                new_pairs(r, 1) = {old_pairs(r, 1)};
                old_idx = old_pairs(r, 2);
                if old_idx <= length(comparison.algorithms)
                    new_pairs(r, 2) = {comparison.algorithms{old_idx}};
                else
                    new_pairs(r, 2) = {'unknown'};
                end
            end
            comparison.completedPairs = new_pairs;
        end

        % 更新元数据以适配扩大范围 (保留已完成数据)
        comparison.seeds = seeds;
        comparison.algorithms = algorithms;
        comparison.numSeeds = length(seeds);
        comparison.numAlgorithms = length(algorithms);

        % 统计已完成数
        n_completed = size(comparison.completedPairs, 1);
        total_pairs = num_seeds * num_algorithms;

        fprintf('断点续跑: 已完成 %d/%d 组合\n', n_completed, total_pairs);
    else
        comparison = init_comparison_struct(seeds, algorithms);
    end

    %% =====================================================
    % 3. 逐 seed × 逐算法运行 (增量保存)
    % ======================================================

    total_tic = tic;
    n_completed_before = size(comparison.completedPairs, 1);

    for alg_idx = 1:num_algorithms
        algorithm = algorithms{alg_idx};

        % 清除 e-LMPC 的 warm start (切换算法时需要重置)
        clear functions;

        for seed_idx = 1:num_seeds
            seed_val = seeds(seed_idx);

            % 检查此 (seed, algorithm) 组合是否已成功完成
            % 只跳过成功的; 失败的重新跑
            if is_pair_completed(comparison.completedPairs, comparison.results, seed_val, algorithm)
                continue;
            end

            case_tic = tic;

            % 加载/生成此 seed 的场景
            [cfg, scen] = scenario_bank(seed_val, force_regen);

            % 设置 config 的算法标签
            cfg.algorithm = algorithm;

            % 对于 proposed 算法, 使用 HPIPM 求解器
            if contains(algorithm, 'proposed')
                cfg.use_hpipm = true;
            end

            % 关闭实时画图
            cfg.output.livePlot = false;
            cfg.output.saveFigures = false;

            fprintf('\n--- [%s] seed=%d (%d/%d seeds) ---\n', algorithm, seed_val, seed_idx, num_seeds);

            try
                summary = run_one_case(cfg, scen);
                summary.elapsedTime = toc(case_tic);
            catch ME
                summary = struct();
                summary.success = false;
                summary.failureReason = getReport(ME);
                summary.scenario = scen;
                summary.config = cfg;
                summary.algorithm = algorithm;
                summary.seed = seed_val;
                summary.elapsedTime = toc(case_tic);

                if isfield(summary, 'metrics')
                    summary.metrics.rmse = NaN;
                    summary.metrics.trajectoryCost = NaN;
                    summary.metrics.successRate = 0;
                else
                    summary.metrics = struct('rmse', NaN, 'trajectoryCost', NaN, 'successRate', 0);
                end

                fprintf('[seed=%d: %s] EXCEPTION: %s\n', seed_val, algorithm, ME.message);
                for si = 1:length(ME.stack)
                    fprintf('  at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
                end
            end

            % 记录 seed
            summary.seed = seed_val;

            % ========== 添加到 results 和 completedPairs ==========
            % 如果是重跑失败的组合，替换旧记录; 否则追加新记录
            existing_idx = find_pair_index(comparison.completedPairs, seed_val, algorithm);
            if ~isempty(existing_idx)
                comparison.results{existing_idx} = summary;
                % completedPairs 不需要更新 (seed_val, algorithm 已存在)
            else
                comparison.results{end+1} = summary;
                comparison.completedPairs(end+1, :) = {seed_val, algorithm};
            end

            % ========== 增量保存 ==========
            save_checkpoint(comparison, checkpoint_file);

            n_done = size(comparison.completedPairs, 1);
            fprintf('  → 已保存进度 (完成 %d/%d 组合)\n', n_done, num_seeds * num_algorithms);
        end

        % ========== 每种算法所有 seed 完成后，生成该算法的汇总 ==========
        alg_names = comparison.completedPairs(:, 2);
        n_alg_done = sum(strcmp(alg_names, algorithm));
        if n_alg_done >= num_seeds
            fprintf('[算法 %s] 全部 seed 完成，生成图像...\n', algorithm);
            plot_one_algorithm(comparison, alg_idx, results_dir);
            save_algorithm_csv(comparison, alg_idx, results_dir);
        end
    end

    total_elapsed = toc(total_tic);

    %% =====================================================
    % 4. 汇总统计
    % ======================================================

    comparison.totalElapsed = total_elapsed;

    for alg_idx = 1:num_algorithms
        algorithm = algorithms{alg_idx};
        alg_results = get_alg_results(comparison, alg_idx);
        n_success = sum([alg_results.success]);

        fprintf('\n[算法 %s] 成功: %d/%d, 失败: %d/%d\n', ...
            algorithm, n_success, num_seeds, num_seeds - n_success, num_seeds);
    end

    %% =====================================================
    % 5. 保存完整结果 + 生成对比图
    % ======================================================

    final_file = fullfile(results_dir, 'comparison_final.mat');
    save(final_file, 'comparison', '-v7.3');
    fprintf('\n完整结果已保存到: %s\n', final_file);

    % 同步更新 checkpoint 文件 (保留以供下次扩大范围续跑)
    save_checkpoint(comparison, checkpoint_file);
    fprintf('Checkpoint 已更新 (下次扩大范围可续跑)\n');

    % 生成论文风格对比图 (调用独立文件 plot_paper_comparison.m)
    plot_paper_comparison(comparison, results_dir);

    % 生成 Table II 格式汇总
    print_table_ii(comparison);

end


%% =========================================================
% 初始化 comparison 结构体
%% ==========================================================

function comparison = init_comparison_struct(seeds, algorithms)
    comparison = struct();
    comparison.seeds = seeds;
    comparison.algorithms = algorithms;
    comparison.numSeeds = length(seeds);
    comparison.numAlgorithms = length(algorithms);

    % results: 所有已完成的结果（cell 数组，每个 cell 是一个 summary struct）
    comparison.results = {};

    % completedPairs: 已完成的 (seed_value, algorithm_name) 记录
    % cell 数组, 每行 = {seed_val, algorithm_name}
    comparison.completedPairs = cell(0, 2);

    comparison.totalElapsed = 0;
end


%% =========================================================
% 判断 (seed, algorithm) 组合是否已完成
%% ==========================================================

function completed = is_pair_completed(completedPairs, results, seed_val, algorithm)
    % 查找 completedPairs 中是否有此 (seed_val, algorithm) 行且结果为成功
    % 只跳过成功的组合; 失败的允许重跑
    if isempty(completedPairs)
        completed = false;
        return;
    end
    seed_match = cell2mat(completedPairs(:, 1)) == seed_val;
    alg_match = strcmp(completedPairs(:, 2), algorithm);
    pair_match = seed_match & alg_match;

    if ~any(pair_match)
        completed = false;
        return;
    end

    % 检查该组合的结果是否成功
    idx = find(pair_match, 1);
    if idx <= length(results) && isfield(results{idx}, 'success') && results{idx}.success
        completed = true;
    else
        completed = false;  % 失败或结果缺失 → 允许重跑
    end
end


function idx = find_pair_index(completedPairs, seed_val, algorithm)
    % 查找 completedPairs 中 (seed_val, algorithm) 的行号
    if isempty(completedPairs)
        idx = [];
        return;
    end
    seed_match = cell2mat(completedPairs(:, 1)) == seed_val;
    alg_match = strcmp(completedPairs(:, 2), algorithm);
    idx = find(seed_match & alg_match, 1);
end


%% =========================================================
% 增量保存断点文件
%% ==========================================================

function save_checkpoint(comparison, checkpoint_file)
    save(checkpoint_file, 'comparison', '-v7.3');
end


%% =========================================================
% 单算法图像生成
%% ==========================================================

function plot_one_algorithm(comparison, alg_idx, results_dir)

    algorithm = comparison.algorithms{alg_idx};
    alg_results = get_alg_results(comparison, alg_idx);

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

    fprintf('[算法 %s] 逐seed轨迹图已保存: %s (%d个)\n', algorithm, per_case_dir, length(success_results));
end


%% =========================================================
% 单算法 CSV 导出
%% ==========================================================

function save_algorithm_csv(comparison, alg_idx, results_dir)

    algorithm = comparison.algorithms{alg_idx};
    alg_results = get_alg_results(comparison, alg_idx);
    n = length(alg_results);

    col_names = { ...
        'seed', 'algorithm', 'success', ...
        'rmse', 'trajectoryCost', ...
        'meanSolveTime', 'maxSolveTime', 'totalSolveTime', ...
        'wheelSpeedViolationRatio', 'maxWheelSpeedViolation', ...
        'steeringRateViolationRatio', 'maxSteeringRateViolation', ...
        'successRate', 'elapsedTime' ...
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
                  'wheelSpeedViolationRatio', 'maxWheelSpeedViolation', ...
                  'steeringRateViolationRatio', 'maxSteeringRateViolation', ...
                  'successRate'};
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

    % 尝试写入 CSV, 如果文件被占用则写到临时文件
    try
        writetable(T, csv_file);
        fprintf('[算法 %s] CSV 已保存: %s\n', algorithm, csv_file);
    catch
        % 文件被 Excel 等占用, 写到临时文件
        csv_tmp = fullfile(results_dir, sprintf('%s_results_new.csv', algorithm));
        writetable(T, csv_tmp);
        fprintf('[算法 %s] 原CSV被占用, 已保存到: %s\n', algorithm, csv_tmp);
        fprintf('[提示] 关闭Excel后, 可手动将 _new.csv 重命名覆盖原文件\n');
    end
end


%% =========================================================
% Table II 格式汇总
%% ==========================================================

function print_table_ii(comparison)

    algorithms = comparison.algorithms;
    num_algorithms = length(algorithms);

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('              Table II: Algorithm Comparison Summary\n');
    fprintf('============================================================\n');
    fprintf('| Algorithm       | Feasible | J_total  | Time(ms) | Constraints |\n');
    fprintf('|-----------------|----------|----------|----------|-------------|\n');

    for alg_idx = 1:num_algorithms
        alg = algorithms{alg_idx};
        alg_results = get_alg_results(comparison, alg_idx);
        n_success = sum([alg_results.success]);

        if n_success > 0
            success_mask = [alg_results.success];
            success_results = alg_results(success_mask);

            traj_costs = arrayfun(@(r) r.metrics.trajectoryCost, success_results);
            mean_J = nanmean(traj_costs);

            solve_times = arrayfun(@(r) r.metrics.meanSolveTime, success_results) * 1000;
            mean_time = nanmean(solve_times);

            ws_viol = nanmean(arrayfun(@(r) r.metrics.wheelSpeedViolationRatio, success_results));
            sr_viol = nanmean(arrayfun(@(r) r.metrics.steeringRateViolationRatio, success_results));

            feas_str = sprintf('%d/%d', n_success, comparison.numSeeds);
            J_str = sprintf('%.2f', mean_J);
            time_str = sprintf('%.1f', mean_time);

            if ws_viol > 0.01 || sr_viol > 0.01
                constr_str = 'Violated';
            else
                constr_str = 'Satisfied';
            end
        else
            feas_str = sprintf('0/%d', comparison.numSeeds);
            J_str = 'N/A';
            time_str = 'N/A';
            constr_str = 'N/A';
        end

        fprintf('| %-15s | %8s | %8s | %8s | %11s |\n', ...
            alg, feas_str, J_str, time_str, constr_str);
    end

    fprintf('============================================================\n');
end


%% =========================================================
% 辅助函数
%% ==========================================================

function alg_results = get_alg_results(comparison, alg_idx_or_name)
    % 从 cell 数组中安全提取某算法的所有 results 为 struct 数组
    % alg_idx_or_name: 可以是索引(int)或算法名字(string)
    if ischar(alg_idx_or_name) || isstring(alg_idx_or_name)
        algorithm = char(alg_idx_or_name);
    else
        algorithm = comparison.algorithms{alg_idx_or_name};
    end

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
