function plot_one_algorithm(comparison, alg_idx, results_dir)
%PLOT_ONE_ALGORITHM  为单个算法生成汇总图 + 逐 seed 轨迹图
%
% 输入:
%   comparison : comparison 结构体
%   alg_idx    : 算法索引
%   results_dir: 结果根目录 (函数会在其下创建 per_algorithm/{算法名}/)
%
% 输出:
%   - {results_dir}/per_algorithm/{算法名}_summary.png
%   - {results_dir}/per_algorithm/{算法名}/seed_XXXX.png  (逐 seed, 通过 replot_per_seed 单独生成)
%
% 注: 逐 seed 轨迹图的批量生成在 print 时容易卡住,
%     此处仅生成 summary 图; 逐 seed 图由 replot_per_seed.m 单独处理。

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
    print(fig, fig_file, '-dpng', '-r150');
    close(fig);

    fprintf('[算法 %s] 汇总图已保存: %s\n', algorithm, fig_file);

    % ===== 逐 seed 轨迹图 =====
    % 注意: 批量生成逐seed图可能导致 print 卡住, 使用独立的 replot_per_seed.m 手动生成
    return;  % 跳过逐seed图, 由 replot_per_seed.m 单独处理
end
