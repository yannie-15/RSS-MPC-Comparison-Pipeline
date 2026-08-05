function plot_paper_comparison(comparison, results_dir)
% PLOT_PAPER_COMPARISON 生成论文风格的算法对比图
%
% 输入:
%   comparison : compare_algorithms 输出的结构体
%   results_dir: 结果保存目录
%
% 输出:
%   在 results_dir/comparison/ 下保存:
%     - fig3_trajectories.png : 轨迹跟踪对比图 (2x2, 论文 Fig.3 风格)
%     - fig4_solve_time.png   : 求解时间 box plot
%     - fig5_constraints.png  : 约束满足对比图
%
% 用法:
%   plot_paper_comparison(comparison, 'D:\Projects\RSS\results');

    if nargin < 2 || isempty(results_dir)
        results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'batch');
    end

    algorithms = comparison.algorithms;
    num_algorithms = length(algorithms);

    fig_dir = fullfile(results_dir, 'comparison');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

    % 检查每种算法是否有成功结果
    alg_has_success = false(1, num_algorithms);
    for alg_idx = 1:num_algorithms
        alg_results = get_alg_results(comparison, alg_idx);
        if ~isempty(alg_results)
            alg_has_success(alg_idx) = any([alg_results.success]);
        end
    end

    viable_indices = find(alg_has_success);

    if isempty(viable_indices)
        fprintf('所有算法都失败了，跳过对比图\n');
        return;
    end

    %% ---- Fig.3: 轨迹跟踪对比 (论文风格 2x2) ----
    plot_fig3_trajectories(comparison, algorithms, viable_indices, fig_dir);

    %% ---- Fig.4: 求解时间 box plot ----
    plot_fig4_solve_time(comparison, algorithms, viable_indices, fig_dir);

    %% ---- Fig.5: 约束满足对比 ----
    plot_fig5_constraints(comparison, algorithms, viable_indices, fig_dir);

    fprintf('论文风格对比图已保存到: %s\n', fig_dir);
end


%% =========================================================
% Fig.3: 轨迹跟踪性能 (2x2 论文风格)
%% =========================================================

function plot_fig3_trajectories(comparison, algorithms, viable_indices, fig_dir)

    n_viable = length(viable_indices);
    if n_viable == 0, return; end

    % 标签名称 (按论文图顺序)
    alg_labels = {'Active-set Algorithm', 'Interior-point Algorithm', 'e-LMPC', 'Proposed Approach with 3 Iterations'};
    subplot_tags = {'(a)', '(b)', '(c)', '(d)'};

    % 颜色
    ref_color = [0.2, 0.2, 0.55];   % 深蓝 (参考轨迹)
    actual_color = [0.75, 0.15, 0.05]; % 暗红 (实际轨迹 x 标记)

    % 字体设置
    font_name = 'Times New Roman';
    label_font_size = 12;
    title_font_size = 13;
    tag_font_size = 14;

    % 布局: 4个算法 → 2x2; 否则 1xn
    if n_viable == 4
        n_rows = 2; n_cols = 2;
        fig_width = 900; fig_height = 900;  % 给底部标签留空间
    else
        n_rows = 1; n_cols = n_viable;
        fig_width = 450 * n_viable; fig_height = 450;
    end

    fig = figure('Position', [50, 50, fig_width, fig_height], 'Visible', 'off');

    for idx = 1:n_viable
        alg_idx = viable_indices(idx);
        alg_results = get_alg_results(comparison, alg_idx);

        success_mask = [alg_results.success];
        success_idx = find(success_mask, 1);
        if isempty(success_idx), continue; end

        r = alg_results(success_idx);

        if ~isfield(r, 'states') || isempty(r.states) || ~isfield(r, 'path') || isempty(r.path)
            continue;
        end

        ax = subplot(n_rows, n_cols, idx);
        hold on;

        % 参考轨迹: 蓝色实线
        plot(r.path(1,:), r.path(2,:), '-', 'Color', ref_color, 'LineWidth', 2.0);

        % 实际轨迹: 红色 × 标记
        plot(r.states(1,:), r.states(2,:), 'x', 'Color', actual_color, ...
             'MarkerSize', 5, 'LineWidth', 1.2);

        % 标题
        title('Trajectory Tracking Performance', ...
              'FontName', font_name, 'FontSize', title_font_size, 'Color', 'k');

        % 坐标轴标签
        xlabel('x_w (m)', 'FontName', font_name, 'FontSize', label_font_size, 'Color', 'k');
        ylabel('y_w (m)', 'FontName', font_name, 'FontSize', label_font_size, 'Color', 'k');

        % 刻度字体
        set(gca, 'FontName', font_name, 'FontSize', label_font_size, 'Color', 'k');

        % 网格
        grid on;
        set(gca, 'GridLineStyle', ':', 'GridColor', [0.5 0.5 0.5], 'GridAlpha', 0.5);

        % 边框颜色
        ax.XColor = 'k';
        ax.YColor = 'k';
        ax.Box = 'on';

        hold off;

        % 子图标签 (放在图下方，论文风格)
        alg_name = alg_labels{min(alg_idx, length(alg_labels))};
        tag_str = sprintf('%s %s', subplot_tags{min(alg_idx, length(subplot_tags))}, alg_name);

        % 使用 annotation 把标签放在子图下方中心
        pos = ax.Position;
        annotation('textbox', [pos(1), pos(2) - 0.06, pos(3), 0.05], ...
                   'String', tag_str, ...
                   'FontName', font_name, 'FontSize', tag_font_size, ...
                   'Color', 'k', 'HorizontalAlignment', 'center', ...
                   'VerticalAlignment', 'top', 'EdgeColor', 'none');
    end

    fig_file = fullfile(fig_dir, 'fig3_trajectories.png');
    saveas(fig, fig_file);
    close(fig);

    fprintf('[Fig.3] 轨迹跟踪对比图已保存: %s\n', fig_file);
end


%% =========================================================
% Fig.4: 求解时间 box plot
%% =========================================================

function plot_fig4_solve_time(comparison, algorithms, viable_indices, fig_dir)

    solve_time_data = [];
    alg_labels = {};

    for idx = 1:length(viable_indices)
        alg_idx = viable_indices(idx);
        alg_results = get_alg_results(comparison, alg_idx);
        success_mask = [alg_results.success];

        times_ms = arrayfun(@(r) safe_metric(r, 'meanSolveTime'), alg_results) * 1000;
        valid_times = times_ms(success_mask);
        valid_times = valid_times(~isnan(valid_times));

        if ~isempty(valid_times)
            solve_time_data = [solve_time_data; valid_times]; %#ok<AGROW>
            alg_labels = [alg_labels; repmat({algorithms{alg_idx}}, length(valid_times), 1)]; %#ok<AGROW>
        end
    end

    if isempty(solve_time_data)
        fprintf('[Fig.4] 无有效求解时间数据，跳过\n');
        return;
    end

    fig = figure('Position', [50, 50, 600, 400], 'Visible', 'off');

    % boxplot 属于 Statistics Toolbox, 协作者可能未安装; 回退到基础 MATLAB 的 boxchart (R2020a+)
    if exist('boxplot', 'file')
        boxplot(solve_time_data, alg_labels);
    else
        fprintf('[Fig.4] Statistics Toolbox 不可用, 改用 boxchart (基础 MATLAB)\n');
        alg_cat = categorical(alg_labels);
        boxchart(alg_cat, solve_time_data);
    end

    font_name = 'Times New Roman';
    xlabel('Algorithm', 'FontName', font_name, 'FontSize', 12, 'Color', 'k');
    ylabel('Mean solve time per MPC step (ms)', 'FontName', font_name, 'FontSize', 12, 'Color', 'k');
    title('Fig.4: Solve time comparison', 'FontName', font_name, 'FontSize', 13, 'Color', 'k');
    set(gca, 'FontName', font_name, 'FontSize', 11, 'XColor', 'k', 'YColor', 'k');
    grid on;

    fig_file = fullfile(fig_dir, 'fig4_solve_time.png');
    saveas(fig, fig_file);
    close(fig);

    fprintf('[Fig.4] 求解时间对比图已保存: %s\n', fig_file);
end


%% =========================================================
% Fig.5: 约束满足对比
%% =========================================================

function plot_fig5_constraints(comparison, algorithms, viable_indices, fig_dir)

    n_viable = length(viable_indices);
    if n_viable == 0, return; end

    fig = figure('Position', [50, 50, 1000, 200 * n_viable], 'Visible', 'off');
    colors = lines(n_viable);
    font_name = 'Times New Roman';

    for idx = 1:n_viable
        alg_idx = viable_indices(idx);
        alg_results = get_alg_results(comparison, alg_idx);
        success_mask = [alg_results.success];

        violation_ratios = arrayfun(@(r) safe_metric(r, 'wheelSpeedViolationRatio'), alg_results);
        violation_ratios(~success_mask) = NaN;
        steer_ratios = arrayfun(@(r) safe_metric(r, 'steeringRateViolationRatio'), alg_results);
        steer_ratios(~success_mask) = NaN;

        subplot(n_viable, 2, 2*idx-1);
        valid_viol = violation_ratios(success_mask);
        bar(valid_viol, 'FaceColor', colors(idx, :));
        ylabel('Wheel speed violation ratio', 'FontName', font_name, 'FontSize', 11, 'Color', 'k');
        title(sprintf('%s: Wheel speed', algorithms{alg_idx}), 'FontName', font_name, 'FontSize', 12, 'Color', 'k');
        set(gca, 'FontName', font_name, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        grid on;

        subplot(n_viable, 2, 2*idx);
        valid_steer = steer_ratios(success_mask);
        bar(valid_steer, 'FaceColor', colors(idx, :));
        ylabel('Steering rate violation ratio', 'FontName', font_name, 'FontSize', 11, 'Color', 'k');
        title(sprintf('%s: Steering rate', algorithms{alg_idx}), 'FontName', font_name, 'FontSize', 12, 'Color', 'k');
        set(gca, 'FontName', font_name, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k');
        grid on;
    end

    fig_file = fullfile(fig_dir, 'fig5_constraints.png');
    saveas(fig, fig_file);
    close(fig);

    fprintf('[Fig.5] 约束满足对比图已保存: %s\n', fig_file);
end


%% =========================================================
% 辅助函数
%% =========================================================

function alg_results = get_alg_results(comparison, alg_idx)
    % 从 cell 数组中安全提取某算法的所有 results 为 struct 数组
    if ~isfield(comparison, 'completedPairs') || ~isfield(comparison, 'results')
        alg_results = struct([]);
        return;
    end

    % completedPairs(:,2) 是 cell 数组，存的是算法名字符串
    if alg_idx >= 1 && alg_idx <= length(comparison.algorithms)
        algorithm = comparison.algorithms{alg_idx};
    else
        alg_results = struct([]);
        return;
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
