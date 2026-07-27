function replot_per_seed(results_file)
% REPLOT_PER_SEED 重新生成所有算法的逐 seed 仿真轨迹图
%
% 风格与旧版 main.m 实时动画绘图一致:
%   子图1: 轨迹跟踪性能 (参考轨迹 + 实际轨迹)
%   子图2: 各轮转向速率 (含约束线)
%   子图3: 各轮输出速度 (含约束线)
%
% 用法:
%   replot_per_seed()                                    % 使用默认路径
%   replot_per_seed('D:\Projects\RSS\results\comparison_checkpoint.mat')
%   replot_per_seed('D:\Projects\RSS\results\comparison_final.mat')
%
% 输出:
%   results/per_algorithm/{算法名}/seed_XXXX.png

    if nargin < 1 || isempty(results_file)
        script_dir = fileparts(mfilename('fullpath'));
        workspace_root = fileparts(script_dir);
        results_dir = fullfile(workspace_root, 'results');

        % 优先使用 final, 其次 checkpoint
        final_file = fullfile(results_dir, 'comparison_final.mat');
        checkpoint_file = fullfile(results_dir, 'comparison_checkpoint.mat');

        if exist(final_file, 'file')
            results_file = final_file;
        elseif exist(checkpoint_file, 'file')
            results_file = checkpoint_file;
        else
            error('未找到结果文件, 请指定路径');
        end
    end

    fprintf('加载结果: %s\n', results_file);
    loaded = load(results_file, 'comparison');
    comparison = loaded.comparison;

    script_dir = fileparts(mfilename('fullpath'));
    workspace_root = fileparts(script_dir);
    results_dir = fullfile(workspace_root, 'results');

    algorithms = comparison.algorithms;
    num_algorithms = length(algorithms);

    for alg_idx = 1:num_algorithms
        algorithm = algorithms{alg_idx};
        alg_results = get_alg_results(comparison, alg_idx);

        n_success = sum([alg_results.success]);
        fprintf('[%s] 成功: %d\n', algorithm, n_success);

        if n_success == 0
            continue;
        end

        fig_dir = fullfile(results_dir, 'per_algorithm', algorithm);
        if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

        success_mask = [alg_results.success];
        success_results = alg_results(success_mask);

        for i = 1:length(success_results)
            r = success_results(i);
            plot_one_case(r, fig_dir);
        end

        fprintf('[%s] 已生成 %d 张逐seed图\n', algorithm, length(success_results));
    end

    fprintf('\n全部完成!\n');
end


%% =========================================================
% 绘制单个 seed 的仿真结果 (匹配 main.m 风格)
% ==========================================================

function plot_one_case(r, fig_dir)

    if ~isfield(r, 'states') || isempty(r.states), return; end
    if ~isfield(r, 'path') || isempty(r.path), return; end

    num_w = r.numWheels;
    dt = r.config.dt;
    num_steps = size(r.states, 2) - 1;
    t_vec = (0:num_steps-1) * dt;

    % 颜色定义 (与 slanCL / lines() 一致)
    wheel_colors = lines(max(num_w, 5));
    ref_color  = wheel_colors(4, :);   % 深蓝 (参考轨迹)
    act_color  = wheel_colors(2, :);   % 橙红 (实际轨迹)
    cons_color = wheel_colors(5, :);   % 浅蓝 (约束线)

    % 字体
    fn = 'Times New Roman';
    fs = 11;

    fig = figure('Position', [100, 100, 1500, 300], 'Visible', 'off');

    % ===== 子图1: 轨迹跟踪性能 =====
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

    % ===== 子图2: 各轮转向速率 =====
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
                'Color', wheel_colors(w, :), 'LineWidth', 1.7, ...
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

    % ===== 子图3: 各轮输出速度 =====
    subplot(1, 3, 3);
    hold on;
    yline(r.config.vimax, '--', 'Color', cons_color, ...
        'LineWidth', 2, 'DisplayName', 'Constraints');

    if isfield(r, 'wheelSpeeds') && ~isempty(r.wheelSpeeds)
        for w = 1:num_w
            plot(t_vec, r.wheelSpeeds(w,:), '-', ...
                'Color', wheel_colors(w, :), 'LineWidth', 1.7, ...
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

    % 保存
    fig_file = fullfile(fig_dir, sprintf('seed_%04d.png', r.seed));
    saveas(fig, fig_file);
    close(fig);
end


%% =========================================================
% 辅助函数
% ==========================================================

function alg_results = get_alg_results(comparison, alg_idx)
    if ~isfield(comparison, 'completedPairs') || ~isfield(comparison, 'results')
        alg_results = struct([]);
        return;
    end

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
