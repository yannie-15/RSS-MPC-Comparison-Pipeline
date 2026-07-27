function [configs, scenarios] = scenario_generator(num_scenarios, method, seed)
% SCENARIO_GENERATOR 生成随机参数场景
%
% 输入：
%   num_scenarios : 场景数量
%   method        : 采样方法，'latin_hypercube'（默认）或 'random'
%   seed          : 随机种子（可选，默认不固定；指定后可复现）
%
% 输出：
%   configs    : 1xnum_scenarios config 结构体数组
%   scenarios  : 1xnum_scenarios scenario 描述结构体数组
%
% 随机化的参数：
%   - 轨迹形状（控制点偏移）
%   - 初始位置偏移
%   - 初始速度
%   - 轮速上限 vimax
%   - 转向速度上限 phidotmax
%   - 仿真时长 t_end
%   - 机器人尺寸 Lx, Ly
%
% [P2修复] 新增 seed 参数用于实验复现

    if nargin < 1 || isempty(num_scenarios)
        num_scenarios = 10;
    end
    if nargin < 2 || isempty(method)
        method = 'latin_hypercube';
    end
    if nargin < 3 || isempty(seed)
        seed = [];  % 不固定种子
    end

    % [P2修复] 固定随机种子用于实验复现
    if ~isempty(seed)
        rng(seed);
        fprintf('随机种子已设置为: %d\n', seed);
    end

    %% =====================================================
    % 参数范围定义
    % ======================================================

    % 轨迹控制点偏移范围（相对于原始控制点的偏移量）
    ctrl_offset_range = 0.3;  % ±0.3m

    % 初始位置偏移范围
    init_pos_range = 0.05;    % ±0.05m

    % 初始速度范围
    init_vel_range = [0.005, 0.05];  % m/s

    % 轮速上限范围
    vimax_range = [3, 8];     % m/s

    % 转向速度上限范围
    phidotmax_range = [3*pi, 8*pi];  % rad/s

    % 仿真时长范围
    t_end_range = [0.5, 2.0];  % s

    % 机器人尺寸范围
    Lx_range = [0.4, 0.9];    % m
    Ly_range = [0.2, 0.5];    % m

    %% =====================================================
    % 采样
    % ======================================================

    % 待采样参数数量
    n_params = 16;  % 5 ctrl_offsets + 2 init_pos + 3 init_vel + vimax + phidotmax + t_end + Lx + Ly + 2 corner

    if strcmp(method, 'latin_hypercube')
        samples = latin_hypercube_sample(n_params, num_scenarios);
    else
        samples = rand(num_scenarios, n_params);
    end

    %% =====================================================
    % 构建场景
    % ======================================================

    % 原始控制点（来自 config.m）
    base_ctrl_pts = [
        0,     0;
        0.750, 0.250;
        0.250, 0.750;
        1.250, -1.000;
        1.000, 0
    ];

    configs = struct([]);
    scenarios = struct([]);

    for i = 1:num_scenarios
        s = samples(i, :);
        idx = 1;

        % 1. 轨迹控制点偏移（5个控制点各2维，但起点和终点不变）
        % 只偏移中间3个控制点
        ctrl_pts = base_ctrl_pts;
        for cp = 2:4
            ctrl_pts(cp, 1) = ctrl_pts(cp, 1) + (s(idx) - 0.5) * 2 * ctrl_offset_range;
            idx = idx + 1;
            ctrl_pts(cp, 2) = ctrl_pts(cp, 2) + (s(idx) - 0.5) * 2 * ctrl_offset_range;
            idx = idx + 1;
        end

        % 2. 初始位置偏移
        init_x = 0.05 + (s(idx) - 0.5) * 2 * init_pos_range; idx = idx + 1;
        init_y = 0.1 + (s(idx) - 0.5) * 2 * init_pos_range; idx = idx + 1;

        % 3. 初始速度
        init_vx = init_vel_range(1) + s(idx) * (init_vel_range(2) - init_vel_range(1)); idx = idx + 1;
        init_vy = init_vel_range(1) + s(idx) * (init_vel_range(2) - init_vel_range(1)); idx = idx + 1;
        init_omega = init_vel_range(1) + s(idx) * (init_vel_range(2) - init_vel_range(1)); idx = idx + 1;

        % 4. 约束参数
        vimax = vimax_range(1) + s(idx) * (vimax_range(2) - vimax_range(1)); idx = idx + 1;
        phidotmax = phidotmax_range(1) + s(idx) * (phidotmax_range(2) - phidotmax_range(1)); idx = idx + 1;

        % 5. 仿真时长
        t_end = t_end_range(1) + s(idx) * (t_end_range(2) - t_end_range(1)); idx = idx + 1;

        % 6. 机器人尺寸
        Lx = Lx_range(1) + s(idx) * (Lx_range(2) - Lx_range(1)); idx = idx + 1;
        Ly = Ly_range(1) + s(idx) * (Ly_range(2) - Ly_range(1)); idx = idx + 1;

        % 构建 config
        cfg = defaultConfig();
        cfg.experimentId = sprintf("batch_%04d", i);
        cfg.caseId = sprintf("case_%04d", i);
        cfg.ctrl_pts = ctrl_pts;
        cfg.vimax = vimax;
        cfg.phidotmax = phidotmax;
        cfg.t_end = t_end;
        cfg.num_steps = round(t_end / cfg.dt);
        cfg.num_path_pts = cfg.num_steps;
        cfg.Lx = Lx;
        cfg.Ly = Ly;
        cfg.a = Lx / 2;
        cfg.b = Ly / 2;
        cfg.wheel_pos = [
            cfg.a,  cfg.b;
            -cfg.a, cfg.b;
            -cfg.a, -cfg.b;
            cfg.a,  -cfg.b
        ];
        % 批量运行时关闭画图
        cfg.output.livePlot = false;
        cfg.output.saveFigures = false;

        configs(i).cfg = cfg;

        % 构建 scenario
        scenarios(i).data = struct( ...
            'name', sprintf('random_%04d', i), ...
            'id', i, ...
            'initialState', [init_x; init_y; 0.2], ...
            'initialVelocity', [init_vx; init_vy; init_omega], ...
            'vimax', vimax, ...
            'phidotmax', phidotmax, ...
            't_end', t_end, ...
            'Lx', Lx, ...
            'Ly', Ly ...
        );
    end

    % 提取为结构体数组
    configs = [configs.cfg];
    scenarios = [scenarios.data];

    fprintf('生成了 %d 个随机场景（方法: %s）\n', num_scenarios, method);
end


%% =========================================================
% 拉丁超方采样
% ==========================================================

function samples = latin_hypercube_sample(n_dims, n_samples)
% 生成 n_samples x n_dims 的拉丁超方采样矩阵，值在 [0, 1]

    samples = zeros(n_samples, n_dims);
    for d = 1:n_dims
        perm = randperm(n_samples);
        for i = 1:n_samples
            samples(i, d) = (perm(i) - 1 + rand()) / n_samples;
        end
    end
end
