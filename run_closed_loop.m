function result = run_closed_loop(cfg)
% RUN_CLOSED_LOOP  RSS proposed 算法完整 100 步闭环仿真
%
% 从 algorithms/RSS_proposed/main.m 无损抽取, 改为函数化接口。
% 保留原始数学公式 (轮角 atan2(vyi,vxi), 状态推进用世界系 new_state_dot)。
%
% 输入:
%   cfg : 配置结构体, 由 config() 合并 JSON 覆盖项生成
%         必须包含 config() 的所有字段, 外加:
%         .live_plot     (logical, 默认 false)
%         .save_figures  (logical, 默认 false)
%         .save_full_log (logical, 默认 false)
%         .output_dir    (string, 默认 '')
%         .case_id       (string, 默认 '')
%         .solver        (string, 默认 'sdpt3')
%
% 输出:
%   result : 结构体, 包含完整仿真数据
%       .reference            : 3xN 参考轨迹
%       .state_history        : Nx3 状态 [x, y, psi]
%       .body_velocity_history: Nx3 车体系速度
%       .world_velocity_history: Nx3 世界系速度 (new_state_dot)
%       .control_history      : Nx3 控制增量 u(:,1)
%       .wheel_speed_history  : Nx4 轮速
%       .wheel_angle_history  : Nx4 轮角
%       .steering_rate_history: Nx4 转向率
%       .solver_time_history  : Nx1 每步求解时间
%       .cvx_status_history   : {Nx1} cell 每步 CVX 状态
%       .cvx_optval_history   : Nx1 每步最优值
%       .objective_history    : Nx1 累计代价 J
%       .time_history         : Nx1 时间戳
%       .trajectory_cost      : scalar 总代价
%       .position_rmse        : scalar 位置 RMSE
%       .max_wheel_speed      : scalar 最大轮速
%       .max_steering_rate    : scalar 最大转向率
%       .completed_steps      : scalar 完成步数
%       .success              : logical 是否成功
%       .failure_reason       : string 失败原因
%       .config               : struct 使用的配置
%       .wall_time            : scalar 总壁钟时间 (s)

    %% 参数默认值
    if nargin < 1 || isempty(cfg)
        cfg = config();
    end

    % 补充控制字段默认值
    if ~isfield(cfg, 'live_plot'),     cfg.live_plot = false;     end
    if ~isfield(cfg, 'save_figures'),  cfg.save_figures = false;  end
    if ~isfield(cfg, 'save_full_log'), cfg.save_full_log = false; end
    if ~isfield(cfg, 'output_dir'),    cfg.output_dir = '';       end
    if ~isfield(cfg, 'case_id'),       cfg.case_id = '';          end
    if ~isfield(cfg, 'solver'),        cfg.solver = 'sdpt3';      end

    wall_tic = tic;

    %% 定位 submodule 路径并 addpath
    script_dir = fileparts(mfilename('fullpath'));
    if isempty(script_dir)
        % 从命令行调用时 mfilename 可能无路径
        script_dir = pwd;
    end
    submodule_dir = fullfile(script_dir, 'algorithms', 'RSS_proposed');

    % 也尝试从 matlab/ 目录定位 (兼容从 matlab/ 调用)
    if ~exist(fullfile(submodule_dir, 'control_RSS.m'), 'file')
        parent = fileparts(script_dir);
        submodule_dir = fullfile(parent, 'algorithms', 'RSS_proposed');
    end

    if ~exist(fullfile(submodule_dir, 'control_RSS.m'), 'file')
        error('run_closed_loop:MissingSubmodule', ...
            '找不到 RSS_proposed submodule: %s', submodule_dir);
    end

    addpath(submodule_dir);
    cleanup_path = onCleanup(@() rmpath(submodule_dir));
    clear functions  % 清除函数缓存, 确保加载 submodule 版本

    %% 参数提取
    num_steps = cfg.num_steps;
    dt = cfg.dt;
    vimax = cfg.vimax;
    phidotmax = cfg.phidotmax;

    %% 生成参考轨迹
    path = bezier_path(cfg.ctrl_pts, cfg.num_path_pts);

    %% 初始状态 (与原 main.m 完全一致)
    q = [0.050, 0.1];           % 初始位置 (1x2)
    psi0 = 0.2;                 % 初始朝向 (rad)
    vx = 0.01; vy = 0.01; omega_b = 0.01;
    last_vel = [vx; vy; omega_b];  % 3x1 车体系速度
    state = [q, psi0];             % 1x3
    prev_phi = zeros(4, 1);
    J = 0;
    u = zeros(3, 1);

    %% 预分配历史数组
    q_history          = zeros(num_steps, 2);
    vi_history         = zeros(num_steps, 4);
    phidot_history     = zeros(num_steps, 4);
    psi_history        = zeros(num_steps, 1);
    t_history          = zeros(num_steps, 1);
    phi_history        = zeros(num_steps, 4);
    u_history          = zeros(num_steps, 3);       % 每步 u(:,1)
    velocity_history   = zeros(num_steps, 3);       % 车体系速度
    state_dot_history  = zeros(num_steps, 3);       % 世界系 new_state_dot
    cvx_optval_history = zeros(num_steps, 1);
    cvx_status_history = cell(num_steps, 1);
    solver_time_history = zeros(num_steps, 1);

    %% 清除全局 solver_time_array (control_RSS 使用)
    global solver_time_array
    solver_time_array = [];

    %% 可选实时绘图
    if cfg.live_plot
        fig = figure('Position', [100, 100, 1500, 300], 'Visible', 'on');
        subplot(1, 3, 1);
        plot(path(1,:), path(2,:), '-', 'LineWidth', 2, 'DisplayName', 'Reference'); hold on;
        h_q = animatedline('Marker', 'x', 'LineStyle', 'none', 'LineWidth', 1.5, 'DisplayName', 'Simulation');
        xlabel('x_w (m)'); ylabel('y_w (m)'); title('Trajectory Tracking');
        legend('FontSize', 6); grid on; axis equal;

        subplot(1, 3, 2); hold on;
        yline(phidotmax, '--', 'LineWidth', 1.5);
        yline(-phidotmax, '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        h_phidot = gobjects(4,1);
        for idx = 1:4
            h_phidot(idx) = animatedline('LineWidth', 1.5, 'DisplayName', ['Wheel ' num2str(idx)]);
        end
        xlabel('Time(s)'); ylabel('Steering rate (rad/s)'); title('Steering Rate');
        legend('FontSize', 6); grid on;

        subplot(1, 3, 3); hold on;
        yline(vimax, '--', 'LineWidth', 1.5);
        h_vi = gobjects(4,1);
        for idx = 1:4
            h_vi(idx) = animatedline('LineWidth', 1.5, 'DisplayName', ['Wheel ' num2str(idx)]);
        end
        xlabel('Time (s)'); ylabel('Velocity (m/s)'); title('Output Velocity');
        legend('FontSize', 6); grid on;
    end

    %% 仿真循环
    success = true;
    failure_reason = '';
    completed_steps = 0;

    for k = 1:num_steps
        try
            t = (k - 1) * dt;
            t_history(k) = t;
            q_history(k, :) = q;
            psi_history(k) = psi0;

            % 跟踪误差代价 (原 main.m 公式)
            e = state' - path(:, k);
            J = J + e' * [30,0,0; 0,30,0; 0,0,1] * e + u(:,1)' * [0.3,0,0; 0,0.3,0; 0,0,0.3] * u(:,1);

            % 调用控制器
            [u, new_state_dot, velocity] = control_RSS(path, k, last_vel, state);

            last_vel = velocity;

            % 记录
            u_history(k, :) = u(:, 1)';
            velocity_history(k, :) = velocity';
            state_dot_history(k, :) = new_state_dot';
            cvx_optval_history(k) = cvx_optval;
            cvx_status_history{k} = cvx_status;

            % 从全局 solver_time_array 提取当前步时间
            global solver_time_array
            if ~isempty(solver_time_array)
                % solver_time_array 索引: 3*step + m - 3, 最后一次迭代是 3*k
                idx = 3*k;
                if idx <= length(solver_time_array)
                    solver_time_history(k) = sum(solver_time_array(3*(k-1)+1:idx));
                end
            end

            % 计算轮速和轮角 (原 main.m 公式: atan2(vyi, vxi))
            phi = zeros(4, 1);
            vi = zeros(4, 1);
            for i = 1:4
                Hj = [1, 0, -cfg.wheel_pos(i,2); 0, 1, cfg.wheel_pos(i,1)];
                zn = Hj * velocity;
                vxi = zn(1);
                vyi = zn(2);
                vi(i) = sqrt(vxi^2 + vyi^2);
                phi(i) = atan2(vyi, vxi);
            end

            vi_history(k, :) = vi';
            phi_history(k, :) = phi';

            % 状态更新 (原 main.m: 用 new_state_dot 世界系)
            vx = new_state_dot(1);
            vy = new_state_dot(2);
            omega_b = new_state_dot(3);
            psi0 = psi0 + omega_b * dt;
            q = q + [vx, vy] * dt;
            state_dot = new_state_dot;
            state = [q, psi0];

            % 实时绘图
            if cfg.live_plot
                addpoints(h_q, q(1), q(2));
                for idx = 1:4
                    addpoints(h_vi(idx), t, vi(idx));
                end
                if k > 1
                    for idx = 1:4
                        d_phi = phi(idx) - prev_phi(idx);
                        d_phi = mod(d_phi + pi, 2*pi) - pi;
                        curr_phidot = d_phi / dt;
                        addpoints(h_phidot(idx), t, curr_phidot);
                    end
                end
                drawnow limitrate;
            end

            prev_phi = phi;
            completed_steps = k;

        catch ME
            success = false;
            failure_reason = getReport(ME, 'extended');
            fprintf('[run_closed_loop] Step %d failed: %s\n', k, ME.message);
            break;
        end
    end

    %% 计算转向率 (原 main.m 方式)
    if completed_steps > 1
        phidot_raw = diff(phi_history(1:completed_steps, :), 1, 1);
        for m = 1:completed_steps-1
            for n = 1:4
                phi_val = phidot_raw(m, n);
                phi_val = mod(phi_val + pi, 2*pi) - pi;
                phidot_history(m, n) = phi_val / dt;
            end
        end
    end

    %% 计算 RMSE
    if completed_steps > 0
        sum1 = 0;
        for i = 1:completed_steps
            sum1 = sum1 + norm(q_history(i, :) - path(1:2, i))^2;
        end
        position_rmse = sqrt(sum1 / completed_steps);
    else
        position_rmse = NaN;
    end

    %% 全局 solver_time_array 汇总
    global solver_time_array
    if exist('solver_time_array', 'var') && ~isempty(solver_time_array)
        total_solve_time = sum(solver_time_array);
    else
        total_solve_time = NaN;
    end

    %% 组装 result
    result = struct();
    result.reference             = path;
    result.state_history         = q_history;
    result.orientation_history   = psi_history;
    result.body_velocity_history = velocity_history;
    result.world_velocity_history = state_dot_history;
    result.control_history       = u_history;
    result.wheel_speed_history   = vi_history;
    result.wheel_angle_history   = phi_history;
    result.steering_rate_history = phidot_history;
    result.solver_time_history   = solver_time_history;
    result.cvx_status_history    = cvx_status_history;
    result.cvx_optval_history    = cvx_optval_history;
    result.objective_history     = zeros(num_steps, 1);  % 累计 J 每步
    % 计算累计 J
    J_cum = 0;
    for i = 1:num_steps
        if i <= completed_steps
            e = [q_history(i,:), psi_history(i)]' - path(:, i);
            u_i = u_history(i, :)';
            J_cum = J_cum + e' * [30,0,0; 0,30,0; 0,0,1] * e + u_i' * [0.3,0,0; 0,0.3,0; 0,0,0.3] * u_i;
        end
        result.objective_history(i) = J_cum;
    end
    result.time_history          = t_history;
    result.trajectory_cost       = J;
    result.position_rmse         = position_rmse;
    result.max_wheel_speed       = max(vi_history(1:completed_steps, :), [], 'all');
    if completed_steps > 1
        result.max_steering_rate = max(abs(phidot_history(1:completed_steps-1, :)), [], 'all');
    else
        result.max_steering_rate = 0;
    end
    result.total_solve_time      = total_solve_time;
    result.completed_steps       = completed_steps;
    result.success               = success;
    result.failure_reason        = failure_reason;
    result.config                = cfg;
    result.wall_time             = toc(wall_tic);

    %% 可选保存图片
    if cfg.save_figures && cfg.live_plot
        if ~isempty(cfg.output_dir) && exist(cfg.output_dir, 'dir')
            saveas(fig, fullfile(cfg.output_dir, 'trajectory_plot.png'));
        end
    end

    %% 打印摘要
    fprintf('[run_closed_loop] Steps: %d/%d | RMSE: %.10f | J: %.10f | SolveTime: %.4fs | Wall: %.2fs\n', ...
        completed_steps, num_steps, position_rmse, J, total_solve_time, result.wall_time);
end
