function freeze_baseline()
% FREEZE_BASELINE  Phase 0: 运行原 algorithms/RSS_proposed/main.m 逻辑, 保存基线数据
%
% 此脚本忠实复现原 main.m 的完整 100 步闭环, 不做任何修改:
%   - 使用 algorithms/RSS_proposed/ 下的 config.m, bezier_path.m, control_RSS.m
%   - 轮角用 atan2(vyi, vxi) (原 main.m 公式, 非 computeWheelOutputs)
%   - 状态推进用 new_state_dot (世界系)
%   - 不开实时绘图
%
% 输出: tests/baseline/rss_cvx_original_main.mat
%       tests/baseline/baseline_meta.json

    %% 定位仓库根目录
    script_dir = fileparts(mfilename('fullpath'));
    repo_root = fileparts(fileparts(script_dir));  % tests/baseline -> tests -> repo_root
    submodule_dir = fullfile(repo_root, 'algorithms', 'RSS_proposed');

    if ~exist(fullfile(submodule_dir, 'control_RSS.m'), 'file')
        error('freeze_baseline:MissingSubmodule', ...
            'RSS_proposed submodule 未初始化: %s', submodule_dir);
    end

    %% addpath submodule (临时, 函数退出自动 rmpath)
    addpath(submodule_dir);
    cleanup_obj = onCleanup(@() rmpath(submodule_dir));
    clear functions  % 确保用到 submodule 版本的函数

    %% 加载配置 (使用 submodule 自己的 config.m)
    params = config();

    %% 生成参考轨迹
    path = bezier_path(params.ctrl_pts, params.num_path_pts);

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
    num_steps = params.num_steps;
    q_history     = zeros(num_steps, 2);
    vi_history    = zeros(num_steps, 4);
    phidot_history = zeros(num_steps, 4);
    psi_history   = zeros(num_steps, 1);
    t_history     = zeros(num_steps, 1);
    phi_history   = zeros(num_steps, 4);
    u_history     = zeros(num_steps, 3);       % 每步 u(:,1) 控制增量
    velocity_history = zeros(num_steps, 3);    % 每步车体系速度
    state_dot_history = zeros(num_steps, 3);   % 每步世界系状态导数
    cvx_optval_history = zeros(num_steps, 1);
    cvx_status_history = cell(num_steps, 1);

    %% 仿真循环 (忠实复现原 main.m)
    for k = 1:num_steps
        t = (k - 1) * params.dt;
        t_history(k) = t;
        q_history(k, :) = q;
        psi_history(k) = psi0;

        % 跟踪误差代价
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

        % 计算轮速和轮角 (原 main.m 公式: atan2(vyi, vxi))
        phi = zeros(4, 1);
        vi = zeros(4, 1);
        for i = 1:4
            Hj = [1, 0, -params.wheel_pos(i,2); 0, 1, params.wheel_pos(i,1)];
            zn = Hj * velocity;
            vxi = zn(1);
            vyi = zn(2);
            vi(i) = sqrt(vxi^2 + vyi^2);
            phi(i) = atan2(vyi, vxi);
        end

        vi_history(k, :) = vi';
        phi_history(k, :) = phi';

        % 状态更新
        vx = new_state_dot(1);
        vy = new_state_dot(2);
        omega_b = new_state_dot(3);
        psi0 = psi0 + omega_b * params.dt;
        q = q + [vx, vy] * params.dt;
        state_dot = new_state_dot;
        state = [q, psi0];

        prev_phi = phi;
    end

    %% 计算转向率 (原 main.m 方式)
    phidot_raw = diff(phi_history, 1, 1);  % (num_steps-1) x 4
    for m = 1:num_steps-1
        for n = 1:4
            phi_val = phidot_raw(m, n);
            phi_val = mod(phi_val + pi, 2*pi) - pi;
            phidot_history(m, n) = phi_val / params.dt;
        end
    end
    % 最后一行保持 0 (原 main.m 行为: diff 产出 num_steps-1 行, 存入 num_steps 数组)

    %% 计算 RMSE
    sum1 = 0;
    for i = 1:num_steps
        sum1 = sum1 + norm(q_history(i, :) - path(1:2, i))^2;
    end
    rmse = sqrt(sum1 / num_steps);

    %% 全局 solver_time_array
    global solver_time_array
    if exist('solver_time_array', 'var') && ~isempty(solver_time_array)
        total_solve_time = sum(solver_time_array);
    else
        total_solve_time = NaN;
        solver_time_array = [];
    end

    %% 组装 baseline 结构体
    baseline = struct();
    baseline.reference = path;
    baseline.q_history = q_history;
    baseline.psi_history = psi_history;
    baseline.t_history = t_history;
    baseline.vi_history = vi_history;
    baseline.phi_history = phi_history;
    baseline.phidot_history = phidot_history;
    baseline.u_history = u_history;
    baseline.velocity_history = velocity_history;
    baseline.state_dot_history = state_dot_history;
    baseline.cvx_optval_history = cvx_optval_history;
    baseline.cvx_status_history = cvx_status_history;
    baseline.solver_time_array = solver_time_array;
    baseline.trajectory_cost = J;
    baseline.rmse = rmse;
    baseline.total_solve_time = total_solve_time;
    baseline.params = params;

    %% 保存
    output_mat = fullfile(script_dir, 'rss_cvx_original_main.mat');
    save(output_mat, 'baseline', '-v7.3');
    fprintf('基线已保存: %s\n', output_mat);

    %% 保存元信息 JSON
    meta = struct();
    meta.git_commit = get_git_commit(repo_root);
    meta.matlab_version = version;
    meta.freeze_date = string(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    try
        [sel_solver, ~] = cvx_solver;
        meta.cvx_solver = sel_solver;
    catch
        meta.cvx_solver = 'unknown';
    end
    meta.num_steps = num_steps;
    meta.rmse = rmse;
    meta.trajectory_cost = J;
    meta.total_solve_time = total_solve_time;

    meta_json = fullfile(script_dir, 'baseline_meta.json');
    fid = fopen(meta_json, 'w');
    fprintf(fid, '%s', jsonencode(meta));
    fclose(fid);
    fprintf('元信息已保存: %s\n', meta_json);

    %% 打印摘要
    fprintf('\n===== 基线摘要 =====\n');
    fprintf('Git commit: %s\n', meta.git_commit);
    fprintf('MATLAB: %s\n', meta.matlab_version);
    fprintf('CVX solver: %s\n', meta.cvx_solver);
    fprintf('Steps: %d\n', num_steps);
    fprintf('RMSE: %.10f\n', rmse);
    fprintf('Trajectory cost J: %.10f\n', J);
    fprintf('Total solve time: %.6f s\n', total_solve_time);
    fprintf('Max wheel speed: %.6f m/s (limit: %.1f)\n', max(vi_history(:)), params.vimax);
    fprintf('Max steering rate: %.6f rad/s (limit: %.6f)\n', max(abs(phidot_history(:))), params.phidotmax);
end


function commit = get_git_commit(repo_root)
    try
        [~, commit] = system(['git -C "' repo_root '" rev-parse HEAD']);
        commit = strtrim(commit);
    catch
        commit = 'unknown';
    end
end
