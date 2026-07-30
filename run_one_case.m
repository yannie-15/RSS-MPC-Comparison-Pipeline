function summary_json = run_one_case(config_json_path)
% RUN_ONE_CASE  Python 面向的 JSON 接口 - 单次 RSS 闭环仿真
%
% 此函数位于仓库根目录, 与 matlab/run_one_case.m (批量对比管道) 不同。
% 本函数接收 JSON 文件路径, 调用 run_closed_loop, 返回 summary JSON 字符串。
%
% 用法 (MATLAB):
%   summary_json = run_one_case('default_python_config.json')
%
% 用法 (Python MATLAB Engine):
%   summary = eng.run_one_case('/path/to/config.json', nargout=1)
%
% 输入:
%   config_json_path : JSON 配置文件路径
%
% 输出:
%   summary_json : JSON 字符串, 包含仿真摘要
%
% 副作用:
%   - 创建 cfg.output_dir 目录 (如不存在)
%   - 保存 result.mat 和 summary.json 到 output_dir

    %% 读取 JSON 配置
    if nargin < 1 || isempty(config_json_path)
        error('run_one_case:NoConfig', '必须提供 config_json_path');
    end

    json_text = fileread(config_json_path);
    user_cfg = jsondecode(json_text);

    %% 合并默认配置
    cfg = merge_case_config(user_cfg);

    %% 创建输出目录
    if ~isempty(cfg.output_dir)
        if ~exist(cfg.output_dir, 'dir')
            mkdir(cfg.output_dir);
        end
    else
        cfg.output_dir = fullfile(pwd, 'results', cfg.case_id);
        if ~exist(cfg.output_dir, 'dir')
            mkdir(cfg.output_dir);
        end
    end

    %% 运行闭环仿真
    result = [];
    success = false;
    failure_reason = '';

    try
        result = run_closed_loop(cfg);
        success = result.success;
        if ~success
            failure_reason = result.failure_reason;
        end
    catch ME
        success = false;
        failure_reason = getReport(ME, 'extended');
        fprintf('[run_one_case] EXCEPTION: %s\n', ME.message);
    end

    %% 构建 summary
    summary = build_summary(result, cfg, success, failure_reason);

    %% 保存 result.mat (如存在)
    if ~isempty(result)
        result_file = fullfile(cfg.output_dir, 'result.mat');
        try
            save(result_file, 'result', 'cfg', '-v7.3');
            summary.result_file = result_file;
        catch ME
            fprintf('[run_one_case] 保存 result.mat 失败: %s\n', ME.message);
            summary.result_file = '';
        end
    else
        summary.result_file = '';
    end

    %% 保存 summary.json
    summary_file = fullfile(cfg.output_dir, 'summary.json');
    try
        fid = fopen(summary_file, 'w');
        fprintf(fid, '%s', jsonencode(summary));
        fclose(fid);
    catch ME
        fprintf('[run_one_case] 保存 summary.json 失败: %s\n', ME.message);
    end

    %% 返回 JSON 字符串
    summary_json = jsonencode(summary);
end


%% =========================================================
%  辅助函数: 合并 JSON 配置与默认 config()
%% ==========================================================
function cfg = merge_case_config(user_cfg)
%MERGE_CASE_CONFIG  将 JSON 用户配置合并到 config() 默认值

    % 定位 submodule 获取默认 config
    script_dir = fileparts(mfilename('fullpath'));
    if isempty(script_dir)
        script_dir = pwd;
    end
    submodule_dir = fullfile(script_dir, 'algorithms', 'RSS_proposed');

    if ~exist(fullfile(submodule_dir, 'config.m'), 'file')
        parent = fileparts(script_dir);
        submodule_dir = fullfile(parent, 'algorithms', 'RSS_proposed');
    end

    addpath(submodule_dir);
    cleanup = onCleanup(@() rmpath(submodule_dir));
    clear functions;

    cfg = config();

    % 从 JSON 覆盖标量参数 (仅覆盖已存在的字段)
    json_fields = fieldnames(user_cfg);
    for i = 1:length(json_fields)
        fname = json_fields{i};
        if isfield(cfg, fname)
            cfg.(fname) = user_cfg.(fname);
        end
    end

    % 处理特殊字段
    if isfield(user_cfg, 'case_id') && ~isempty(user_cfg.case_id)
        cfg.case_id = user_cfg.case_id;
    else
        cfg.case_id = 'default';
    end

    if isfield(user_cfg, 'seed')
        cfg.seed = user_cfg.seed;
        % 设置随机种子 (用于未来随机场景)
        rng(user_cfg.seed);
    else
        cfg.seed = 0;
    end

    if isfield(user_cfg, 'live_plot')
        cfg.live_plot = logical(user_cfg.live_plot);
    else
        cfg.live_plot = false;
    end

    if isfield(user_cfg, 'save_figures')
        cfg.save_figures = logical(user_cfg.save_figures);
    else
        cfg.save_figures = false;
    end

    if isfield(user_cfg, 'save_full_log')
        cfg.save_full_log = logical(user_cfg.save_full_log);
    else
        cfg.save_full_log = false;
    end

    if isfield(user_cfg, 'output_dir') && ~isempty(user_cfg.output_dir)
        cfg.output_dir = user_cfg.output_dir;
    else
        cfg.output_dir = '';
    end

    if isfield(user_cfg, 'solver') && ~isempty(user_cfg.solver)
        cfg.solver = user_cfg.solver;
    else
        cfg.solver = 'sdpt3';
    end

    % 如果 num_steps 被覆盖, 同步 num_path_pts
    if isfield(user_cfg, 'num_steps')
        cfg.num_path_pts = cfg.num_steps;
    end

    % 确保有 output_dir 字段
    if ~isfield(cfg, 'output_dir')
        cfg.output_dir = '';
    end
end


%% =========================================================
%  辅助函数: 构建 summary 结构体
%% ==========================================================
function summary = build_summary(result, cfg, success, failure_reason)
%BUILD_SUMMARY  从 result 构建可 JSON 序列化的 summary

    summary = struct();
    summary.case_id = cfg.case_id;
    summary.success = success;
    summary.failure_reason = failure_reason;

    if ~isempty(result) && success
        summary.completed_steps = result.completed_steps;
        summary.position_rmse = result.position_rmse;
        summary.orientation_rmse = 0;  % TODO: 从 result 计算
        summary.trajectory_cost = result.trajectory_cost;
        summary.max_wheel_speed = result.max_wheel_speed;
        summary.max_steering_rate = result.max_steering_rate;
        summary.max_constraint_violation = max(0, ...
            max(result.max_wheel_speed - cfg.vimax, 0));
        summary.mean_solver_time = mean(result.solver_time_history( ...
            result.solver_time_history ~= 0));
        summary.total_solver_time = result.total_solve_time;
        summary.wall_time = result.wall_time;

        % CVX 状态汇总
        statuses = result.cvx_status_history(1:result.completed_steps);
        unique_statuses = unique(statuses);
        summary.cvx_status_summary = '';
        for i = 1:length(unique_statuses)
            if i > 1
                summary.cvx_status_summary = [summary.cvx_status_summary, '; '];
            end
            cnt = sum(strcmp(statuses, unique_statuses{i}));
            summary.cvx_status_summary = sprintf('%s%s (%d)', ...
                summary.cvx_status_summary, unique_statuses{i}, cnt);
        end
    else
        summary.completed_steps = 0;
        summary.position_rmse = 0;
        summary.orientation_rmse = 0;
        summary.trajectory_cost = 0;
        summary.max_wheel_speed = 0;
        summary.max_steering_rate = 0;
        summary.max_constraint_violation = 0;
        summary.mean_solver_time = 0;
        summary.total_solver_time = 0;
        summary.wall_time = 0;
        summary.cvx_status_summary = 'N/A';
    end

    summary.output_dir = cfg.output_dir;
    summary.num_steps = cfg.num_steps;
    summary.solver = cfg.solver;
end
