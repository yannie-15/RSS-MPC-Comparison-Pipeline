function [worldVelocity, bodyVelocity, solve_time, iter_num, u_col, log_text] = ...
        matlab_control_bridge(k, lastBodyVelocity, state)
%MATLAB_CONTROL_BRIDGE MATLAB Engine 每步求解桥 (e-lmpc / interior-point / active-set)
%
% 架构 (pipeline 统一框架): 轨迹生成 / 参数 / 原始动力学闭环 / 评估 全在 Python
% (trajectory_generator.py + params.py + simulator.py); MATLAB 只在每步调用本函数
% 求解控制。由 pipeline/matlab_algorithm.py 在常驻 MATLAB Engine 会话中逐步调用:
%   [wv, bv, st, it, u, txt] = matlab_control_bridge(k, lastBodyVelocity, state)
%
% 初始化数据从 base workspace 读取 (Python 端 eng.workspace 一次性写入):
%   PIPELINE_ALGORITHM : 'e-lmpc' | 'interior-point' | 'active-set'
%   PIPELINE_SEED      : 场景 seed (0=paper_fixed, >=1=scenario_bank)
%   PIPELINE_PATH      : 3xN 参考轨迹 (Python trajectory_generator 生成, 与 MATLAB 同源)
%
% 首次调用完成一次性初始化 (persistent, 会话内只执行一次):
%   - addpath core/ (defaultConfig) + batch_simulation/ (scenario_bank)
%   - 组装 config: seed=0 -> defaultConfig; seed>=1 -> scenario_bank(seed)
%   - addpath 算法目录 (algorithms/RSS_sqp | RSS_fmincon | RSS_active_set,
%     均为随主仓库 checkout 的普通目录, 原 submodule 已并入主仓库)
%   - e-lmpc / active-set: 写临时 config.m 覆盖算法目录 config()
%     (二者内部调 config() 取场景参数; interior-point 直接接收 config 第 5 参)
%
% 每步: 按算法分发 control_RSS (接口差异与 run_one_case 一致):
%   e-lmpc:          [wv, bv, st, it] = control_RSS(path, k, vel, state)
%   interior-point:  [wv, bv, st, it] = control_RSS(path, k, vel, state, config)
%   active-set:      [wv, bv, st, it] = control_RSS(path, k, vel, state, config)
% evalc 捕获算法的 exitflag/status 打印, 回传 Python 统一输出。

    persistent ready algorithm pass_config config path_ref
    if isempty(ready)
        this_dir = fileparts(mfilename('fullpath'));   % <repo>/pipeline
        repo_root = fileparts(this_dir);
        addpath(fullfile(repo_root, 'core'));          % defaultConfig
        addpath(fullfile(repo_root, 'batch_simulation'));  % scenario_bank

        algorithm = lower(evalin('base', 'PIPELINE_ALGORITHM'));
        seed = double(evalin('base', 'PIPELINE_SEED'));
        path_ref = evalin('base', 'PIPELINE_PATH');

        alg_map = containers.Map( ...
            {'e-lmpc', 'interior-point', 'active-set'}, ...
            {'RSS_sqp', 'RSS_fmincon', 'RSS_active_set'});
        if ~isKey(alg_map, algorithm)
            error('matlab_control_bridge:BadAlgorithm', ...
                '本桥仅服务 MATLAB 算法 (e-lmpc/interior-point/active-set), 收到: %s', algorithm);
        end
        submodule_dir = fullfile(repo_root, 'algorithms', alg_map(algorithm));
        if ~exist(fullfile(submodule_dir, 'control_RSS.m'), 'file')
            error('matlab_control_bridge:MissingAlgorithm', ...
                '缺少算法目录: %s (普通目录, 请检查仓库完整性)', submodule_dir);
        end

        % 组装算法 config (与 others/batch_simulation/run_one_case.m 同源)
        if seed == 0
            config = defaultConfig();
        else
            [config, ~] = scenario_bank(seed);
        end
        config.algorithm = algorithm;

        addpath(submodule_dir);
        pass_config = ~strcmp(algorithm, 'e-lmpc');
        if ~pass_config || strcmp(algorithm, 'active-set')
            % e-lmpc / active-set 的控制器内部调 config(): 写临时 config.m 覆盖
            override_dir = setup_config_override(config);
            addpath(override_dir);
        end
        ready = true;
    end

    if pass_config
        call_expr = 'control_RSS(path_ref, k, lastBodyVelocity, state, config)';
    else
        call_expr = 'control_RSS(path_ref, k, lastBodyVelocity, state)';
    end
    % evalc: 输出变量落在当前工作区; 捕获 exitflag 打印回传 Python
    log_text = evalc(sprintf( ...
        '[worldVelocity, bodyVelocity, solve_time, iter_num] = %s;', call_expr));

    % u(:,1) = bodyVelocity - lastBodyVelocity (与 run_one_case 反推口径一致)
    u_col = bodyVelocity - lastBodyVelocity;
end


%% =========================================================
% 辅助函数: 创建临时 config.m 覆盖算法目录的 config
% (与 others/batch_simulation/run_one_case.m 的同名辅助函数保持一致)
%% ==========================================================

function override_dir = setup_config_override(cfg)
%SETUP_CONFIG_OVERRIDE 创建临时 config.m, 用场景 config 覆盖算法目录的 config
%
% e-lmpc / active-set 的控制器内部调 config() 取参数, 不接收外部 config。
% 本函数在临时目录写一个 config.m, 返回场景 cfg, 通过 addpath 覆盖。

    override_dir = fullfile(tempdir, 'rss_config_override');
    if ~exist(override_dir, 'dir'), mkdir(override_dir); end

    % 保存 cfg 到 mat (每次调用更新, 因为不同 seed 的 config 不同)
    cfg_file = fullfile(override_dir, 'config_data.mat');
    save(cfg_file, 'cfg', '-v7');

    % 写 config.m (只写一次, 复用)
    config_m = fullfile(override_dir, 'config.m');
    if ~exist(config_m, 'file')
        fid = fopen(config_m, 'w');
        fprintf(fid, 'function params = config()\n');
        fprintf(fid, '    f = fullfile(fileparts(mfilename(''fullpath'')), ''config_data.mat'');\n');
        fprintf(fid, '    loaded = load(f, ''cfg'');\n');
        fprintf(fid, '    params = loaded.cfg;\n');
        fprintf(fid, 'end\n');
        fclose(fid);
    end
end
