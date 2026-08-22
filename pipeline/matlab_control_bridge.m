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
%   - 优先从 Python 注入的 PIPELINE_CONFIG (base workspace) 组装 config
%   - 缺失时兜底: 旧架构 defaultConfig/scenario_bank (已归档至 others/,
%     探测存在才 addpath; 找不到时报错说明应经 matlab_algorithm.py 启动)
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
        % 旧架构 core/ (defaultConfig) 与 batch_simulation/ (scenario_bank)
        % 已归档至 others/; 仅在目录存在时加入 (PIPELINE_CONFIG 缺失时的
        % 兜底需要, 正常流程不需要)
        cand = { ...
            fullfile(repo_root, 'core'), ...
            fullfile(repo_root, 'batch_simulation'), ...
            fullfile(repo_root, 'others', 'batch_simulation', 'batch_simulation')};
        dd = dir(fullfile(repo_root, 'others', 'batch_results_*', 'core'));
        for ii = 1:numel(dd)
            cand{end+1} = dd(ii).folder;   %#ok<SAGROW>
        end
        for ii = 1:numel(cand)
            if exist(cand{ii}, 'dir') == 7
                addpath(cand{ii});
            end
        end

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
        % 优先使用 Python 注入的完整 PIPELINE_CONFIG（若存在）以保证参数一致性
        if evalin('base', 'exist(''PIPELINE_CONFIG'', ''var'')') == 1
            try
                pc = evalin('base', 'PIPELINE_CONFIG');
                % 重建为算法期望的 config 结构体（兼容原 MATLAB config.m 字段）
                config = struct();
                % 车辆几何与约束
                if isfield(pc, 'algorithm_params') && isfield(pc.algorithm_params, 'vehicle')
                    v = pc.algorithm_params.vehicle;
                    if isfield(v, 'Lx'), config.Lx = v.Lx; else config.Lx = 0.655; end
                    if isfield(v, 'Ly'), config.Ly = v.Ly; else config.Ly = 0.335; end
                    if isfield(v, 'wheel_pos'), config.wheel_pos = v.wheel_pos; else config.wheel_pos = [config.Lx/2, config.Ly/2; -config.Lx/2, config.Ly/2; -config.Lx/2, -config.Ly/2; config.Lx/2, -config.Ly/2]; end
                    if isfield(v, 'vimax'), config.vimax = v.vimax; else config.vimax = 5; end
                    if isfield(v, 'phidotmax'), config.phidotmax = v.phidotmax; else config.phidotmax = 5 * pi; end
                else
                    % 无 vehicle 信息，回退到旧架构配置
                    config = legacy_config(seed);
                end

                % 控制器与正则化参数（尝试从顶层或 algorithm_params.weights 映射, 否则默认）
                if isfield(pc, 'k1'), config.k1 = pc.k1; end
                if ~isfield(config, 'k1') && isfield(pc, 'algorithm_params') && isfield(pc.algorithm_params, 'weights') && isfield(pc.algorithm_params.weights, 'w_pos')
                    config.k1 = pc.algorithm_params.weights.w_pos;
                end
                if ~isfield(config, 'k1'), config.k1 = 0.15; end

                if isfield(pc, 'k2'), config.k2 = pc.k2; end
                if ~isfield(config, 'k2') && isfield(pc, 'algorithm_params') && isfield(pc.algorithm_params, 'weights') && isfield(pc.algorithm_params.weights, 'w_psi')
                    config.k2 = pc.algorithm_params.weights.w_psi;
                end
                if ~isfield(config, 'k2'), config.k2 = 0.15; end

                if isfield(pc, 'k3'), config.k3 = pc.k3; end
                if ~isfield(config, 'k3') && isfield(pc, 'algorithm_params') && isfield(pc.algorithm_params, 'weights') && isfield(pc.algorithm_params.weights, 'w_control')
                    config.k3 = pc.algorithm_params.weights.w_control;
                end
                if ~isfield(config, 'k3'), config.k3 = 0.1; end

                if isfield(pc, 'eps'), config.eps = pc.eps; end
                if ~isfield(config, 'eps') && isfield(pc, 'algorithm_params') && isfield(pc.algorithm_params, 'weights') && isfield(pc.algorithm_params.weights, 'rho')
                    config.eps = pc.algorithm_params.weights.rho;
                end
                if ~isfield(config, 'eps'), config.eps = 0.001; end

                % 时间/步数/路径点
                if isfield(pc, 'dt'), config.dt = pc.dt; elseif isfield(pc, 'algorithm_params') && isfield(pc.algorithm_params, 'dt'), config.dt = pc.algorithm_params.dt; else config.dt = 0.01; end
                if isfield(pc, 'num_steps'), config.num_steps = pc.num_steps; elseif isfield(pc, 'num_steps'), config.num_steps = pc.num_steps; else config.num_steps = round(1.0 / config.dt); end
                if isfield(pc, 'num_path_pts'), config.num_path_pts = pc.num_path_pts; elseif isfield(pc, 'algorithm_params') && isfield(pc.algorithm_params, 'num_path_pts'), config.num_path_pts = pc.algorithm_params.num_path_pts; else config.num_path_pts = config.num_steps; end
                if isfield(pc, 'ctrl_pts'), config.ctrl_pts = pc.ctrl_pts; elseif isfield(pc, 'algorithm_params') && isfield(pc.algorithm_params, 'ctrl_pts'), config.ctrl_pts = pc.algorithm_params.ctrl_pts; end

                % 预测时域 K
                if isfield(pc, 'K'), config.K = pc.K; elseif isfield(pc, 'algorithm_params') && isfield(pc.algorithm_params, 'K'), config.K = pc.algorithm_params.K; end

            catch
                % 解析失败，回退到旧架构配置
                config = legacy_config(seed);
            end
        else
            % Python 未注入 PIPELINE_CONFIG, 回退到旧架构配置
            config = legacy_config(seed);
        end
        % 如果 Python 端传入 PIPELINE_K，则覆盖 config 中的 K 字段
        if evalin('base', 'exist(''PIPELINE_K'', ''var'')') == 1
            try
                config.K = double(evalin('base', 'PIPELINE_K'));
            catch
                % 忽略转换错误，保留原 config
            end
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

    % ---- 每步 rho 覆盖 (Python --rho 逐步序列经 PIPELINE_RHO 注入;
    %      未注入时保持 config 默认, 不改变基线行为) ----
    if evalin('base', 'exist(''PIPELINE_RHO'', ''var'')') == 1
        try
            config.rho = double(evalin('base', 'PIPELINE_RHO'));
            if ~pass_config
                % e-lmpc 内部经 config() 读参数 (override 的 config.m 每次调用
                % 都重新加载 config_data.mat): 重写使其读到当前步 rho
                setup_config_override(config);
            end
        catch
            % 覆盖失败则保持默认
        end
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
% 辅助函数: PIPELINE_CONFIG 缺失/解析失败时的兜底配置
%% ==========================================================

function cfg = legacy_config(seed)
%LEGACY_CONFIG 旧架构 defaultConfig/scenario_bank 兜底 (正常流程不触发)
%
% core/ 与 batch_simulation/ 已归档至 others/ (初始化时探测 addpath);
% 仍未找到时给出明确错误 (正常流程应经 pipeline/matlab_algorithm.py
% 启动并注入 PIPELINE_CONFIG, 不依赖旧架构文件)。

    if exist('defaultConfig', 'file') ~= 2
        error('matlab_control_bridge:NoLegacyConfig', ...
            ['未找到 defaultConfig/scenario_bank (旧架构 core/ 与 batch_simulation/ 已归档至 others/).\n' ...
             '请经 pipeline/matlab_algorithm.py 启动 (会注入 PIPELINE_CONFIG),\n' ...
             '或手动 addpath others/ 下对应目录后再调用本桥.']);
    end
    if seed == 0
        cfg = defaultConfig();
    else
        [cfg, ~] = scenario_bank(seed);
    end
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
