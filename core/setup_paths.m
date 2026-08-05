% SETUP_PATHS 设置项目所需的所有 MATLAB 路径
% 在运行任何仿真脚本前执行此文件
%
% 注: 四个算法目录 (algorithms/RSS_proposed, RSS_sqp, RSS_fmincon, RSS_active_set)
%     不在此处统一 addpath, 因为它们的 control_RSS.m 和 config.m 同名会冲突。
%     run_one_case.m / run_paper_baseline_case.m 会在调用时临时 addpath/rmpath 切换。
%     RSS_proposed 已转为普通目录 (非 submodule), 其余三个仍为 submodule。
%
% 用法：
%   cd d:\PROJECT\RSS_V2\core
%   setup_paths
% 或从任意功能目录调用前先 addpath('core')。

function setup_paths()
    script_dir = fileparts(mfilename('fullpath'));
    workspace_root = fileparts(script_dir);

    % 添加核心工具与批量仿真目录
    addpath(fullfile(workspace_root, 'core'));
    addpath(fullfile(workspace_root, 'batch_simulation'));
    addpath(fullfile(workspace_root, 'algorithms'));

    % 检查算法目录是否就绪 (RSS_proposed 为普通目录, 其余为 submodule)
    submodule_names = {'RSS_proposed', 'RSS_sqp', 'RSS_fmincon', 'RSS_active_set'};
    missing = {};
    for i = 1:length(submodule_names)
        sub_path = fullfile(workspace_root, 'algorithms', submodule_names{i});
        if ~exist(fullfile(sub_path, 'control_RSS.m'), 'file')
            missing{end+1} = submodule_names{i};
        end
    end
    if ~isempty(missing)
        warning('setup_paths:MissingAlgorithm', ...
            ['以下算法目录未就绪: %s\n', ...
             'RSS_proposed 为普通目录 (已随主仓库 checkout), 其余为 submodule, 请运行: git submodule update --init --recursive'], ...
            strjoin(missing, ', '));
    end

    fprintf('MATLAB 路径已设置：\n');
    fprintf('  - core/\n');
    fprintf('  - batch_simulation/\n');
    fprintf('  - algorithms/ (RSS_proposed 为普通目录; RSS_sqp/RSS_fmincon/RSS_active_set 为 submodule)\n');
end
