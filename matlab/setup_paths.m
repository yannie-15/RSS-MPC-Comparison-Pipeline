% SETUP_PATHS 设置项目所需的所有 MATLAB 路径
% 在运行任何仿真脚本前执行此文件
%
% 注: 四个算法 submodule (algorithms/RSS_proposed, RSS_sqp, RSS_fmincon, RSS_active_set)
%     不在此处统一 addpath, 因为它们的 control_RSS.m 和 config.m 同名会冲突。
%     run_one_case.m / run_paper_baseline_case.m 会在调用时临时 addpath/rmpath 切换。
%
% 用法：
%   cd d:\PROJECT\RSS_V2\matlab
%   setup_paths

function setup_paths()
    script_dir = fileparts(mfilename('fullpath'));
    workspace_root = fileparts(script_dir);

    % 添加核心仿真目录
    addpath(fullfile(workspace_root, 'matlab'));
    addpath(fullfile(workspace_root, 'algorithms'));

    % 检查 submodule 是否已初始化
    submodule_names = {'RSS_proposed', 'RSS_sqp', 'RSS_fmincon', 'RSS_active_set'};
    missing = {};
    for i = 1:length(submodule_names)
        sub_path = fullfile(workspace_root, 'algorithms', submodule_names{i});
        if ~exist(fullfile(sub_path, 'control_RSS.m'), 'file')
            missing{end+1} = submodule_names{i};
        end
    end
    if ~isempty(missing)
        warning('setup_paths:MissingSubmodule', ...
            ['以下 submodule 未初始化: %s\n', ...
             '请运行: git submodule update --init --recursive'], ...
            strjoin(missing, ', '));
    end

    fprintf('MATLAB 路径已设置：\n');
    fprintf('  - matlab/\n');
    fprintf('  - algorithms/ (4 个 submodule, 由 run_one_case.m / run_paper_baseline_case.m 动态切换)\n');
end
