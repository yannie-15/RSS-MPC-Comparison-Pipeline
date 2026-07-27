% SETUP_PATHS 设置项目所需的所有 MATLAB 路径
% 在运行任何仿真脚本前执行此文件
%
% 用法：
%   cd d:\Projects\RSS\matlab
%   setup_paths

function setup_paths()
    script_dir = fileparts(mfilename('fullpath'));
    workspace_root = fileparts(script_dir);

    % 添加核心仿真目录
    addpath(fullfile(workspace_root, 'matlab'));
    addpath(fullfile(workspace_root, 'RSS_proposed'));

    fprintf('MATLAB 路径已设置：\n');
    fprintf('  - matlab/\n');
    fprintf('  - RSS_proposed/\n');
end
