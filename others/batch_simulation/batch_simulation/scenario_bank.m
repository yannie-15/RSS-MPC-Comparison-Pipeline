function [config, scenario] = scenario_bank(seed, force_regen)
% SCENARIO_BANK  每个 seed 生成 1 个可复现场景，保存为单独 .mat
%
% 用法:
%   [config, scenario] = scenario_bank(7)           % 加载/生成 seed=7 的场景
%   [config, scenario] = scenario_bank(7, true)     % 强制重新生成
%
% 设计:
%   - rng(seed) 固定种子 → 生成 1 个场景
%   - 保存到 scenario_bank/scenario_seed007.mat
%   - 已有文件直接加载（除非 force_regen=true）
%   - seed 之间互不干扰，增删不影响已有结果
%
% 输入:
%   seed        : 随机种子值（也是场景 ID）
%   force_regen : 是否强制重新生成（默认 false）
%
% 输出:
%   config    : 该场景的配置结构体
%   scenario  : 该场景的描述结构体

    if nargin < 1 || isempty(seed), seed = 1; end
    if nargin < 2 || isempty(force_regen), force_regen = false; end

    %% =====================================================
    % 确定保存路径
    % ======================================================

    script_dir = fileparts(mfilename('fullpath'));
    workspace_root = fileparts(script_dir);
    bank_dir = fullfile(workspace_root, 'scenario_bank');

    if ~exist(bank_dir, 'dir')
        mkdir(bank_dir);
    end

    mat_filename = fullfile(bank_dir, sprintf('scenario_seed%d.mat', seed));

    %% =====================================================
    % 如果已存在且不强制重生成，直接加载
    % ======================================================

    if exist(mat_filename, 'file') && ~force_regen
        fprintf('加载已有场景: seed=%d (%s)\n', seed, mat_filename);
        loaded = load(mat_filename, 'config', 'scenario');
        config = loaded.config;
        scenario = loaded.scenario;
        return;
    end

    %% =====================================================
    % 生成新场景: rng(seed) → scenario_generator(1, 'random', seed)
    % ======================================================

    fprintf('生成新场景: seed=%d\n', seed);

    % 用 rng(seed) 保证可复现
    rng(seed);

    % 调用 scenario_generator，每个 seed 生成 1 个场景
    [configs, scenarios] = scenario_generator(1, 'random', seed);

    config = configs;
    scenario = scenarios;

    % 标记 seed ID
    scenario.seed = seed;
    scenario.name = sprintf('seed_%d', seed);
    scenario.id = seed;

    %% =====================================================
    % 保存到磁盘
    % ======================================================

    save(mat_filename, 'config', 'scenario', '-v7.3');
    fprintf('场景已保存: %s\n', mat_filename);

end
