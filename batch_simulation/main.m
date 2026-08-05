function comparison = main(varargin)
%MAIN  批量仿真主入口 - 按步骤编排
%
% 这是整个项目的主入口文件。按照清晰的步骤顺序编排各模块调用,
% 每一步对应一个独立的功能模块, 便于阅读、维护和扩展。
%
% 用法:
%   main                                                  % 默认: 1:10 seeds, 全部 4 种算法
%   main('seeds', 1:100)                                  % 自定义 seed 范围
%   main('algorithms', {'proposed-3iter'})                % 只跑一种算法
%   main('seeds', 1:50, 'algorithms', {'proposed-3iter'}, 'forceRegen', false)
%   comparison = main(...)                                % 返回 comparison 结构体
%
% 步骤概览:
%   Step 0: 路径设置
%   Step 1: 生成 seed 列表
%   Step 2: 选择算法
%   Step 3: 检查断点续跑 (加载 / 初始化 comparison)
%   Step 4: 运行批量仿真 (含增量保存, 每算法完成后即时出图/CSV)
%   Step 5: 保存最终结果 (comparison_final.mat + 同步 checkpoint)
%   Step 6: 生成论文对比图
%   Step 7: 打印 Table II 汇总
%
% 输入 (name-value):
%   'seeds'       : 种子范围, 默认 1:10
%   'algorithms'  : 算法 cell 数组, 默认 4 种
%   'forceRegen'  : 是否强制重新生成场景文件, 默认 false

    %% ====== 参数解析 ======
    opts = parse_args(varargin{:});

    %% =====================================================
    %  Step 0: 路径设置
    %  ======================================================
    fprintf('\n========================================\n');
    fprintf('Step 0: 路径设置\n');
    fprintf('========================================\n');
    % 定位 core/ 目录 (setup_paths.m 所在), 兼容从 batch_simulation/ 或其他目录调用
    script_dir = fileparts(mfilename('fullpath'));
    repo_root = fileparts(script_dir);
    core_dir = fullfile(repo_root, 'core');
    if exist(core_dir, 'dir')
        addpath(core_dir);
    end
    setup_paths();

    %% =====================================================
    %  Step 1: 生成 seed 列表
    %  ======================================================
    fprintf('\n========================================\n');
    fprintf('Step 1: 生成 seed 列表\n');
    fprintf('========================================\n');
    seeds = opts.seeds(:)';   % 确保为行向量
    fprintf('种子范围: %d-%d (共 %d 个)\n', seeds(1), seeds(end), length(seeds));

    %% =====================================================
    %  Step 2: 选择算法
    %  ======================================================
    fprintf('\n========================================\n');
    fprintf('Step 2: 选择算法\n');
    fprintf('========================================\n');
    algorithms = opts.algorithms;
    fprintf('算法: %s\n', strjoin(algorithms, ', '));
    fprintf('总组合数: %d\n', length(seeds) * length(algorithms));

    %% =====================================================
    %  Step 3: 检查断点续跑 (加载 / 初始化 comparison)
    %  ======================================================
    fprintf('\n========================================\n');
    fprintf('Step 3: 检查断点续跑\n');
    fprintf('========================================\n');

    workspace_root = fileparts(fileparts(mfilename('fullpath')));
    results_dir = fullfile(workspace_root, 'results', 'batch');
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end

    checkpoint_file = fullfile(results_dir, 'comparison_checkpoint.mat');
    comparison = comparison_load(checkpoint_file, seeds, algorithms, opts.forceRegen);

    %% =====================================================
    %  Step 4: 运行批量仿真 (含增量保存)
    %  ======================================================
    fprintf('\n========================================\n');
    fprintf('Step 4: 运行批量仿真\n');
    fprintf('========================================\n');

    comparison = run_batch_simulation( ...
        comparison, seeds, algorithms, opts.forceRegen, checkpoint_file);

    %% =====================================================
    %  Step 5: 保存最终结果 (comparison_final.mat + 同步 checkpoint)
    %  注: 每算法的 summary 图与 CSV 已在 Step 4 中即时生成
    %  ======================================================
    fprintf('\n========================================\n');
    fprintf('Step 5: 保存最终结果\n');
    fprintf('========================================\n');
    save_all_artifacts(comparison, results_dir);

    %% =====================================================
    %  Step 6: 生成论文对比图
    %  ======================================================
    fprintf('\n========================================\n');
    fprintf('Step 6: 生成论文对比图\n');
    fprintf('========================================\n');
    plot_paper_comparison(comparison, results_dir);

    %% =====================================================
    %  Step 7: 打印 Table II 汇总
    %  ======================================================
    fprintf('\n========================================\n');
    fprintf('Step 7: 打印 Table II 汇总\n');
    fprintf('========================================\n');
    print_table_ii(comparison);

    %% ====== 收尾 ======
    fprintf('\n========================================\n');
    fprintf('全部步骤完成!\n');
    fprintf('  结果目录: %s\n', results_dir);
    fprintf('  Checkpoint: %s\n', checkpoint_file);
    fprintf('========================================\n');
end


%% =========================================================
%  参数解析
%% =========================================================

function opts = parse_args(varargin)
%PARSE_ARGS  解析 name-value 参数
%
% 支持的键:
%   'seeds'       : 种子范围 (默认 1:10)
%   'algorithms'  : 算法 cell 数组 (默认 4 种)
%   'forceRegen'  : 是否强制重新生成场景 (默认 false)

    opts = struct();
    opts.seeds = 1:10;
    opts.algorithms = {'e-lmpc', 'active-set', 'interior-point', 'proposed-3iter'};
    opts.forceRegen = false;

    % 解析 name-value 对
    i = 1;
    while i <= length(varargin)
        if ischar(varargin{i}) || isstring(varargin{i})
            key = lower(char(varargin{i}));
            if i + 1 <= length(varargin)
                value = varargin{i+1};
                switch key
                    case 'seeds'
                        opts.seeds = value;
                    case 'algorithms'
                        opts.algorithms = value;
                    case 'forceregen'
                        opts.forceRegen = logical(value);
                    otherwise
                        warning('main:unknownArg', '未知参数: %s', key);
                end
                i = i + 2;
            else
                error('main:missingValue', '参数 %s 缺少值', key);
            end
        else
            error('main:invalidArg', '参数名必须是字符串, 得到 %s', class(varargin{i}));
        end
    end
end


%% =========================================================
%  Step 5 的具体实现: 保存最终 comparison_final.mat + 同步 checkpoint
%% =========================================================

function save_all_artifacts(comparison, results_dir)
%SAVE_ALL_ARTIFACTS  保存最终 comparison_final.mat 并同步更新 checkpoint
%
% 注意: 每个算法的 summary 图与 CSV 已在 run_batch_simulation 内部,
%       该算法所有 seed 完成后即时生成。
%       逐 seed 轨迹图由 replot_per_seed.m 单独处理。
%       此函数只负责收尾: 写 final.mat + 更新 checkpoint (不修改 checkpoint 逻辑)。

    % 1. 保存最终完整结果
    final_file = fullfile(results_dir, 'comparison_final.mat');
    save(final_file, 'comparison', '-v7.3');
    fprintf('完整结果已保存到: %s\n', final_file);

    % 2. 同步更新 checkpoint 文件 (保留以供下次扩大范围续跑)
    checkpoint_file = fullfile(results_dir, 'comparison_checkpoint.mat');
    comparison_save(comparison, checkpoint_file);
    fprintf('Checkpoint 已更新 (下次扩大范围可续跑)\n');
end
