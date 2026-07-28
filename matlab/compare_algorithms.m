function comparison = compare_algorithms(seed_range, algorithms, force_regen)
%COMPARE_ALGORITHMS  [兼容性包装] 批量对比主入口 (已迁移至 main.m)
%
% 此函数保留是为了兼容旧调用方式。新代码请直接使用:
%   comparison = main('seeds', 1:100, 'algorithms', {...}, 'forceRegen', false)
%
% 用法 (与原版一致):
%   comparison = compare_algorithms(1:10, {'e-lmpc', 'active-set', 'interior-point', 'proposed-3iter'})
%   comparison = compare_algorithms(1:100, {'proposed-3iter'})
%   comparison = compare_algorithms(1:200, {...}, true)   % 强制重新生成场景
%
% 实际逻辑已拆分到以下模块 (架构清晰):
%   - main.m                      主入口, 按步骤编排
%   - comparison_init.m           初始化 comparison 结构体
%   - comparison_load.m           加载 / 迁移 checkpoint
%   - comparison_save.m           增量保存 checkpoint
%   - run_batch_simulation.m      核心 (seed × algorithm) 仿真循环
%   - pair_is_completed.m         检查 (seed, alg) 是否完成
%   - find_pair_index.m           查找 (seed, alg) 行号
%   - get_alg_results.m           提取某算法的所有 results
%   - safe_metric.m               安全读取 metrics 字段
%   - plot_one_algorithm.m        单算法汇总图
%   - save_algorithm_csv.m        单算法 CSV 导出
%   - print_table_ii.m            Table II 汇总
%   - plot_paper_comparison.m     论文风格对比图 (已存在)

    if nargin < 1 || isempty(seed_range), seed_range = 1:10; end
    if nargin < 2 || isempty(algorithms)
        algorithms = {'e-lmpc', 'active-set', 'interior-point', 'proposed-3iter'};
    end
    if nargin < 3 || isempty(force_regen), force_regen = false; end

    comparison = main( ...
        'seeds', seed_range, ...
        'algorithms', algorithms, ...
        'forceRegen', force_regen);
end
