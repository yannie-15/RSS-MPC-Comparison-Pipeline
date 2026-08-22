function comparison = comparison_load(checkpoint_file, seeds, algorithms, force_regen)
%COMPARISON_LOAD  加载已有 checkpoint, 兼容旧格式; 不存在则初始化新结构
%
% 与原 compare_algorithms.m 的 "检查断点续跑" 段落保持完全等价的逻辑,
% 仅做名字上的拆分以便 main.m 调用。
%
% 输入:
%   checkpoint_file : checkpoint 文件全路径
%   seeds           : 当前期望的种子行向量 (用于扩大范围时更新元数据)
%   algorithms      : 当前期望的算法 cell 数组
%   force_regen     : 是否强制重新生成 (true → 即使有 checkpoint 也初始化新的)
%
% 输出:
%   comparison : 加载 (并迁移) 后的 comparison 结构体; 无文件时为新初始化

    if exist(checkpoint_file, 'file') && ~force_regen
        fprintf('发现已有结果文件: %s\n', checkpoint_file);
        loaded = load(checkpoint_file, 'comparison');
        comparison = loaded.comparison;

        % 兼容旧格式: 旧版 completedPairs 是数值矩阵 [seed, alg_idx]
        % 需要转换为 cell 数组 {seed, algorithm_name}
        if ~iscell(comparison.completedPairs)
            old_pairs = comparison.completedPairs;
            new_pairs = cell(size(old_pairs, 1), 2);
            for r = 1:size(old_pairs, 1)
                new_pairs(r, 1) = {old_pairs(r, 1)};
                old_idx = old_pairs(r, 2);
                if old_idx <= length(comparison.algorithms)
                    new_pairs(r, 2) = {comparison.algorithms{old_idx}};
                else
                    new_pairs(r, 2) = {'unknown'};
                end
            end
            comparison.completedPairs = new_pairs;
        end

        % 更新元数据以适配扩大范围 (保留已完成数据)
        comparison.seeds = seeds;
        comparison.algorithms = algorithms;
        comparison.numSeeds = length(seeds);
        comparison.numAlgorithms = length(algorithms);

        % 统计已完成数
        n_completed = size(comparison.completedPairs, 1);
        total_pairs = length(seeds) * length(algorithms);

        fprintf('断点续跑: 已完成 %d/%d 组合\n', n_completed, total_pairs);
    else
        comparison = comparison_init(seeds, algorithms);
    end
end
