function comparison = comparison_init(seeds, algorithms)
%COMPARISON_INIT  初始化一个空的 comparison 结构体
%
% 输入:
%   seeds      : 种子值行向量 (如 1:10)
%   algorithms : 算法名 cell 数组 (如 {'proposed-3iter', 'e-lmpc'})
%
% 输出:
%   comparison : 结构体, 字段如下:
%     .seeds          : 种子行向量
%     .algorithms     : 算法 cell 数组
%     .numSeeds       : 种子数
%     .numAlgorithms  : 算法数
%     .results        : cell 数组 (空), 每个 cell 是一个 summary struct
%     .completedPairs : cell 数组 (0x2), 每行 = {seed, algorithm_name}
%     .totalElapsed   : 0

    comparison = struct();
    comparison.seeds = seeds;
    comparison.algorithms = algorithms;
    comparison.numSeeds = length(seeds);
    comparison.numAlgorithms = length(algorithms);

    % results: 所有已完成的结果（cell 数组，每个 cell 是一个 summary struct）
    comparison.results = {};

    % completedPairs: 已完成的 (seed_value, algorithm_name) 记录
    % cell 数组, 每行 = {seed_val, algorithm_name}
    comparison.completedPairs = cell(0, 2);

    comparison.totalElapsed = 0;
end
