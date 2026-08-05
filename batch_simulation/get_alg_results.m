function alg_results = get_alg_results(comparison, alg_idx_or_name)
%GET_ALG_RESULTS  从 comparison 中提取某算法的所有 results 为 struct 数组
%
% 输入:
%   comparison       : compare_algorithms / main 输出的结构体
%   alg_idx_or_name  : 算法索引 (int) 或算法名 (string)
%
% 输出:
%   alg_results : struct 数组; 无匹配时为空 struct

    if ~isfield(comparison, 'completedPairs') || ~isfield(comparison, 'results')
        alg_results = struct([]);
        return;
    end

    if ischar(alg_idx_or_name) || isstring(alg_idx_or_name)
        algorithm = char(alg_idx_or_name);
    else
        algorithm = comparison.algorithms{alg_idx_or_name};
    end

    alg_mask = strcmp(comparison.completedPairs(:, 2), algorithm);
    n_match = sum(alg_mask);

    if n_match == 0
        alg_results = struct([]);
    elseif n_match == 1
        alg_results = comparison.results{find(alg_mask, 1)};
    else
        alg_results = [comparison.results{alg_mask}];
    end
end
