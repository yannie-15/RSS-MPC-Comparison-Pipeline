function idx = find_pair_index(completedPairs, seed_val, algorithm)
%FIND_PAIR_INDEX  在 completedPairs 中查找 (seed_val, algorithm) 的行号
%
% 输入:
%   completedPairs : cell 数组, 每行 = {seed_val, algorithm_name}
%   seed_val       : 种子值
%   algorithm      : 算法名
%
% 输出:
%   idx : 行号; 未找到返回 []

    if isempty(completedPairs)
        idx = [];
        return;
    end
    seed_match = cell2mat(completedPairs(:, 1)) == seed_val;
    alg_match = strcmp(completedPairs(:, 2), algorithm);
    idx = find(seed_match & alg_match, 1);
end
