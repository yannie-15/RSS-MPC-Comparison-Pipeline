function completed = pair_is_completed(completedPairs, results, seed_val, algorithm)
%PAIR_IS_COMPLETED  判断 (seed, algorithm) 组合是否已成功完成
%
% 仅当该组合存在且对应 result.success=true 时返回 true。
% 失败的组合允许重跑, 因此返回 false。
%
% 输入:
%   completedPairs : cell 数组 {seed_val, algorithm_name}
%   results        : cell 数组, 与 completedPairs 行对齐
%   seed_val       : 种子值
%   algorithm      : 算法名
%
% 输出:
%   completed : true 表示该组合已成功, 应跳过; false 表示需要 (重) 跑

    if isempty(completedPairs)
        completed = false;
        return;
    end
    seed_match = cell2mat(completedPairs(:, 1)) == seed_val;
    alg_match = strcmp(completedPairs(:, 2), algorithm);
    pair_match = seed_match & alg_match;

    if ~any(pair_match)
        completed = false;
        return;
    end

    % 检查该组合的结果是否成功
    idx = find(pair_match, 1);
    if idx <= length(results) && isfield(results{idx}, 'success') && results{idx}.success
        completed = true;
    else
        completed = false;  % 失败或结果缺失 → 允许重跑
    end
end
