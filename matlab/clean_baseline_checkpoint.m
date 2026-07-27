% 清除三种基线算法的 checkpoint 记录, 保留 proposed-3iter
% 用法: cd('d:\Projects\RSS\matlab'); setup_paths; clean_baseline_checkpoint;

checkpoint_file = fullfile('..', 'results', 'comparison_checkpoint.mat');

if ~exist(checkpoint_file, 'file')
    fprintf('未找到 checkpoint 文件\n');
    return;
end

loaded = load(checkpoint_file, 'comparison');
comparison = loaded.comparison;

baseline_algs = {'e-lmpc', 'active-set', 'interior-point'};

% 找出需要保留的索引 (proposed-3iter)
keep_mask = true(size(comparison.completedPairs, 1), 1);
for i = 1:size(comparison.completedPairs, 1)
    alg_name = comparison.completedPairs{i, 2};
    if any(strcmp(alg_name, baseline_algs))
        keep_mask(i) = false;
    end
end

n_removed = sum(~keep_mask);
n_kept = sum(keep_mask);

fprintf('checkpoint 中原有 %d 条记录\n', length(keep_mask));
fprintf('将删除 %d 条基线算法记录\n', n_removed);
fprintf('将保留 %d 条 proposed 记录\n', n_kept);

% 过滤
comparison.completedPairs = comparison.completedPairs(keep_mask, :);
comparison.results = comparison.results(keep_mask);

% 保存
save(checkpoint_file, 'comparison', '-v7.3');
fprintf('已更新 checkpoint: %s\n', checkpoint_file);

% 同时处理 comparison_final.mat
final_file = fullfile('..', 'results', 'comparison_final.mat');
if exist(final_file, 'file')
    loaded2 = load(final_file, 'comparison');
    comp2 = loaded2.comparison;
    keep_mask2 = true(size(comp2.completedPairs, 1), 1);
    for i = 1:size(comp2.completedPairs, 1)
        alg_name = comp2.completedPairs{i, 2};
        if any(strcmp(alg_name, baseline_algs))
            keep_mask2(i) = false;
        end
    end
    comp2.completedPairs = comp2.completedPairs(keep_mask2, :);
    comp2.results = comp2.results(keep_mask2);
    save(final_file, 'comparison', '-v7.3');
    fprintf('已更新 final: %s (删除 %d 条)\n', final_file, sum(~keep_mask2));
end
