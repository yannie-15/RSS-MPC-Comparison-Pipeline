function comparison = run_batch_simulation(comparison, seeds, algorithms, force_regen, checkpoint_file)
%RUN_BATCH_SIMULATION  核心 (seed × algorithm) 批量仿真循环
%
% 不修改 checkpoint 逻辑, 仅是把原 compare_algorithms.m 中第 3 段 "逐 seed × 逐算法运行"
% 抽出为独立函数, 便于 main.m 编排调用。
%
% 设计要点 (与原行为完全一致):
%   - 外层遍历算法, 内层遍历 seed
%   - 切换算法时 clear functions 重置 e-LMPC warm start
%   - 已成功完成的 (seed, alg) 跳过; 失败的允许重跑
%   - 每完成 1 个组合 → 立即增量保存 checkpoint
%   - 某算法全部 seed 完成后, 即时生成该算法的汇总图与 CSV
%
% 输入:
%   comparison      : 已加载/初始化的 comparison 结构体 (会被原地更新)
%   seeds           : 种子行向量
%   algorithms      : 算法名 cell 数组
%   force_regen     : 是否强制重新生成场景
%   checkpoint_file : checkpoint 文件路径 (用于增量保存)
%
% 输出:
%   comparison : 更新后的结构体 (含本次新增的 results / completedPairs)

    num_seeds = length(seeds);
    num_algorithms = length(algorithms);

    total_tic = tic;

    for alg_idx = 1:num_algorithms
        algorithm = algorithms{alg_idx};

        % 清除 e-LMPC 的 warm start (切换算法时需要重置)
        clear functions;

        for seed_idx = 1:num_seeds
            seed_val = seeds(seed_idx);

            % 检查此 (seed, algorithm) 组合是否已成功完成
            % 只跳过成功的; 失败的重新跑
            if pair_is_completed(comparison.completedPairs, comparison.results, seed_val, algorithm)
                continue;
            end

            case_tic = tic;

            % 加载/生成此 seed 的场景
            [cfg, scen] = scenario_bank(seed_val, force_regen);

            % 设置 config 的算法标签
            cfg.algorithm = algorithm;

            % 关闭实时画图
            cfg.output.livePlot = false;
            cfg.output.saveFigures = false;

            fprintf('\n--- [%s] seed=%d (%d/%d seeds) ---\n', algorithm, seed_val, seed_idx, num_seeds);

            try
                summary = run_one_case(cfg, scen);
                summary.elapsedTime = toc(case_tic);
            catch ME
                summary = struct();
                summary.success = false;
                summary.failureReason = getReport(ME);
                summary.scenario = scen;
                summary.config = cfg;
                summary.algorithm = algorithm;
                summary.seed = seed_val;
                summary.elapsedTime = toc(case_tic);

                if isfield(summary, 'metrics')
                    summary.metrics.rmse = NaN;
                    summary.metrics.trajectoryCost = NaN;
                    summary.metrics.successRate = 0;
                else
                    summary.metrics = struct('rmse', NaN, 'trajectoryCost', NaN, 'successRate', 0);
                end

                fprintf('[seed=%d: %s] EXCEPTION: %s\n', seed_val, algorithm, ME.message);
                for si = 1:length(ME.stack)
                    fprintf('  at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
                end
            end

            % 记录 seed
            summary.seed = seed_val;

            % ========== 添加到 results 和 completedPairs ==========
            % 如果是重跑失败的组合，替换旧记录; 否则追加新记录
            existing_idx = find_pair_index(comparison.completedPairs, seed_val, algorithm);
            if ~isempty(existing_idx)
                comparison.results{existing_idx} = summary;
                % completedPairs 不需要更新 (seed_val, algorithm 已存在)
            else
                comparison.results{end+1} = summary;
                comparison.completedPairs(end+1, :) = {seed_val, algorithm};
            end

            % ========== 增量保存 ==========
            comparison_save(comparison, checkpoint_file);

            n_done = size(comparison.completedPairs, 1);
            fprintf('  → 已保存进度 (完成 %d/%d 组合)\n', n_done, num_seeds * num_algorithms);
        end

        % ========== 每种算法所有 seed 完成后，生成该算法的汇总 ==========
        alg_names = comparison.completedPairs(:, 2);
        n_alg_done = sum(strcmp(alg_names, algorithm));
        if n_alg_done >= num_seeds
            results_dir = fileparts(checkpoint_file);
            fprintf('[算法 %s] 全部 seed 完成，生成图像...\n', algorithm);
            try
                plot_one_algorithm(comparison, alg_idx, results_dir);
            catch ME_plot
                fprintf('[算法 %s] 绘图失败(跳过): %s\n', algorithm, ME_plot.message);
            end
            save_algorithm_csv(comparison, alg_idx, results_dir);
        end
    end

    comparison.totalElapsed = toc(total_tic);
end
