function save_algorithm_csv(comparison, alg_idx, results_dir)
%SAVE_ALGORITHM_CSV  将单个算法的所有 result 导出为 CSV
%
% 输入:
%   comparison : comparison 结构体
%   alg_idx    : 算法索引
%   results_dir: 结果根目录
%
% 输出:
%   {results_dir}/{算法名}_results.csv
%   (若原 CSV 被占用, 写到 {算法名}_results_new.csv)

    algorithm = comparison.algorithms{alg_idx};
    alg_results = get_alg_results(comparison, alg_idx);
    n = length(alg_results);

    col_names = { ...
        'seed', 'algorithm', 'success', ...
        'rmse', 'trajectoryCost', ...
        'meanSolveTime', 'maxSolveTime', 'totalSolveTime', ...
        'wheelSpeedViolationRatio', 'maxWheelSpeedViolation', ...
        'steeringRateViolationRatio', 'maxSteeringRateViolation', ...
        'successRate', 'elapsedTime' ...
    };

    T = table();
    for c = 1:length(col_names)
        if strcmp(col_names{c}, 'algorithm')
            T.(col_names{c}) = cell(n, 1);
        else
            T.(col_names{c}) = zeros(n, 1);
        end
    end

    for i = 1:n
        r = alg_results(i);
        m = r.metrics;

        T.seed(i) = r.seed;
        T.algorithm{i} = algorithm;
        T.success(i) = r.success;

        fnames = {'rmse', 'trajectoryCost', ...
                  'meanSolveTime', 'maxSolveTime', 'totalSolveTime', ...
                  'wheelSpeedViolationRatio', 'maxWheelSpeedViolation', ...
                  'steeringRateViolationRatio', 'maxSteeringRateViolation', ...
                  'successRate'};
        for f = 1:length(fnames)
            if isfield(m, fnames{f}) && ~isempty(m.(fnames{f})) && isnumeric(m.(fnames{f}))
                T.(fnames{f})(i) = m.(fnames{f});
            else
                T.(fnames{f})(i) = NaN;
            end
        end

        T.elapsedTime(i) = r.elapsedTime;
    end

    csv_file = fullfile(results_dir, sprintf('%s_results.csv', algorithm));

    % 尝试写入 CSV, 如果文件被占用则写到临时文件
    try
        writetable(T, csv_file);
        fprintf('[算法 %s] CSV 已保存: %s\n', algorithm, csv_file);
    catch
        % 文件被 Excel 等占用, 写到临时文件
        csv_tmp = fullfile(results_dir, sprintf('%s_results_new.csv', algorithm));
        writetable(T, csv_tmp);
        fprintf('[算法 %s] 原CSV被占用, 已保存到: %s\n', algorithm, csv_tmp);
        fprintf('[提示] 关闭Excel后, 可手动将 _new.csv 重命名覆盖原文件\n');
    end
end
