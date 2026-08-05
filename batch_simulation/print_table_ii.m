function print_table_ii(comparison)
%PRINT_TABLE_II  打印 Table II 风格的算法对比汇总表
%
% 输入:
%   comparison : comparison 结构体

    algorithms = comparison.algorithms;
    num_algorithms = length(algorithms);

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('              Table II: Algorithm Comparison Summary\n');
    fprintf('============================================================\n');
    fprintf('| Algorithm       | Feasible | J_total  | Time(ms) | Constraints |\n');
    fprintf('|-----------------|----------|----------|----------|-------------|\n');

    for alg_idx = 1:num_algorithms
        alg = algorithms{alg_idx};
        alg_results = get_alg_results(comparison, alg_idx);
        n_success = sum([alg_results.success]);

        if n_success > 0
            success_mask = [alg_results.success];
            success_results = alg_results(success_mask);

            traj_costs = arrayfun(@(r) r.metrics.trajectoryCost, success_results);
            mean_J = nanmean(traj_costs);

            solve_times = arrayfun(@(r) r.metrics.meanSolveTime, success_results) * 1000;
            mean_time = nanmean(solve_times);

            ws_viol = nanmean(arrayfun(@(r) r.metrics.wheelSpeedViolationRatio, success_results));
            sr_viol = nanmean(arrayfun(@(r) r.metrics.steeringRateViolationRatio, success_results));

            feas_str = sprintf('%d/%d', n_success, comparison.numSeeds);
            J_str = sprintf('%.2f', mean_J);
            time_str = sprintf('%.1f', mean_time);

            if ws_viol > 0.01 || sr_viol > 0.01
                constr_str = 'Violated';
            else
                constr_str = 'Satisfied';
            end
        else
            feas_str = sprintf('0/%d', comparison.numSeeds);
            J_str = 'N/A';
            time_str = 'N/A';
            constr_str = 'N/A';
        end

        fprintf('| %-15s | %8s | %8s | %8s | %11s |\n', ...
            alg, feas_str, J_str, time_str, constr_str);
    end

    fprintf('============================================================\n');
end
