function test_ocp_qcqp_golden_100steps()
% TEST_OCP_QCQP_GOLDEN_100STEPS
% Test D: 100 步闭环 Golden 测试 (Dense QCQP vs OCP QCQP)
%
% 使用相同配置分别运行:
%   cfg.solver = 'denseqcqp'  (Golden oracle)
%   cfg.solver = 'ocpqcqp'    (新原生 OCP QCQP)
%
% 验收:
%   validSteps = 100/100 (两者)
%   strict/approx/none = 100/0/0 (两者)
%   solver_call_count = 3 per step (两者)
%   max original wheel violation < 1e-8 (两者)
%   max original cone violation < 1e-8 (两者)
%   abs(RMSE_ocp - RMSE_dense) < 1e-6
%   abs(J_total_ocp - J_total_dense) < 1e-6
%   max closed-loop state difference < 1e-6
%   max executed-control difference < 1e-6
%
% 若出现第一处分叉, 必须报告 step 和误差, 不得放宽阈值.
%
% 用法:
%   cd d:\PROJECT\RSS-MPC-Comparison-Pipeline-rss_hpipm
%   rmpath('D:\PROJECT\RSS_V2\matlab'); rmpath('D:\Projects\RSS\matlab');
%   clear all; close all; clc;
%   run('core/setup_paths.m'); addpath('paper_reproduction');
%   addpath('algorithms/RSS_proposed/tests');
%   setenv('HPIPM_OCP_QCQP_MODE','speed');
%   test_ocp_qcqp_golden_100steps

    % 测试文件位于 algorithms/RSS_proposed/tests/, 需添加父目录以访问 control/config 等
    addpath(fileparts(fileparts(mfilename('fullpath'))));

    fprintf('============== Test D: 100 步闭环 Golden 对齐 ==============\n');

    %% ===== 1. Dense QCQP (Golden oracle) =====
    fprintf('\n========== [1/2] Dense QCQP (Golden oracle) ==========\n');
    cfg_dense = defaultConfig();
    cfg_dense.algorithm = 'proposed-3iter';
    cfg_dense.solver = 'denseqcqp';
    scenario = struct('name', 'paper_fixed', 'id', 0);

    summary_dense = run_paper_baseline_case(cfg_dense, scenario);

    %% ===== 2. OCP QCQP =====
    fprintf('\n========== [2/2] OCP QCQP ==========\n');
    cfg_ocp = defaultConfig();
    cfg_ocp.algorithm = 'proposed-3iter';
    cfg_ocp.solver = 'ocpqcqp';

    summary_ocp = run_paper_baseline_case(cfg_ocp, scenario);

    %% ===== 3. 对比 =====
    fprintf('\n========== 对比结果 ==========\n');

    % --- 基本指标 ---
    dense_valid = summary_dense.solverInfo.nValidSteps;
    ocp_valid = summary_ocp.solverInfo.nValidSteps;
    dense_success = summary_dense.success;
    ocp_success = summary_ocp.success;

    % strict/approx/none
    dense_cnt_strict = sum(strcmp(summary_dense.constraintInfo.incumbentTypePerStep, 'strict'));
    dense_cnt_approx = sum(strcmp(summary_dense.constraintInfo.incumbentTypePerStep, 'approximate'));
    dense_cnt_none = sum(strcmp(summary_dense.constraintInfo.incumbentTypePerStep, 'none'));
    ocp_cnt_strict = sum(strcmp(summary_ocp.constraintInfo.incumbentTypePerStep, 'strict'));
    ocp_cnt_approx = sum(strcmp(summary_ocp.constraintInfo.incumbentTypePerStep, 'approximate'));
    ocp_cnt_none = sum(strcmp(summary_ocp.constraintInfo.incumbentTypePerStep, 'none'));

    % solver_call_count
    dense_calls = summary_dense.solverInfo.solverCallCounts;
    ocp_calls = summary_ocp.solverInfo.solverCallCounts;
    dense_all_3 = all(dense_calls == 3);
    ocp_all_3 = all(ocp_calls == 3);

    % 约束违反
    dense_max_wheel = summary_dense.constraintInfo.maxWheelViol;
    dense_max_cone = summary_dense.constraintInfo.maxConeViol;
    ocp_max_wheel = summary_ocp.constraintInfo.maxWheelViol;
    ocp_max_cone = summary_ocp.constraintInfo.maxConeViol;

    % RMSE 和 J_total
    dense_RMSE = summary_dense.metrics.rmse;
    ocp_RMSE = summary_ocp.metrics.rmse;
    dense_J = summary_dense.metrics.J_total;
    ocp_J = summary_ocp.metrics.J_total;

    % 闭环轨迹和控制差异
    n_compare = min(dense_valid, ocp_valid);
    if n_compare > 0
        state_diff = max(abs(summary_ocp.states(:, 1:n_compare+1) - summary_dense.states(:, 1:n_compare+1)), [], 'all');
        u_diff = max(abs(summary_ocp.executedU(:, 1:n_compare) - summary_dense.executedU(:, 1:n_compare)), [], 'all');
    else
        state_diff = inf;
        u_diff = inf;
    end

    % 性能
    dense_median_time = summary_dense.metrics.medianSolveTime;
    ocp_median_time = summary_ocp.metrics.medianSolveTime;

    % --- 打印 ---
    fprintf('\n--- 基本指标 ---\n');
    fprintf('                    Dense QCQP          OCP QCQP\n');
    fprintf('success:            %-20s%-20s\n', num2str(dense_success), num2str(ocp_success));
    fprintf('validSteps:         %-20d%-20d\n', dense_valid, ocp_valid);
    fprintf('strict/approx/none: %d/%d/%-13d%6d/%d/%d\n', ...
        dense_cnt_strict, dense_cnt_approx, dense_cnt_none, ...
        ocp_cnt_strict, ocp_cnt_approx, ocp_cnt_none);
    fprintf('solver_call all==3: %-20s%-20s\n', num2str(dense_all_3), num2str(ocp_all_3));
    fprintf('\n--- 约束违反 ---\n');
    fprintf('maxWheel:           %-20.6e%-20.6e\n', dense_max_wheel, ocp_max_wheel);
    fprintf('maxCone:            %-20.6e%-20.6e\n', dense_max_cone, ocp_max_cone);
    fprintf('\n--- RMSE / J_total ---\n');
    fprintf('RMSE:               %-20.6f%-20.6f\n', dense_RMSE, ocp_RMSE);
    fprintf('J_total:            %-20.6f%-20.6f\n', dense_J, ocp_J);
    fprintf('\n--- 对齐差异 ---\n');
    fprintf('abs(RMSE diff):     %.6e  (target < 1e-6)\n', abs(ocp_RMSE - dense_RMSE));
    fprintf('abs(J_total diff):  %.6e  (target < 1e-6)\n', abs(ocp_J - dense_J));
    fprintf('max state diff:     %.6e  (target < 1e-6)\n', state_diff);
    fprintf('max control diff:   %.6e  (target < 1e-6)\n', u_diff);
    fprintf('\n--- 性能 ---\n');
    fprintf('median solve time:  %-20.6f%-20.6f (seconds)\n', dense_median_time, ocp_median_time);

    % --- 找第一处分叉 ---
    if n_compare > 0
        first_diverge = 0;
        for k = 1:n_compare
            if max(abs(summary_ocp.states(:, k+1) - summary_dense.states(:, k+1))) > 1e-8
                first_diverge = k;
                break;
            end
        end
        if first_diverge > 0
            fprintf('\n--- 第一处轨迹分叉 ---\n');
            fprintf('step %d: state_diff=%.6e\n', first_diverge, ...
                max(abs(summary_ocp.states(:, first_diverge+1) - summary_dense.states(:, first_diverge+1))));
        end
    end

    % --- 验收 ---
    fprintf('\n========== 验收 ==========\n');
    pass = true;
    fail_reasons = {};

    if dense_valid ~= 100
        pass = false; fail_reasons{end+1} = sprintf('Dense validSteps=%d != 100', dense_valid);
    end
    if ocp_valid ~= 100
        pass = false; fail_reasons{end+1} = sprintf('OCP validSteps=%d != 100', ocp_valid);
    end
    if ~(ocp_cnt_strict == 100 && ocp_cnt_approx == 0 && ocp_cnt_none == 0)
        pass = false; fail_reasons{end+1} = sprintf('OCP strict/approx/none=%d/%d/%d != 100/0/0', ...
            ocp_cnt_strict, ocp_cnt_approx, ocp_cnt_none);
    end
    if ~ocp_all_3
        pass = false; fail_reasons{end+1} = sprintf('OCP solver_call_count not all 3');
    end
    if ocp_max_wheel >= 1e-8
        pass = false; fail_reasons{end+1} = sprintf('OCP maxWheel=%.6e >= 1e-8', ocp_max_wheel);
    end
    if ocp_max_cone >= 1e-8
        pass = false; fail_reasons{end+1} = sprintf('OCP maxCone=%.6e >= 1e-8', ocp_max_cone);
    end
    if abs(ocp_RMSE - dense_RMSE) >= 1e-6
        pass = false; fail_reasons{end+1} = sprintf('|RMSE diff|=%.6e >= 1e-6', abs(ocp_RMSE - dense_RMSE));
    end
    if abs(ocp_J - dense_J) >= 1e-6
        pass = false; fail_reasons{end+1} = sprintf('|J_total diff|=%.6e >= 1e-6', abs(ocp_J - dense_J));
    end
    if state_diff >= 1e-6
        pass = false; fail_reasons{end+1} = sprintf('max state diff=%.6e >= 1e-6', state_diff);
    end
    if u_diff >= 1e-6
        pass = false; fail_reasons{end+1} = sprintf('max control diff=%.6e >= 1e-6', u_diff);
    end

    if pass
        fprintf('结论: PASS (100 步 Golden 对齐)\n');
    else
        fprintf('结论: FAIL\n');
        for i = 1:numel(fail_reasons)
            fprintf('  原因: %s\n', fail_reasons{i});
        end
    end
    fprintf('==============================================================\n');
end
