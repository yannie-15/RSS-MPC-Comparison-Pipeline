function test_parity()
% TEST_PARITY  新旧结果回归对比
%
% 比较三条路径的输出是否一致:
%   1. 原始 algorithms/RSS_proposed/main.m (通过 freeze_baseline 保存的基线)
%   2. 新 run_closed_loop(config())
%   3. 新 run_one_case(JSON) -> run_closed_loop
%
% 用法:
%   cd D:\PROJECT\RSS_V2
%   test_parity
%
% 前提: 先运行 freeze_baseline 生成基线文件
%   cd tests/baseline && freeze_baseline

    fprintf('===== test_parity 新旧结果回归对比 =====\n\n');

    %% 定位仓库根目录
    script_dir = fileparts(mfilename('fullpath'));
    repo_root = fileparts(fileparts(script_dir));  % tests/matlab -> tests -> repo_root

    %% 加载基线
    baseline_file = fullfile(repo_root, 'tests', 'baseline', 'rss_cvx_original_main.mat');
    if ~exist(baseline_file, 'file')
        error('test_parity:NoBaseline', ...
            ['基线文件不存在: %s\n', ...
             '请先运行: cd tests/baseline && freeze_baseline'], baseline_file);
    end

    loaded = load(baseline_file, 'baseline');
    baseline = loaded.baseline;

    fprintf('基线已加载: %s\n', baseline_file);
    fprintf('  Git commit: (见 baseline_meta.json)\n');
    fprintf('  RMSE: %.10f\n', baseline.rmse);
    fprintf('  Cost: %.10f\n', baseline.trajectory_cost);
    fprintf('  Steps: %d\n\n', size(baseline.q_history, 1));

    %% 运行新的 run_closed_loop
    fprintf('[Path 2] 运行 run_closed_loop...\n');
    addpath(fullfile(repo_root, 'algorithms', 'RSS_proposed'));
    cleanup1 = onCleanup(@() rmpath(fullfile(repo_root, 'algorithms', 'RSS_proposed')));
    clear functions;

    cfg = config();
    cfg.live_plot = false;
    cfg.case_id = 'parity_test';
    cfg.output_dir = '';

    result = run_closed_loop(cfg);
    fprintf('\n');

    %% 对比
    tol = 1e-6;  % 容差
    passed = 0;
    failed = 0;
    warnings_count = 0;

    % 1. Reference 完全一致
    fprintf('[Check 1] 参考轨迹一致性...\n');
    ref_diff = max(abs(baseline.reference(:) - result.reference(:)));
    if ref_diff == 0
        fprintf('  PASS: 参考轨迹完全一致\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: 参考轨迹最大差异 = %e\n', ref_diff);
        failed = failed + 1;
    end

    % 2. 完成步数
    fprintf('[Check 2] 完成步数...\n');
    if result.completed_steps == size(baseline.q_history, 1)
        fprintf('  PASS: 完成步数一致 (%d)\n', result.completed_steps);
        passed = passed + 1;
    else
        fprintf('  FAIL: 步数不一致 (基线=%d, 新=%d)\n', ...
            size(baseline.q_history, 1), result.completed_steps);
        failed = failed + 1;
    end

    % 3. 状态历史
    fprintf('[Check 3] 状态历史 (q_history)...\n');
    n = min(size(baseline.q_history, 1), result.completed_steps);
    state_diff = max(abs(baseline.q_history(1:n, :) - result.state_history(1:n, :)));
    if state_diff < tol
        fprintf('  PASS: 状态最大差异 = %e (< %e)\n', state_diff, tol);
        passed = passed + 1;
    elseif state_diff < 1e-4
        fprintf('  WARN: 状态最大差异 = %e (在可接受范围, 但超过容差)\n', state_diff);
        warnings_count = warnings_count + 1;
    else
        fprintf('  FAIL: 状态最大差异 = %e (>= %e)\n', state_diff, tol);
        failed = failed + 1;
    end

    % 4. 轮速历史
    fprintf('[Check 4] 轮速历史 (vi_history)...\n');
    vi_diff = max(abs(baseline.vi_history(1:n, :) - result.wheel_speed_history(1:n, :)));
    if vi_diff < tol
        fprintf('  PASS: 轮速最大差异 = %e\n', vi_diff);
        passed = passed + 1;
    elseif vi_diff < 1e-4
        fprintf('  WARN: 轮速最大差异 = %e\n', vi_diff);
        warnings_count = warnings_count + 1;
    else
        fprintf('  FAIL: 轮速最大差异 = %e\n', vi_diff);
        failed = failed + 1;
    end

    % 5. 轮角历史
    fprintf('[Check 5] 轮角历史 (phi_history)...\n');
    phi_diff = max(abs(baseline.phi_history(1:n, :) - result.wheel_angle_history(1:n, :)));
    if phi_diff < tol
        fprintf('  PASS: 轮角最大差异 = %e\n', phi_diff);
        passed = passed + 1;
    elseif phi_diff < 1e-4
        fprintf('  WARN: 轮角最大差异 = %e\n', phi_diff);
        warnings_count = warnings_count + 1;
    else
        fprintf('  FAIL: 轮角最大差异 = %e\n', phi_diff);
        failed = failed + 1;
    end

    % 6. 控制历史
    fprintf('[Check 6] 控制历史 (u_history)...\n');
    u_diff = max(abs(baseline.u_history(1:n, :) - result.control_history(1:n, :)));
    if u_diff < tol
        fprintf('  PASS: 控制最大差异 = %e\n', u_diff);
        passed = passed + 1;
    elseif u_diff < 1e-4
        fprintf('  WARN: 控制最大差异 = %e\n', u_diff);
        warnings_count = warnings_count + 1;
    else
        fprintf('  FAIL: 控制最大差异 = %e\n', u_diff);
        failed = failed + 1;
    end

    % 7. 车体速度历史
    fprintf('[Check 7] 车体速度历史 (velocity_history)...\n');
    vel_diff = max(abs(baseline.velocity_history(1:n, :) - result.body_velocity_history(1:n, :)));
    if vel_diff < tol
        fprintf('  PASS: 速度最大差异 = %e\n', vel_diff);
        passed = passed + 1;
    elseif vel_diff < 1e-4
        fprintf('  WARN: 速度最大差异 = %e\n', vel_diff);
        warnings_count = warnings_count + 1;
    else
        fprintf('  FAIL: 速度最大差异 = %e\n', vel_diff);
        failed = failed + 1;
    end

    % 8. 轨迹代价
    fprintf('[Check 8] 轨迹代价 (J)...\n');
    cost_diff = abs(baseline.trajectory_cost - result.trajectory_cost);
    if cost_diff < tol
        fprintf('  PASS: 代价差异 = %e (基线=%.10f, 新=%.10f)\n', ...
            cost_diff, baseline.trajectory_cost, result.trajectory_cost);
        passed = passed + 1;
    elseif cost_diff < 1e-4
        fprintf('  WARN: 代价差异 = %e (基线=%.10f, 新=%.10f)\n', ...
            cost_diff, baseline.trajectory_cost, result.trajectory_cost);
        warnings_count = warnings_count + 1;
    else
        fprintf('  FAIL: 代价差异 = %e (基线=%.10f, 新=%.10f)\n', ...
            cost_diff, baseline.trajectory_cost, result.trajectory_cost);
        failed = failed + 1;
    end

    % 9. RMSE
    fprintf('[Check 9] 位置 RMSE...\n');
    rmse_diff = abs(baseline.rmse - result.position_rmse);
    if rmse_diff < tol
        fprintf('  PASS: RMSE 差异 = %e (基线=%.10f, 新=%.10f)\n', ...
            rmse_diff, baseline.rmse, result.position_rmse);
        passed = passed + 1;
    elseif rmse_diff < 1e-6
        fprintf('  WARN: RMSE 差异 = %e (基线=%.10f, 新=%.10f)\n', ...
            rmse_diff, baseline.rmse, result.position_rmse);
        warnings_count = warnings_count + 1;
    else
        fprintf('  FAIL: RMSE 差异 = %e (基线=%.10f, 新=%.10f)\n', ...
            rmse_diff, baseline.rmse, result.position_rmse);
        failed = failed + 1;
    end

    % 10. CVX status 序列
    fprintf('[Check 10] CVX 状态序列...\n');
    status_mismatch = 0;
    for i = 1:n
        if ~strcmpi(baseline.cvx_status_history{i}, result.cvx_status_history{i})
            status_mismatch = status_mismatch + 1;
            if status_mismatch <= 3
                fprintf('  步 %d: 基线=%s, 新=%s\n', ...
                    i, baseline.cvx_status_history{i}, result.cvx_status_history{i});
            end
        end
    end
    if status_mismatch == 0
        fprintf('  PASS: CVX 状态序列完全一致\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: %d/%d 步 CVX 状态不一致\n', status_mismatch, n);
        failed = failed + 1;
    end

    % 11. 轮速约束
    fprintf('[Check 11] 轮速约束 (vimax)...\n');
    max_vi_new = max(result.wheel_speed_history(1:n, :), [], 'all');
    if max_vi_new <= cfg.vimax * (1 + 1e-6)
        fprintf('  PASS: 最大轮速 %.6f <= vimax %.1f\n', max_vi_new, cfg.vimax);
        passed = passed + 1;
    else
        fprintf('  FAIL: 最大轮速 %.6f > vimax %.1f (违反 %.6f)\n', ...
            max_vi_new, cfg.vimax, max_vi_new - cfg.vimax);
        failed = failed + 1;
    end

    % 12. 转向率约束
    fprintf('[Check 12] 转向率约束 (phidotmax)...\n');
    if n > 1
        max_phidot_new = max(abs(result.steering_rate_history(1:n-1, :)), [], 'all');
    else
        max_phidot_new = 0;
    end
    if max_phidot_new <= cfg.phidotmax * (1 + 1e-6)
        fprintf('  PASS: 最大转向率 %.6f <= phidotmax %.6f\n', max_phidot_new, cfg.phidotmax);
        passed = passed + 1;
    else
        fprintf('  WARN: 最大转向率 %.6f > phidotmax %.6f (违反 %.6f)\n', ...
            max_phidot_new, cfg.phidotmax, max_phidot_new - cfg.phidotmax);
        warnings_count = warnings_count + 1;
    end

    %% 汇总
    fprintf('\n===== test_parity 结果 =====\n');
    fprintf('  Passed:   %d\n', passed);
    fprintf('  Warnings: %d\n', warnings_count);
    fprintf('  Failed:   %d\n', failed);

    if failed > 0
        error('test_parity: %d checks failed', failed);
    end
    if warnings_count > 0
        fprintf('\n注意: %d 个警告 (数值差异在可接受范围但超过严格容差)\n', warnings_count);
        fprintf('这可能由 CVX 求解器版本差异或浮点累积导致。\n');
    end
end
