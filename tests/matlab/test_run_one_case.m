function test_run_one_case()
% TEST_RUN_ONE_CASE  MATLAB smoke test for run_one_case / run_closed_loop
%
% 用法:
%   cd D:\PROJECT\RSS_V2
%   test_run_one_case
%
% 或从 tests/matlab/ 目录:
%   cd tests/matlab
%   test_run_one_case

    fprintf('===== test_run_one_case 开始 =====\n');
    passed = 0;
    failed = 0;

    %% Test 1: run_closed_loop 基本功能
    fprintf('\n[Test 1] run_closed_loop 默认配置...\n');
    try
        cfg = config();
        cfg.live_plot = false;
        cfg.case_id = 'test_smoke';
        cfg.output_dir = '';

        result = run_closed_loop(cfg);

        % 检查必需字段
        required_fields = {'reference', 'state_history', 'body_velocity_history', ...
            'control_history', 'wheel_speed_history', 'wheel_angle_history', ...
            'steering_rate_history', 'solver_time_history', 'cvx_status_history', ...
            'trajectory_cost', 'position_rmse', 'completed_steps', 'success', ...
            'config', 'wall_time'};

        for i = 1:length(required_fields)
            if ~isfield(result, required_fields{i})
                error('MissingField: %s', required_fields{i});
            end
        end

        % 检查完成步数
        assert(result.completed_steps == 100, ...
            sprintf('completed_steps should be 100, got %d', result.completed_steps));

        % 检查 success
        assert(result.success == true, 'success should be true');

        % 检查 state_history 大小
        [rows, cols] = size(result.state_history);
        assert(rows == 100 && cols == 2, ...
            sprintf('state_history should be 100x2, got %dx%d', rows, cols));

        fprintf('  PASS: run_closed_loop 基本功能正常\n');
        fprintf('  RMSE: %.10f\n', result.position_rmse);
        fprintf('  Cost: %.10f\n', result.trajectory_cost);
        passed = passed + 1;
    catch ME
        fprintf('  FAIL: %s\n', ME.message);
        failed = failed + 1;
    end

    %% Test 2: run_one_case JSON 接口
    fprintf('\n[Test 2] run_one_case JSON 接口...\n');
    try
        % 创建临时 JSON 配置
        cfg_json = struct();
        cfg_json.case_id = 'test_json';
        cfg_json.solver = 'sdpt3';
        cfg_json.num_steps = 100;
        cfg_json.live_plot = false;
        cfg_json.save_figures = false;
        cfg_json.save_full_log = false;
        cfg_json.output_dir = fullfile(tempdir, 'rss_test_json');

        json_path = fullfile(tempdir, 'test_config.json');
        fid = fopen(json_path, 'w');
        fprintf(fid, '%s', jsonencode(cfg_json));
        fclose(fid);

        summary_json = run_one_case(json_path);
        summary = jsondecode(summary_json);

        % 检查必需字段
        assert(isfield(summary, 'success'), 'summary missing success');
        assert(isfield(summary, 'completed_steps'), 'summary missing completed_steps');
        assert(isfield(summary, 'position_rmse'), 'summary missing position_rmse');
        assert(isfield(summary, 'result_file'), 'summary missing result_file');

        % 检查 result.mat 存在
        assert(exist(summary.result_file, 'file') == 2, ...
            sprintf('result.mat not found: %s', summary.result_file));

        % 检查 summary.json 存在
        summary_path = fullfile(cfg_json.output_dir, 'summary.json');
        assert(exist(summary_path, 'file') == 2, ...
            'summary.json not found');

        % 检查 success
        assert(summary.success == true, 'summary.success should be true');
        assert(summary.completed_steps == 100, ...
            sprintf('completed_steps should be 100, got %d', summary.completed_steps));

        fprintf('  PASS: run_one_case JSON 接口正常\n');
        fprintf('  Output dir: %s\n', cfg_json.output_dir);
        passed = passed + 1;
    catch ME
        fprintf('  FAIL: %s\n', ME.message);
        failed = failed + 1;
    end

    %% Test 3: 失败时返回合法 JSON
    fprintf('\n[Test 3] 失败场景 JSON 合法性...\n');
    try
        % 用一个会导致失败的配置 (num_steps = 0)
        cfg_bad = struct();
        cfg_bad.case_id = 'test_fail';
        cfg_bad.num_steps = 0;
        cfg_bad.live_plot = false;
        cfg_bad.output_dir = fullfile(tempdir, 'rss_test_fail');

        json_path = fullfile(tempdir, 'test_config_bad.json');
        fid = fopen(json_path, 'w');
        fprintf(fid, '%s', jsonencode(cfg_bad));
        fclose(fid);

        summary_json = run_one_case(json_path);

        % 验证返回的是合法 JSON
        summary = jsondecode(summary_json);
        assert(isfield(summary, 'success'), 'failure summary missing success');
        assert(isfield(summary, 'failure_reason'), 'failure summary missing failure_reason');

        fprintf('  PASS: 失败时返回合法 JSON\n');
        passed = passed + 1;
    catch ME
        fprintf('  FAIL: %s\n', ME.message);
        failed = failed + 1;
    end

    %% Test 4: control_RSS diagnostics 输出
    fprintf('\n[Test 4] control_RSS diagnostics 输出...\n');
    try
        % 直接调用 control_RSS 检查第 4 输出
        addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'algorithms', 'RSS_proposed'));
        cleanup = onCleanup(@() rmpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'algorithms', 'RSS_proposed')));
        clear functions;

        params = config();
        path = bezier_path(params.ctrl_pts, params.num_path_pts);
        state = [0.05, 0.1, 0.2];
        last_vel = [0.01; 0.01; 0.01];

        [u, new_state_dot, velocity, diag] = control_RSS(path, 1, last_vel, state);

        assert(isstruct(diag), 'diagnostics should be a struct');
        assert(isfield(diag, 'iterations'), 'diagnostics missing iterations');
        assert(isfield(diag, 'has_feasible'), 'diagnostics missing has_feasible');
        assert(isfield(diag, 'total_solve_time'), 'diagnostics missing total_solve_time');

        fprintf('  PASS: control_RSS diagnostics 输出正常\n');
        fprintf('  has_feasible: %d\n', diag.has_feasible);
        fprintf('  total_solve_time: %.6f\n', diag.total_solve_time);
        passed = passed + 1;
    catch ME
        fprintf('  FAIL: %s\n', ME.message);
        failed = failed + 1;
    end

    %% Summary
    fprintf('\n===== test_run_one_case 结果: %d passed, %d failed =====\n', passed, failed);
    if failed > 0
        error('test_run_one_case: %d tests failed', failed);
    end
end
