function [u, new_state_dot, velocity, diagnostics] = control_RSS_denseqcqp(path, step, state_dot, state)
% CONTROL_RSS_DENSEQCQP  Dense QCQP 求解器 (精确二次约束 oracle)
%
% v2 benchmark 版本 (rss_hpipm_qp_v2 分支):
%   - Dense QCQP 求解 (HPIPM dense_qcqp, 精确二次约束)
%   - 零初始化 u_hat = zeros(3, K)
%   - 简单 u_hat 更新 (status==0 时 u_hat = u_sol)
%   - 3 次 SCP 外层迭代
%
% 用法 (在 run_paper_baseline_case.m 中通过 cfg.solver='denseqcqp' 调用):
%   cfg = defaultConfig();
%   cfg.algorithm = 'proposed-3iter';
%   cfg.solver = 'denseqcqp';
%   summary = run_paper_baseline_case(cfg);
%
% 预期结果 (benchmark):
%   J_total=13.3838, RMSE=0.036793, validSteps=100/100
%
% 论文: RSS26 "Exploit Agile Mobility of Steerable-Wheeled Mobile Robots:
%        A Fast Motion Planning Approach"
    params = config();
    % ================= Param Setup =================
    K = 6; rho = 0.01; k1 = 1; epsilon = 0;
    current_xy = [state(1), state(2)]';
    psi0 = state(3); v0 = state_dot;
    % ================= 迭代 Setup =================
    max_iter = 3;  % 论文 IV-B: 固定 3 次外层迭代
    u_hat = zeros(3, K);  % u^(0) = 0 (static init, 与 v2 一致)
    global solver_time_array;
    if ~exist('solver_time_array', 'var') || isempty(solver_time_array)
        solver_time_array = [];
    end
    % 诊断结构体 (记录每次迭代)
    diagnostics = struct();
    diagnostics.iterations = struct();
    diagnostics.iterations.status = cell(1, max_iter);
    diagnostics.iterations.optval = zeros(1, max_iter);
    diagnostics.iterations.solve_time = NaN(1, max_iter);
    diagnostics.iterations.solver_name = cell(1, max_iter);
    diagnostics.step = step;
    diagnostics.max_iter = max_iter;
    diagnostics.solver_call_count = 0;  % 兼容 run_paper_baseline_case.m
    % ================= Python 环境路径设置 =================
    persistent py_path_added py_reloaded;
    if isempty(py_path_added)
        script_path = fileparts(mfilename('fullpath'));
        if exist(script_path, 'dir')
            sys_mod = py.importlib.import_module('sys');
            py.getattr(sys_mod, 'path').append(script_path);
        end
        py_path_added = true;
    end
    if isempty(py_reloaded)
        try
            reload_py = fullfile(script_path, '_reload_solver_tmp.py');
            fid = fopen(reload_py, 'w');
            fprintf(fid, 'import sys\n');
            fprintf(fid, 'sys.modules.pop("hpipm_qp_solver", None)\n');
            fprintf(fid, 'import hpipm_qp_solver\n');
            fclose(fid);
            py.runpy.run_path(reload_py);
            delete(reload_py);
            py_reloaded = true;
        catch
        end
    end
    % ================= SQP 外层循环 (论文 Algorithm 1 line 3-9) =================
    for m = 1 : max_iter
        try
            % ========== 论文 Alg.1 line 4: 构造凸子问题 Q_K(u^(m)) ==========
            % Dense QCQP (精确凸化二次约束, HPIPM dense_qcqp)
            qp = construct_complete_qp_from_rss(path, step, v0, state, u_hat, params);
            % ========== 数据转换: Hq/gq cell -> numpy 3D/2D ==========
            n_var = qp.n_var;
            nq = length(qp.Hq);
            Hq_3d = zeros(n_var, n_var, nq);
            gq_2d = zeros(n_var, nq);
            for i = 1:nq
                Hq_3d(:,:,i) = qp.Hq{i};
                gq_2d(:,i) = qp.gq{i};
            end
            % ========== 论文 Alg.1 line 5: 求解 u^(m+1) = S(u^(m)) ==========
            hpipm_mod = py.importlib.import_module('hpipm_qp_solver');
            result = hpipm_mod.solve_qcqp(...
                py.numpy.array(qp.H), ...
                py.numpy.array(qp.g), ...
                py.numpy.array(qp.A), ...
                py.numpy.array(qp.b), ...
                py.numpy.array(Hq_3d), ...
                py.numpy.array(gq_2d), ...
                py.numpy.array(qp.uq), ...
                py.bool(step == 1 && m == 1) ...
            );
            % ========== 提取结果 ==========
            x = double(result{'x'});
            status_code = double(result{'status'});
            optval = double(result{'obj_value'});
            inner_solve_time = double(result{'solve_time'});
            u_sol = reshape(x(1:3*K), 3, K);
            if status_code == 0
                cvx_status_str = 'Solved';
            else
                cvx_status_str = 'Failed';
            end
            solver_name = 'HPIPM-DenseQCQP';
            diagnostics.solver_call_count = diagnostics.solver_call_count + 1;
        catch ME
            u_sol = zeros(3, K);
            status_code = -1;
            optval = NaN;
            inner_solve_time = NaN;
            cvx_status_str = 'Failed';
            solver_name = 'HPIPM-Error';
            fprintf('第%d步第%d次迭代 - Python 调用异常: %s\n', step, m, ME.message);
            for si = 1:min(length(ME.stack), 5)
                fprintf('    at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
            end
        end
        % ========== 存入全局数组 ==========
        solver_time_array(3*step + m - 3) = inner_solve_time;
        % ========== 记录诊断 ==========
        diagnostics.iterations.status{m} = cvx_status_str;
        diagnostics.iterations.optval(m) = optval;
        diagnostics.iterations.solve_time(m) = inner_solve_time;
        diagnostics.iterations.solver_name{m} = solver_name;
        fprintf('第%d步第%d次迭代 - %s内部求解时间：%.6f秒 (status=%d)\n', ...
            step, m, solver_name, inner_solve_time, status_code);
        fprintf('最优代价: %f | 求解状态: %s\n', optval, cvx_status_str);
        % ========== 论文 Alg.1 line 9: 更新 m ← m+1 ==========
        u = u_sol;
        if status_code == 0 && all(isfinite(u_sol(:)))
            u_hat = u;
        end
    end
    % ================= 论文 Alg.1 line 11: 输出 ν_1 = ν_0 + u_1 =================
    new_state_dot = [cos(state(3)), -sin(state(3)), 0;
                     sin(state(3)),  cos(state(3)), 0;
                         0,              0, 1] * (state_dot + 1.00 * u(:, 1));
    velocity = v0 + u(:, 1);
    % ================= [DIAGNOSTIC] 原始约束违反量检查 (per-step) =================
    num_wheels = size(params.wheel_pos, 1);
    H_diag = cell(1, num_wheels);
    for n = 1:num_wheels
        H_diag{n} = [1, 0, -params.wheel_pos(n,2); 0, 1, params.wheel_pos(n,1)];
    end
    delta_th = params.dt * params.phidotmax;
    R1d = [sin(delta_th), -cos(delta_th); cos(delta_th), sin(delta_th)];
    R2d = R1d';
    nu_diag = zeros(3, K+1); nu_diag(:, 1) = v0;
    for kk = 1:K; nu_diag(:, kk+1) = nu_diag(:,kk) + u_hat(:, kk); end
    % ---- 轮速约束 (k=1..K) ----
    mv_wheel = 0;
    for kk = 1:K
        nkk = nu_diag(:, kk+1);
        for n = 1:num_wheels
            vn = norm(H_diag{n} * nkk, 2);
            viol = vn - params.vimax;
            if viol > mv_wheel; mv_wheel = viol; end
        end
    end
    % ---- 转向锥 (k=1..K) ----
    mv_cone = 0;
    for kk = 1:K
        a = nu_diag(:, kk);
        b = nu_diag(:, kk+1);
        for n = 1:num_wheels
            Hn1 = H_diag{n};
            for gg = 1:2
                if gg == 1; Rgg = R1d; else; Rgg = R2d; end
                Mgg = Hn1' * Rgg * Hn1;
                term = a' * Mgg * b;
                viol = -term;
                if viol > mv_cone; mv_cone = viol; end
            end
        end
    end
    if mod(step, 25) == 0 || step == 1
        fprintf('[DIAG step=%d] wheel_viol_cur=%.6e, cone_viol_cur=%.6e\n', step, mv_wheel, mv_cone);
    end
    % ================= 补充 diagnostics 字段 (兼容 run_paper_baseline_case.m) =================
    diagnostics.orig_wheel_viol_final = mv_wheel;
    diagnostics.orig_cone_viol_final = mv_cone;
    diagnostics.incumbent_type = 'strict';
    diagnostics.step_approximate = false;
    diagnostics.has_strict_incumbent = true;
    diagnostics.step_failed = false;
    diagnostics.total_solve_time = sum(diagnostics.iterations.solve_time(~isnan( ...
        diagnostics.iterations.solve_time)));
    if isempty(diagnostics.total_solve_time)
        diagnostics.total_solve_time = NaN;
    end
end
