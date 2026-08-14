function [u, new_state_dot, velocity, diagnostics] = control_RSS_ocpqcqp(path, step, state_dot, state)
% CONTROL_RSS_OCPQCQP  HPIPM 原生 OCP QCQP 求解器 (精确凸二次约束)
%
% Golden 对齐模式: 严格复制 control_RSS_denseqcqp.m 的 outer 语义.
%   - u_hat = zeros(3, K) (零初始化, 无 warm start, 无 cut_points)
%   - 3 次 outer 迭代, 每次恰好 1 次 solve_ocp_qcqp (solver_call_count=1)
%   - status==0 且解有限: u_hat = u_sol
%   - 否则: 保留上一轮 u_hat, 不推进
%   - 所有成功步标记为 strict (无 approximate)
%
% 与 Dense QCQP 的区别:
%   - Dense QCQP: dense 化后用 hpipm_dense_qcqp 求解 (36 维, 72 约束)
%   - OCP QCQP:   利用 OCP 块三对角结构, 用 hpipm_ocp_qcqp 求解 (7 stages, 72 约束)
%   - 两者数学等价, 数值容差内应得到相同的 u, v, obj, 约束可行性
%
% 用法 (在 run_paper_baseline_case.m 中通过 cfg.solver='ocpqcqp' 调用):
%   cfg = defaultConfig();
%   cfg.algorithm = 'proposed-3iter';
%   cfg.solver = 'ocpqcqp';
%   summary = run_paper_baseline_case(cfg);

    params = config();

    % ================= Param Setup =================
    K = 6; rho = 0.01; k1 = 1; epsilon = 0;
    current_xy = [state(1), state(2)]';
    psi0 = state(3); v0 = state_dot;

    % ================= 迭代 Setup =================
    max_iter = 3;  % 论文 IV-B: 固定 3 次外层迭代
    u_hat = zeros(3, K);  % u^(0) = 0 (static init, 与 Dense QCQP 一致)

    global solver_time_array;
    if ~exist('solver_time_array', 'var') || isempty(solver_time_array)
        solver_time_array = [];
    end

    % 诊断结构体 (记录每次 outer 迭代)
    diagnostics = struct();
    diagnostics.iterations = struct();
    diagnostics.iterations.status = cell(1, max_iter);
    diagnostics.iterations.optval = zeros(1, max_iter);
    diagnostics.iterations.solve_time = NaN(1, max_iter);
    diagnostics.iterations.solver_name = cell(1, max_iter);
    diagnostics.iterations.hpipm_iters = zeros(1, max_iter);
    diagnostics.iterations.nq = zeros(1, max_iter);
    diagnostics.iterations.max_viol = NaN(1, max_iter);
    diagnostics.iterations.u_diff_to_prev = NaN(1, max_iter);
    diagnostics.iterations.qk_wheel_viol = NaN(1, max_iter);
    diagnostics.iterations.qk_cone_viol = NaN(1, max_iter);
    diagnostics.iterations.orig_wheel_viol = NaN(1, max_iter);
    diagnostics.iterations.orig_cone_viol = NaN(1, max_iter);
    diagnostics.step = step;
    diagnostics.max_iter = max_iter;
    diagnostics.step_failed = false;
    diagnostics.solver_call_count = 0;

    % ================= Python 环境路径设置 =================
    persistent py_path_added py_reloaded script_path;
    if isempty(py_path_added)
        script_path = fileparts(mfilename('fullpath'));
        if exist(script_path, 'dir')
            sys_mod = py.importlib.import_module('sys');
            py.getattr(sys_mod, 'path').insert(0, script_path);
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

    total_solve_time = 0;

    % ================= Outer 循环 (严格 3 次, 无 inner loop) =================
    for outer = 1:max_iter
        u_anchor = u_hat;  % RSS 凸化锚点

        % 诊断默认值 (失败时保留)
        outer_status = 'Failed';
        outer_optval = NaN;
        outer_solve_time = NaN;
        outer_solver_name = 'HPIPM-OCPQCQP';
        outer_hpipm_iters = 0;
        outer_nq = 0;
        outer_max_viol = NaN;
        outer_u_diff_val = NaN;
        solver_calls_this = 0;

        try
            % ========== 论文 Alg.1 line 4: 构造精确凸子问题 Q_K(u_anchor) ==========
            % OCP QCQP (精确凸二次约束, HPIPM ocp_qcqp)
            ocp = construct_ocp_qcqp_from_rss(path, step, v0, state, u_anchor, params);
            outer_nq = sum(ocp.nq_per_stage);

            % ========== 数据格式转换: cell → 2D numpy ==========
            b_stack = zeros(6, ocp.K);
            r_stack = zeros(3, ocp.K);
            for n = 1:ocp.K
                b_stack(:, n) = ocp.b{n};
                r_stack(:, n) = ocp.r{n};
            end

            % ========== 论文 Alg.1 line 5: 求解 u^(m+1) = S(u^(m)) ==========
            hpipm_mod = py.importlib.import_module('hpipm_qp_solver');
            result = hpipm_mod.solve_ocp_qcqp(...
                py.numpy.array(ocp.A), ...
                py.numpy.array(ocp.B), ...
                py.numpy.array(b_stack), ...
                py.numpy.array(ocp.Q_eff), ...
                py.numpy.array(ocp.S_eff), ...
                py.numpy.array(ocp.R_eff), ...
                py.numpy.array(ocp.q_stack), ...
                py.numpy.array(r_stack), ...
                py.numpy.array(int32(ocp.nx)), ...
                py.numpy.array(int32(ocp.nu)), ...
                py.numpy.array(int32(ocp.nq)), ...
                py.numpy.array(int32(ocp.nbx)), ...
                py.numpy.array(int32(ocp.nq_per_stage)), ...
                py.numpy.array(ocp.Qq_stack), ...
                py.numpy.array(ocp.Sq_stack), ...
                py.numpy.array(ocp.Rq_stack), ...
                py.numpy.array(ocp.qq_stack), ...
                py.numpy.array(ocp.rq_stack), ...
                py.numpy.array(ocp.uq_stack), ...
                py.numpy.array(ocp.x0), ...
                py.numpy.array(int32(ocp.idxbx)), ...
                py.float(ocp.const), ...
                py.bool(false) ...                    % verbose
            );

            % ========== 提取结果 ==========
            x = double(result{'x'});
            status_code = double(result{'status'});
            outer_optval = double(result{'obj_value_full'});
            outer_solve_time = double(result{'solve_time'});
            outer_hpipm_iters = double(result{'iters'});
            solver_calls_this = double(result{'solver_call_count'});
            outer_solver_name = 'HPIPM-OCPQCQP';

            % x_out = [u(:); nu(:)] (与 Dense QCQP 格式一致)
            u_sol = reshape(x(1:3*K), 3, K);

        catch ME
            u_sol = zeros(3, K);
            status_code = -1;
            outer_optval = NaN;
            outer_solve_time = NaN;
            outer_hpipm_iters = 0;
            solver_calls_this = 0;
            outer_solver_name = 'HPIPM-Error';
            fprintf('第%d步 outer=%d - Python 调用异常: %s\n', step, outer, ME.message);
            for si = 1:min(length(ME.stack), 5)
                fprintf('    at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
            end
        end

        % 累计真实 solver 调用次数
        diagnostics.solver_call_count = diagnostics.solver_call_count + solver_calls_this;
        total_solve_time = total_solve_time + outer_solve_time;
        if ~isempty(solver_time_array)
            solver_time_array(end+1) = outer_solve_time;
        else
            solver_time_array = outer_solve_time;
        end

        if status_code == 0
            outer_status = 'Solved';
        else
            outer_status = sprintf('Failed(%d)', status_code);
        end

        fprintf('第%d步第%d次迭代 - HPIPM-OCPQCQP内部求解时间：%.6f秒 (status=%d)\n', ...
            step, outer, outer_solve_time, status_code);
        fprintf('最优代价: %f | 求解状态: %s\n', outer_optval, outer_status);

        % ========== 严格接受条件: status==0 且解有限 ==========
        % status=1 (MAX_ITER) 不视为成功
        if status_code == 0 && all(isfinite(u_sol(:)))
            outer_u_diff = max(abs(u_sol(:) - u_hat(:)));
            outer_u_diff_val = outer_u_diff;

            % 精确约束检查 (固定 Q_K(u_anchor))
            [max_viol, wheel_viol, cone_viol] = check_qk_violation(u_sol, u_anchor, v0, params);
            outer_max_viol = max_viol;

            % 原始物理约束违反量 (双线性约束, 与 Q_K 凸化约束不同)
            [orig_wheel_viol, orig_cone_viol] = compute_original_violation(u_sol, v0, params);
            diagnostics.iterations.orig_wheel_viol(outer) = orig_wheel_viol;
            diagnostics.iterations.orig_cone_viol(outer) = orig_cone_viol;
            diagnostics.iterations.qk_wheel_viol(outer) = wheel_viol;
            diagnostics.iterations.qk_cone_viol(outer) = cone_viol;

            % 更新 u_hat (与 Dense QCQP 一致: status==0 即更新)
            u_hat = u_sol;
        end

        % 记录 outer 诊断
        diagnostics.iterations.status{outer} = outer_status;
        diagnostics.iterations.optval(outer) = outer_optval;
        diagnostics.iterations.solve_time(outer) = outer_solve_time;
        diagnostics.iterations.solver_name{outer} = outer_solver_name;
        diagnostics.iterations.hpipm_iters(outer) = outer_hpipm_iters;
        diagnostics.iterations.nq(outer) = outer_nq;
        if ~isnan(outer_max_viol)
            diagnostics.iterations.max_viol(outer) = outer_max_viol;
        end
        if ~isnan(outer_u_diff_val)
            diagnostics.iterations.u_diff_to_prev(outer) = outer_u_diff_val;
        end
    end  % end outer loop

    % ================= 输出: 仅当 3 次 outer 全部成功才标记 strict =================
    u = u_hat;
    all_outer_solved = true;
    for outer = 1:max_iter
        if ~strcmp(diagnostics.iterations.status{outer}, 'Solved')
            all_outer_solved = false;
            break;
        end
    end
    if all_outer_solved
        diagnostics.incumbent_type = 'strict';
        diagnostics.step_approximate = false;
        diagnostics.has_strict_incumbent = true;
        diagnostics.step_failed = false;
    else
        diagnostics.incumbent_type = 'none';
        diagnostics.step_approximate = false;
        diagnostics.has_strict_incumbent = false;
        diagnostics.step_failed = true;
    end

    % ================= 论文 Alg.1 line 11: 输出 ν_1 = ν_0 + u_1 =================
    new_state_dot = [cos(state(3)), -sin(state(3)), 0;
                     sin(state(3)),  cos(state(3)), 0;
                         0,              0, 1] * (state_dot + 1.00 * u(:, 1));
    velocity = v0 + u(:, 1);

    % ================= [DIAGNOSTIC] 原始约束违反量检查 (per-step) =================
    [orig_wheel_viol_final, orig_cone_viol_final] = compute_original_violation(u, v0, params);
    diagnostics.orig_wheel_viol_final = orig_wheel_viol_final;
    diagnostics.orig_cone_viol_final = orig_cone_viol_final;
    if mod(step, 25) == 0 || step == 1
        fprintf('[DIAG step=%d] wheel_viol_cur=%.6e, cone_viol_cur=%.6e\n', ...
            step, orig_wheel_viol_final, orig_cone_viol_final);
    end

    % 汇总诊断
    diagnostics.total_solve_time = total_solve_time;
end


function [wheel_viol, cone_viol] = compute_original_violation(u, v0, params)
% COMPUTE_ORIGINAL_VIOLATION 原始物理约束违反量 (双线性, 非 Q_K 凸化)
%   轮速: ||H_n * nu_k|| <= vimax
%   转向锥: a'*M*b >= 0 (a=nu_{k-1}, b=nu_k)
    K = 6;
    num_wheels = size(params.wheel_pos, 1);
    Hn = cell(1, num_wheels);
    for n = 1:num_wheels
        Hn{n} = [1, 0, -params.wheel_pos(n,2); 0, 1, params.wheel_pos(n,1)];
    end
    delta_th = params.dt * params.phidotmax;
    R1d = [sin(delta_th), -cos(delta_th); cos(delta_th), sin(delta_th)];
    R2d = R1d';
    nu = zeros(3, K+1); nu(:, 1) = v0;
    for kk = 1:K; nu(:, kk+1) = nu(:,kk) + u(:, kk); end
    % 轮速约束 (k=1..K)
    wheel_viol = 0;
    for kk = 1:K
        nkk = nu(:, kk+1);
        for n = 1:num_wheels
            vn = norm(Hn{n} * nkk, 2);
            viol = vn - params.vimax;
            if viol > wheel_viol; wheel_viol = viol; end
        end
    end
    % 转向锥 (k=1..K)
    cone_viol = 0;
    for kk = 1:K
        a = nu(:, kk);
        b = nu(:, kk+1);
        for n = 1:num_wheels
            Hn1 = Hn{n};
            for gg = 1:2
                if gg == 1; Rgg = R1d; else; Rgg = R2d; end
                Mgg = Hn1' * Rgg * Hn1;
                term = a' * Mgg * b;
                viol = -term;
                if viol > cone_viol; cone_viol = viol; end
            end
        end
    end
end
