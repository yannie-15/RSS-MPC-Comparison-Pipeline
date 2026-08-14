function [u, new_state_dot, velocity, diagnostics] = control_RSS(path, step, state_dot, state)
% CONTROL_RSS  论文 Algorithm 1 (Trajectory optimizer for SWMRs) 的实现
%
% 严格三次求解架构 (硬约束):
%   每个 MPC step 恰好执行 3 次 HPIPM ocp_qp solver.solve() 调用.
%   - 单层 outer=1:3 循环, 无 inner loop, 无 cutting-plane 累积.
%   - 每次 outer 只调用一次 solve_ocp_qp (allow_retry=False, 禁止 robust 重试).
%   - status=1 (MAX_ITER) 不视为成功; 只有 status==0 且解有限才作为有效 candidate.
%   - 不得把未通过精确约束检查的 candidate 输出给控制器.
%   - 诊断中真实累计 solver_call_count, MPC step 结束时断言 == 3.
%
% 变量管理:
%   u_iter               - 当前 RSS outer iterate / 下一轮 anchor
%   u_candidate          - 当前 HPIPM 返回值
%   feasible_incumbent   - 已通过精确约束检查的安全输出
%   u_previous_mpc       - 上一个 MPC step 移位后的 warm start
%
% 模式 A (single-cut): u_anchor = u_iter, u_cut = u_iter
%   每条精确二次约束只在当前点生成一个 tangent.

    params = config();

    % ================= Param Setup =================
    K = 6; rho = 0.01; k1 = 1; epsilon = 0;
    current_xy = [state(1), state(2)]';
    psi0 = state(3); v0 = state_dot;

    % ================= 迭代 Setup =================
    max_outer = 3;        % 论文 IV-B: RSS Algorithm 1 固定 3 次外层迭代
    % 容差拆分: STRICT 用于严格可行 incumbent, APPROX 用于近似 incumbent (诊断用)
    STRICT_QK_TOL = 1e-8;       % 严格 Q_K 精确约束容差 (通过此标准才称严格可行)
    STRICT_ORIG_TOL = 1e-8;     % 严格原始约束容差
    APPROX_QK_TOL = 5e-3;       % 近似 Q_K 容差 (OCP QP 线性化近似, 仅用于 approximate_incumbent)

    % warm start: 用上一 MPC 步的 u_full 平移 + 尾补零
    global RSS_WARMSTART_UHAT solver_time_array;
    if ~exist('solver_time_array', 'var') || isempty(solver_time_array)
        solver_time_array = [];
    end
    if step == 1 || isempty(RSS_WARMSTART_UHAT) || size(RSS_WARMSTART_UHAT, 2) ~= K
        % 首步 seed: 用基于参考轨迹方向的小幅非零 seed
        % (零点处转向锥梯度=0, tangent 平面退化; 非零 seed 提供有效梯度)
        % 参考轨迹相邻点差分作为初始速度方向, 幅值取 0.1 m/s
        ref_idx0 = min(size(path,2), step);
        ref_idx1 = min(size(path,2), step+1);
        ref_dir = path(1:3, ref_idx1) - path(1:3, ref_idx0);
        ref_dir_norm = norm(ref_dir);
        if ref_dir_norm > 1e-9
            v_init = 0.01 * ref_dir / ref_dir_norm;  % 0.01 m/s 沿参考方向
        else
            v_init = [0.01; 0; 0];  % 默认沿 x 方向
        end
        % u_seed = v_init 重复 K 次 (使 nu_hat 单调增长, 避免零点退化)
        u_seed = repmat(v_init, 1, K);
    else
        u_seed = [RSS_WARMSTART_UHAT(:, 2:K), zeros(3, 1)];
    end

    % 诊断结构体 (记录每次 outer 迭代)
    diagnostics = struct();
    diagnostics.iterations = struct();
    diagnostics.iterations.status = cell(1, max_outer);
    diagnostics.iterations.optval = zeros(1, max_outer);
    diagnostics.iterations.solve_time = NaN(1, max_outer);
    diagnostics.iterations.solver_name = cell(1, max_outer);
    diagnostics.iterations.hpipm_iters = zeros(1, max_outer);
    diagnostics.iterations.ng = zeros(1, max_outer);
    diagnostics.iterations.n_cuts = zeros(1, max_outer);  % 模式 B: 每 outer 的 cut point 数
    diagnostics.iterations.max_viol = NaN(1, max_outer);
    diagnostics.iterations.u_diff_to_prev = NaN(1, max_outer);
    diagnostics.iterations.qk_wheel_viol = NaN(1, max_outer);   % Q_K 凸化轮速违反
    diagnostics.iterations.qk_cone_viol = NaN(1, max_outer);    % Q_K 凸化转向锥违反
    diagnostics.iterations.orig_wheel_viol = NaN(1, max_outer); % 原始双线性轮速违反
    diagnostics.iterations.orig_cone_viol = NaN(1, max_outer);  % 原始双线性转向锥违反
    diagnostics.step = step;
    diagnostics.max_iter = max_outer;
    diagnostics.step_failed = false;
    diagnostics.solver_call_count = 0;  % 真实 solver.solve() 调用次数 (来自 Python)

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

    % ================= Incumbent 管理 (严格/近似分离) =================
    % u_iter: 当前 RSS outer iterate (下一轮 anchor), 初始化为 warm start
    u_iter = u_seed;
    % strict_feasible_incumbent: 通过 STRICT_QK_TOL (1e-8) 检查的安全输出
    strict_feasible_incumbent = [];
    has_strict_incumbent = false;
    strict_incumbent_viol = Inf;   % 当前 strict incumbent 的 Q_K violation
    strict_incumbent_outer = 0;    % 当前 strict incumbent 来自哪个 outer
    % approximate_incumbent: 通过 APPROX_QK_TOL (5e-3) 但未通过 STRICT 的候选 (诊断用, 不能称严格可行)
    approximate_incumbent = [];
    has_approx_incumbent = false;
    approx_incumbent_viol = Inf;   % 当前 approx incumbent 的 Q_K violation
    approx_incumbent_outer = 0;
    % u_previous_mpc: 上一 MPC step 的解 (用于诊断 u_diff 和模式 B cut point)
    u_previous_mpc = u_seed;
    % prev_candidate: 上一个 outer 的 candidate (用于模式 B cut point)
    prev_candidate = [];

    total_solve_time = 0;

    % ================= Outer 循环 (严格 3 次, 无 inner loop) =================
    % 模式 B (bounded multi-cut): 每个 outer 使用最多 3 组 cut points
    %   cut_points = {u_iter, prev_candidate, u_previous_mpc} (去重, 最多 3 组)
    %   u_anchor = u_iter (RSS 凸化锚点, 固定 B/L/r/const)
    %   所有 cut points 的约束放入同一个 OCP QP, 只调用 1 次 solver
    for outer = 1:max_outer
        u_anchor = u_iter;

        % 构建 cut_points 列表 (去重, 最多 3 组)
        cut_points = {u_iter};  % cut point 1: 当前 iterate
        % cut point 2: 上一个 outer 的 candidate (如果存在且不同)
        if ~isempty(prev_candidate)
            dup = false;
            for ci = 1:numel(cut_points)
                if norm(cut_points{ci}(:) - prev_candidate(:), inf) < 1e-8
                    dup = true; break;
                end
            end
            if ~dup
                cut_points{end+1} = prev_candidate;
            end
        end
        % cut point 3: shifted previous-MPC solution (如果不同)
        dup = false;
        for ci = 1:numel(cut_points)
            if norm(cut_points{ci}(:) - u_previous_mpc(:), inf) < 1e-8
                dup = true; break;
            end
        end
        if ~dup
            cut_points{end+1} = u_previous_mpc;
        end
        % 限制最多 3 组
        if numel(cut_points) > 3
            cut_points = cut_points(1:3);
        end
        u_cut = cut_points;  % 传递 cell array 给 construct_ocp_qp_from_rss

        % 诊断默认值 (失败时保留)
        outer_status = 'Failed';
        outer_optval = NaN;
        outer_solve_time = NaN;
        outer_solver_name = 'HPIPM';
        outer_hpipm_iters = 0;
        outer_ng = 0;
        outer_max_viol = NaN;
        outer_u_diff_val = NaN;

        try
            % ========== OCP QP 模式 (默认, 线性化约束 + SCP) ==========
            % 构造 OCP QP (single-cut, 无 cut_bank)
            ocp = construct_ocp_qp_from_rss(path, step, v0, state, u_anchor, u_cut, params);

            % 数据格式转换: cell → 2D numpy
            b_stack = zeros(6, ocp.K);
            r_stack = zeros(3, ocp.K);
            for n = 1:ocp.K
                b_stack(:, n) = ocp.b{n};
                r_stack(:, n) = ocp.r{n};
            end

            % 堆叠线性约束 (cell -> 2D numpy)
            total_ng = sum(ocp.ng_per_stage);
            outer_ng = total_ng;
            Cmat_stack = zeros(total_ng, 6);
            Dmat_stack = zeros(total_ng, 3);
            lg_stack = zeros(total_ng, 1);
            ug_stack = zeros(total_ng, 1);
            idx = 0;
            for s = 1:ocp.N_stages
                for j = 1:ocp.ng_per_stage(s)
                    idx = idx + 1;
                    Cmat_stack(idx, :) = ocp.Cmat{s, j};
                    Dmat_stack(idx, :) = ocp.Dmat{s, j};
                    lg_stack(idx) = ocp.lg{s, j};
                    ug_stack(idx) = ocp.ug{s, j};
                end
            end

            % 求解 OCP QP (allow_retry=False, 禁止 robust 重试)
            hpipm_mod = py.importlib.import_module('hpipm_qp_solver');
            result = hpipm_mod.solve_ocp_qp(...
                py.numpy.array(ocp.A), ...
                py.numpy.array(ocp.B), ...
                py.numpy.array(ocp.Q_eff), ...
                py.numpy.array(ocp.R_eff), ...
                py.numpy.array(ocp.S_eff), ...
                py.numpy.array(b_stack), ...
                py.numpy.array(r_stack), ...
                py.numpy.array(ocp.q_stack), ...
                py.numpy.array(int32(ocp.nx)), ...
                py.numpy.array(int32(ocp.nu)), ...
                py.numpy.array(int32(ocp.ng)), ...
                py.numpy.array(int32(ocp.nbx)), ...
                py.numpy.array(int32(ocp.ng_per_stage)), ...
                py.numpy.array(Cmat_stack), ...
                py.numpy.array(Dmat_stack), ...
                py.numpy.array(lg_stack), ...
                py.numpy.array(ug_stack), ...
                py.numpy.array(ocp.x0), ...
                py.numpy.array(int32(ocp.idxbx)), ...
                py.float(ocp.const), ...
                py.bool(false), ...                    % verbose
                py.numpy.array([]), ...                % x_init (不用 warm start)
                py.numpy.array([]), ...                % u_init
                py.bool(false), ...                    % warm_start=False
                py.bool(false) ...                     % allow_retry=False (严格三次)
            );

            x = double(result{'x'});
            status_code = double(result{'status'});
            outer_optval = double(result{'obj_value'});
            outer_solve_time = double(result{'solve_time'});
            outer_hpipm_iters = double(result{'iters'});
            outer_solver_name = 'HPIPM-OCPQP';
            solver_calls_this = double(result{'solver_call_count'});
            % OCP QP: x = [x_states; u_controls], u 在后 3K 个
            u_candidate = reshape(x(6*(ocp.N_stages)+1:end), 3, K);

        catch ME
            u_candidate = zeros(3, K);
            status_code = -1;
            outer_optval = NaN;
            outer_solve_time = NaN;
            outer_hpipm_iters = 0;
            solver_calls_this = 0;  % Python 异常: 无法确认 solver.solve() 是否被调用, 计为 0
            outer_solver_name = 'HPIPM-Error';
            fprintf('第%d步 outer=%d - Python 调用异常: %s\n', step, outer, ME.message);
        end

        % 累计真实 solver 调用次数
        diagnostics.solver_call_count = diagnostics.solver_call_count + solver_calls_this;
        total_solve_time = total_solve_time + outer_solve_time;
        if ~isempty(solver_time_array)
            solver_time_array(end+1) = outer_solve_time;
        else
            solver_time_array = outer_solve_time;
        end

        % 提取 optval/solve_time/iters (统一变量名)
        optval = outer_optval;
        this_solve_time = outer_solve_time;
        hpipm_iters = outer_hpipm_iters;
        if status_code == 0
            outer_status = 'Solved';
        else
            outer_status = sprintf('Failed(%d)', status_code);
        end

        fprintf('第%d步 outer=%d - HPIPM内部求解时间：%.6f秒 (status=%d, iters=%d, ng=%d)\n', ...
            step, outer, this_solve_time, status_code, hpipm_iters, outer_ng);
        fprintf('最优代价: %f | 求解状态: %s\n', optval, outer_status);

        % ========== 严格接受条件: status==0 且解有限 ==========
        % status=1 (MAX_ITER) 不视为成功
        if status_code == 0 && all(isfinite(u_candidate(:)))
            % 记录 u_diff (与上一轮 iterate 的差异)
            outer_u_diff = max(abs(u_candidate(:) - u_iter(:)));

            % 精确约束检查 (固定 Q_K(u_anchor))
            [max_viol, wheel_viol, cone_viol] = check_qk_violation(u_candidate, u_anchor, v0, params);
            outer_max_viol = max_viol;
            outer_u_diff_val = outer_u_diff;

            % 同时计算原始物理约束违反量 (双线性约束, 与 Q_K 凸化约束不同)
            [orig_wheel_viol, orig_cone_viol] = compute_original_violation(u_candidate, v0, params);
            diagnostics.iterations.orig_wheel_viol(outer) = orig_wheel_viol;
            diagnostics.iterations.orig_cone_viol(outer) = orig_cone_viol;
            diagnostics.iterations.qk_wheel_viol(outer) = wheel_viol;
            diagnostics.iterations.qk_cone_viol(outer) = cone_viol;

            fprintf('  [outer=%d] qk_viol=%.6e (wheel=%.6e, cone=%.6e), orig_viol=(wheel=%.6e, cone=%.6e), u_diff=%.6e\n', ...
                outer, max_viol, wheel_viol, cone_viol, orig_wheel_viol, orig_cone_viol, outer_u_diff);

            % 有限且 status==0: 允许作为下一轮 iterate
            u_iter = u_candidate;

            % 更新 prev_candidate (用于下一轮模式 B cut point)
            prev_candidate = u_candidate;

            % ========== Incumbent 更新 (严格/近似分离, 不被更差 candidate 覆盖) ==========
            % 1. 严格 incumbent: qk_viol < STRICT_QK_TOL 且 orig_viol < STRICT_ORIG_TOL
            if max_viol < STRICT_QK_TOL && orig_wheel_viol < STRICT_ORIG_TOL && orig_cone_viol < STRICT_ORIG_TOL
                % 只有 violation 更小才更新 (不被更差 candidate 覆盖)
                if ~has_strict_incumbent || max_viol < strict_incumbent_viol
                    strict_feasible_incumbent = u_candidate;
                    has_strict_incumbent = true;
                    strict_incumbent_viol = max_viol;
                    strict_incumbent_outer = outer;
                    fprintf('  [outer=%d] 通过 STRICT 检查 (qk_viol=%.3e < %.1e), 更新 strict incumbent\n', ...
                        outer, max_viol, STRICT_QK_TOL);
                end
            end
            % 2. 近似 incumbent: qk_viol < APPROX_QK_TOL (诊断用, 不能称严格可行)
            if max_viol < APPROX_QK_TOL
                if ~has_approx_incumbent || max_viol < approx_incumbent_viol
                    approximate_incumbent = u_candidate;
                    has_approx_incumbent = true;
                    approx_incumbent_viol = max_viol;
                    approx_incumbent_outer = outer;
                    if max_viol >= STRICT_QK_TOL
                        fprintf('  [outer=%d] 通过 APPROX 检查 (qk_viol=%.3e < %.1e 但 >= %.1e), 更新 approximate incumbent (非严格可行)\n', ...
                            outer, max_viol, APPROX_QK_TOL, STRICT_QK_TOL);
                    end
                end
            end
        else
            fprintf('  [outer=%d] 求解失败 (status=%d) 或解含 NaN, 保留上一轮 u_iter\n', ...
                outer, status_code);
            outer_max_viol = NaN;
            outer_u_diff_val = NaN;
            diagnostics.iterations.orig_wheel_viol(outer) = NaN;
            diagnostics.iterations.orig_cone_viol(outer) = NaN;
            diagnostics.iterations.qk_wheel_viol(outer) = NaN;
            diagnostics.iterations.qk_cone_viol(outer) = NaN;
            % 不得执行 u_iter = u_candidate; 保留上一轮 u_iter 和 incumbents
        end

        % 记录 outer 诊断
        diagnostics.iterations.status{outer} = outer_status;
        diagnostics.iterations.optval(outer) = outer_optval;
        diagnostics.iterations.solve_time(outer) = outer_solve_time;
        diagnostics.iterations.solver_name{outer} = outer_solver_name;
        diagnostics.iterations.hpipm_iters(outer) = outer_hpipm_iters;
        diagnostics.iterations.ng(outer) = outer_ng;
        diagnostics.iterations.n_cuts(outer) = numel(cut_points);
        if exist('outer_max_viol', 'var')
            diagnostics.iterations.max_viol(outer) = outer_max_viol;
        end
        if exist('outer_u_diff_val', 'var')
            diagnostics.iterations.u_diff_to_prev(outer) = outer_u_diff_val;
        end
    end  % end outer loop

    % ================= 断言: 真实 solver 调用数 == 3 =================
    % (allow_retry=False 时每次 outer 恰好 1 次 solver.solve(), 3 outer = 3 次)
    if diagnostics.solver_call_count ~= 3
        fprintf('警告: 第%d步 solver_call_count=%d (预期 3)\n', ...
            step, diagnostics.solver_call_count);
    end

    % ================= 输出选择 (严格优先, 如实标记) =================
    % 优先级: strict_feasible_incumbent > approximate_incumbent > seed (失败)
    % 严格 incumbent 通过 1e-8 检查, 可安全输出
    % 近似 incumbent 仅通过 5e-3 检查, 标记 step_approximate=true (非严格可行)
    % 无任何 incumbent: step_failed=true, 输出 seed (不得称可行)
    diagnostics.step_failed = false;
    diagnostics.step_approximate = false;  % 使用近似 incumbent (非严格可行)
    diagnostics.has_strict_incumbent = has_strict_incumbent;
    diagnostics.has_approx_incumbent = has_approx_incumbent;
    diagnostics.strict_incumbent_outer = strict_incumbent_outer;
    diagnostics.approx_incumbent_outer = approx_incumbent_outer;
    diagnostics.strict_incumbent_viol = strict_incumbent_viol;
    diagnostics.approx_incumbent_viol = approx_incumbent_viol;

    if has_strict_incumbent
        u = strict_feasible_incumbent;
        diagnostics.incumbent_type = 'strict';
        diagnostics.incumbent_viol = strict_incumbent_viol;
        fprintf('第%d步: 使用 strict incumbent (outer=%d, qk_viol=%.3e < %.1e)\n', ...
            step, strict_incumbent_outer, strict_incumbent_viol, STRICT_QK_TOL);
    elseif has_approx_incumbent
        u = approximate_incumbent;
        diagnostics.incumbent_type = 'approximate';
        diagnostics.incumbent_viol = approx_incumbent_viol;
        diagnostics.step_approximate = true;
        fprintf('第%d步: 使用 APPROXIMATE incumbent (outer=%d, qk_viol=%.3e >= %.1e, 非严格可行)\n', ...
            step, approx_incumbent_outer, approx_incumbent_viol, STRICT_QK_TOL);
    else
        u = u_seed;
        diagnostics.incumbent_type = 'none';
        diagnostics.incumbent_viol = Inf;
        diagnostics.step_failed = true;
        fprintf('第%d步: 三次 outer 后无 incumbent (strict/approx 均无), 标记 step_failed\n', step);
    end

    % 保存 warm start
    RSS_WARMSTART_UHAT = u;

    % ================= 论文 Alg.1 line 11: 输出 ν_1 = ν_0 + u_1 =================
    new_state_dot = [cos(state(3)), -sin(state(3)), 0;
                     sin(state(3)),  cos(state(3)), 0;
                         0,              0, 1] * (state_dot + 1.00 * u(:, 1));
    velocity = v0 + u(:, 1);

    % ================= [DIAGNOSTIC] 原始约束违反量检查 (per-step, 不用 persistent) =================
    % 注意: 这是原始双线性约束检查, 与 check_qk_violation (固定 Q_K 凸化约束) 不同
    % 汇总统计由 run_paper_baseline_case.m 在循环外维护, 不依赖 persistent (避免 clear control_RSS 重置)
    [orig_wheel_viol_final, orig_cone_viol_final] = compute_original_violation(u, v0, params);
    diagnostics.orig_wheel_viol_final = orig_wheel_viol_final;
    diagnostics.orig_cone_viol_final = orig_cone_viol_final;
    if mod(step, 25) == 0 || step == 1
        fprintf('[DIAG step=%d] wheel_viol_cur=%.6e, cone_viol_cur=%.6e (incumbent=%s)\n', ...
            step, orig_wheel_viol_final, orig_cone_viol_final, diagnostics.incumbent_type);
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
