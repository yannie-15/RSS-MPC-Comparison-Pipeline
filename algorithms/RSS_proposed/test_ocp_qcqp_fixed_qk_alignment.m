function test_ocp_qcqp_fixed_qk_alignment()
% TEST_OCP_QCQP_FIXED_QK_ALIGNMENT
% Test C: 固定 Q_K 求解对齐 (Dense QCQP vs OCP QCQP)
%
% 阶段 1: step 1 / outer 1 对齐 (确保 Dense 和 OCP 使用相同输入)
% 阶段 2: 三次 outer 对齐 (逐步更新 u_anchor)
%
% 验收:
%   max |u_ocpqcqp - u_dense|       < 1e-6
%   max |v_ocpqcqp - v_dense|       < 1e-6
%   max dynamics residual            < 1e-10
%   max exact Q_K violation          < 1e-8
%   normalized/full objective gap    < 1e-7
%   status_dense = 0, status_ocp = 0
%
% 用法:
%   addpath('core'); setup_paths;
%   addpath('algorithms/RSS_proposed');
%   test_ocp_qcqp_fixed_qk_alignment

    addpath(fileparts(mfilename('fullpath')));
    params = config();

    K = 6; nx = 6; nu = 3;

    % 生成参考轨迹
    path = generateReference(params, params.num_path_pts);

    % 共同输入 (与 control_RSS_denseqcqp.m 一致)
    v0 = [0.01; 0.01; 0.01];
    state = [0.05, 0.1, 0.2]';
    step = 1;

    fprintf('============== Test C: OCP QCQP vs Dense QCQP 求解对齐 ==============\n');

    % ========== Python 环境设置 ==========
    sys_mod = py.importlib.import_module('sys');
    script_path = fileparts(mfilename('fullpath'));
    py.getattr(sys_mod, 'path').insert(0, script_path);
    try
        py.importlib.reload(py.importlib.import_module('hpipm_qp_solver'));
    catch
    end
    hpipm_mod = py.importlib.import_module('hpipm_qp_solver');

    %% ===== Phase 1: step 1 / outer 1 对齐 =====
    fprintf('\n--- Phase 1: step 1 / outer 1 ---\n');
    u_anchor = zeros(nu, K);  % u_hat = 0 (与 Dense QCQP 初始化一致)

    % --- Dense QCQP ---
    qp = construct_complete_qp_from_rss(path, step, v0, state, u_anchor, params);
    n_qcqp = length(qp.Hq);
    Hq_stack = zeros(qp.n_var, qp.n_var, n_qcqp);
    gq_stack = zeros(qp.n_var, n_qcqp);
    uq_stack = zeros(n_qcqp, 1);
    for i = 1:n_qcqp
        Hq_stack(:, :, i) = qp.Hq{i};
        gq_stack(:, i) = qp.gq{i};
        uq_stack(i) = qp.uq(i);
    end

    result_dense = hpipm_mod.solve_qcqp(...
        py.numpy.array(qp.H), py.numpy.array(qp.g), ...
        py.numpy.array(qp.A), py.numpy.array(qp.b), ...
        py.numpy.array(Hq_stack), py.numpy.array(gq_stack), ...
        py.numpy.array(uq_stack), py.bool(false));

    status_dense = double(result_dense{'status'});
    iters_dense = double(result_dense{'iters'});
    obj_dense_reduced = double(result_dense{'obj_value'});
    x_dense = double(result_dense{'x'}); x_dense = x_dense(:);  % 强制列向量 (36,1)
    u_dense = reshape(x_dense(1:nu*K), nu, K);
    nu_dense = reshape(x_dense(nu*K+1:end), nu, K);  % ν_1..ν_K

    % --- OCP QCQP ---
    ocp = construct_ocp_qcqp_from_rss(path, step, v0, state, u_anchor, params);
    b_stack = zeros(nx, ocp.K);
    r_stack = zeros(nu, ocp.K);
    for n = 1:ocp.K
        b_stack(:, n) = ocp.b{n};
        r_stack(:, n) = ocp.r{n};
    end

    result_ocp = hpipm_mod.solve_ocp_qcqp(...
        py.numpy.array(ocp.A), py.numpy.array(ocp.B), ...
        py.numpy.array(b_stack), ...
        py.numpy.array(ocp.Q_eff), py.numpy.array(ocp.S_eff), ...
        py.numpy.array(ocp.R_eff), py.numpy.array(ocp.q_stack), ...
        py.numpy.array(r_stack), ...
        py.numpy.array(int32(ocp.nx)), py.numpy.array(int32(ocp.nu)), ...
        py.numpy.array(int32(ocp.nq)), py.numpy.array(int32(ocp.nbx)), ...
        py.numpy.array(int32(ocp.nq_per_stage)), ...
        py.numpy.array(ocp.Qq_stack), py.numpy.array(ocp.Sq_stack), ...
        py.numpy.array(ocp.Rq_stack), py.numpy.array(ocp.qq_stack), ...
        py.numpy.array(ocp.rq_stack), py.numpy.array(ocp.uq_stack), ...
        py.numpy.array(ocp.x0), py.numpy.array(int32(ocp.idxbx)), ...
        py.float(ocp.const), py.bool(false));

    status_ocp = double(result_ocp{'status'});
    iters_ocp = double(result_ocp{'iters'});
    obj_ocp_full = double(result_ocp{'obj_value_full'});
    obj_ocp_reduced = double(result_ocp{'obj_value_reduced'});
    x_ocp = double(result_ocp{'x'}); x_ocp = x_ocp(:);  % 强制列向量 (36,1)
    u_ocp = reshape(x_ocp(1:nu*K), nu, K);
    nu_ocp = reshape(x_ocp(nu*K+1:end), nu, K);  % ν_1..ν_K

    % --- 比较 ---
    u_diff = max(abs(u_ocp(:) - u_dense(:)));
    nu_diff = max(abs(nu_ocp(:) - nu_dense(:)));
    obj_dense_full = obj_dense_reduced + qp.objective_constant;
    obj_full_gap = abs(obj_ocp_full - obj_dense_full);
    obj_reduced_gap = abs(obj_ocp_reduced - obj_dense_reduced);

    % 约束检查
    [dense_viol, dense_wheel, dense_cone] = check_qk_violation(u_dense, u_anchor, v0, params);
    [ocp_viol, ocp_wheel, ocp_cone] = check_qk_violation(u_ocp, u_anchor, v0, params);

    % 动力学残差 (Dense)
    dyn_res_dense = norm(qp.A * x_dense - qp.b, inf);

    fprintf('Dense: status=%d, iters=%d, obj_reduced=%.10f, obj_full=%.10f\n', ...
        status_dense, iters_dense, obj_dense_reduced, obj_dense_full);
    fprintf('OCP:   status=%d, iters=%d, obj_reduced=%.10f, obj_full=%.10f\n', ...
        status_ocp, iters_ocp, obj_ocp_reduced, obj_ocp_full);
    fprintf('\n');
    fprintf('max |u_ocp - u_dense|         = %.6e  (target < 1e-6)\n', u_diff);
    fprintf('max |v_ocp - v_dense|         = %.6e  (target < 1e-6)\n', nu_diff);
    fprintf('obj_full gap                  = %.6e  (target < 1e-7)\n', obj_full_gap);
    fprintf('obj_reduced gap               = %.6e  (target < 1e-7)\n', obj_reduced_gap);
    fprintf('Dense dynamics residual       = %.6e  (target < 1e-10)\n', dyn_res_dense);
    fprintf('Dense Q_K violation           = %.6e  (wheel=%.2e, cone=%.2e)\n', ...
        dense_viol, dense_wheel, dense_cone);
    fprintf('OCP Q_K violation             = %.6e  (wheel=%.2e, cone=%.2e)\n', ...
        ocp_viol, ocp_wheel, ocp_cone);

    phase1_pass = status_dense == 0 && status_ocp == 0 && ...
                  u_diff < 1e-6 && nu_diff < 1e-6 && ...
                  obj_full_gap < 1e-7 && dyn_res_dense < 1e-10 && ...
                  dense_viol < 1e-8 && ocp_viol < 1e-8;
    if phase1_pass
        fprintf('Phase 1 结论: PASS\n');
    else
        fprintf('Phase 1 结论: FAIL\n');
    end

    %% ===== Phase 2: 三次 outer 对齐 =====
    fprintf('\n--- Phase 2: 三次 outer ---\n');
    u_hat_dense = zeros(nu, K);  % Dense 起点
    u_hat_ocp = zeros(nu, K);    % OCP 起点
    solver_calls_dense = 0;
    solver_calls_ocp = 0;
    obj_dense_outer = NaN;  % 预初始化
    obj_ocp_outer = NaN;    % 预初始化

    phase2_pass = true;
    for outer = 1:3
        fprintf('\n[outer=%d]\n', outer);

        % --- Dense QCQP ---
        u_anchor_d = u_hat_dense;
        qp_d = construct_complete_qp_from_rss(path, step, v0, state, u_anchor_d, params);
        Hq_d = zeros(qp_d.n_var, qp_d.n_var, length(qp_d.Hq));
        gq_d = zeros(qp_d.n_var, length(qp_d.Hq));
        uq_d = zeros(length(qp_d.Hq), 1);
        for i = 1:length(qp_d.Hq)
            Hq_d(:, :, i) = qp_d.Hq{i};
            gq_d(:, i) = qp_d.gq{i};
            uq_d(i) = qp_d.uq(i);
        end
        rd = hpipm_mod.solve_qcqp(...
            py.numpy.array(qp_d.H), py.numpy.array(qp_d.g), ...
            py.numpy.array(qp_d.A), py.numpy.array(qp_d.b), ...
            py.numpy.array(Hq_d), py.numpy.array(gq_d), ...
            py.numpy.array(uq_d), py.bool(false));
        sd = double(rd{'status'});
        od = double(rd{'obj_value'});
        xd = double(rd{'x'}); xd = xd(:);  % 强制列向量
        solver_calls_dense = solver_calls_dense + 1;
        if sd == 0 && all(isfinite(xd))
            u_hat_dense = reshape(xd(1:nu*K), nu, K);
            obj_dense_outer = od + qp_d.objective_constant;
        end
        [dv, ~, ~] = check_qk_violation(u_hat_dense, u_anchor_d, v0, params);

        % --- OCP QCQP ---
        u_anchor_o = u_hat_ocp;
        ocp_o = construct_ocp_qcqp_from_rss(path, step, v0, state, u_anchor_o, params);
        b_o = zeros(nx, ocp_o.K);
        r_o = zeros(nu, ocp_o.K);
        for n = 1:ocp_o.K
            b_o(:, n) = ocp_o.b{n};
            r_o(:, n) = ocp_o.r{n};
        end
        ro = hpipm_mod.solve_ocp_qcqp(...
            py.numpy.array(ocp_o.A), py.numpy.array(ocp_o.B), ...
            py.numpy.array(b_o), ...
            py.numpy.array(ocp_o.Q_eff), py.numpy.array(ocp_o.S_eff), ...
            py.numpy.array(ocp_o.R_eff), py.numpy.array(ocp_o.q_stack), ...
            py.numpy.array(r_o), ...
            py.numpy.array(int32(ocp_o.nx)), py.numpy.array(int32(ocp_o.nu)), ...
            py.numpy.array(int32(ocp_o.nq)), py.numpy.array(int32(ocp_o.nbx)), ...
            py.numpy.array(int32(ocp_o.nq_per_stage)), ...
            py.numpy.array(ocp_o.Qq_stack), py.numpy.array(ocp_o.Sq_stack), ...
            py.numpy.array(ocp_o.Rq_stack), py.numpy.array(ocp_o.qq_stack), ...
            py.numpy.array(ocp_o.rq_stack), py.numpy.array(ocp_o.uq_stack), ...
            py.numpy.array(ocp_o.x0), py.numpy.array(int32(ocp_o.idxbx)), ...
            py.float(ocp_o.const), py.bool(false));
        so = double(ro{'status'});
        oo_full = double(ro{'obj_value_full'});
        xo = double(ro{'x'}); xo = xo(:);  % 强制列向量
        solver_calls_ocp = solver_calls_ocp + 1;
        if so == 0 && all(isfinite(xo))
            u_hat_ocp = reshape(xo(1:nu*K), nu, K);
            obj_ocp_outer = oo_full;
        end
        [ov, ~, ~] = check_qk_violation(u_hat_ocp, u_anchor_o, v0, params);

        % 比较
        u_diff_o = max(abs(u_hat_ocp(:) - u_hat_dense(:)));
        obj_gap_o = abs(obj_ocp_outer - obj_dense_outer);
        fprintf('  Dense: status=%d, obj_full=%.10f, qk_viol=%.6e\n', sd, obj_dense_outer, dv);
        fprintf('  OCP:   status=%d, obj_full=%.10f, qk_viol=%.6e\n', so, obj_ocp_outer, ov);
        fprintf('  u_diff=%.6e, obj_gap=%.6e\n', u_diff_o, obj_gap_o);

        if sd ~= 0 || so ~= 0 || u_diff_o >= 1e-6 || obj_gap_o >= 1e-7 || isnan(obj_gap_o)
            phase2_pass = false;
        end
    end

    fprintf('\n--- Phase 2 Summary ---\n');
    fprintf('solver_calls: Dense=%d, OCP=%d (expected 3 each)\n', solver_calls_dense, solver_calls_ocp);
    fprintf('final u_diff = %.6e  (target < 1e-6)\n', max(abs(u_hat_ocp(:) - u_hat_dense(:))));
    if phase2_pass && solver_calls_dense == 3 && solver_calls_ocp == 3
        fprintf('Phase 2 结论: PASS\n');
    else
        fprintf('Phase 2 结论: FAIL\n');
    end

    %% ===== Overall =====
    fprintf('\n========== Overall ==========\n');
    if phase1_pass && phase2_pass
        fprintf('结论: PASS (OCP QCQP 与 Dense QCQP 在固定 Q_K 下对齐)\n');
    else
        fprintf('结论: FAIL\n');
    end
    fprintf('==============================================================\n');
end
