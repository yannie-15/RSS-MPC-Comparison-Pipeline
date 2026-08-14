function test_ocp_qcqp_construction_equivalence()
% TEST_OCP_QCQP_CONSTRUCTION_EQUIVALENCE
% Test A: 构造级逐点恒等测试 (Dense QCQP vs OCP QCQP)
%
% 不调用求解器, 验证两个构造在数学上等价:
%   1. 动力学等式: Dense A*x=b vs OCP x_{n+1}=A*x_n+B*u_n+b_n
%   2. 完整目标: obj_dense (含常数) vs obj_ocp (含常数)
%   3. 24 条 wheel 约束: 逐条 dense_value == ocp_value
%   4. 48 条 steering 约束: 逐条 dense_value == 2 * ocp_value
%      (Dense 代码将 A/B/L 各项整体 ×2, 但约束 2*C<=0 与 C<=0 等价)
%
% 验收: max error < 1e-10
%
% 用法:
%   addpath('core'); setup_paths;
%   addpath('algorithms/RSS_proposed');
%   addpath('algorithms/RSS_proposed/tests');
%   addpath('paper_reproduction');
%   test_ocp_qcqp_construction_equivalence

    % 测试文件位于 algorithms/RSS_proposed/tests/, 需添加父目录以访问 control/config 等
    addpath(fileparts(fileparts(mfilename('fullpath'))));
    params = config();

    K = 6;
    N_stages = K + 1;
    nx = 6;
    nu = 3;
    num_wheels = size(params.wheel_pos, 1);

    % 生成参考轨迹
    path = generateReference(params, params.num_path_pts);

    fprintf('============== Test A: OCP QCQP vs Dense QCQP 构造恒等 ==============\n');

    % 测试多个随机 u_anchor 和测试状态
    n_test_points = 5;
    max_dyn_dense_err = 0;
    max_dyn_ocp_err = 0;
    max_obj_err = 0;
    max_wheel_err = 0;
    max_steer_err = 0;
    n_constraints_checked = 0;

    rng(42);  % 可复现
    for tp = 1:n_test_points
        % 随机 u_anchor 和测试状态
        u_anchor = 0.01 * randn(nu, K);
        v0 = 0.01 * randn(3, 1);  % 3×1 当前车体系速度
        state = [0.05 * randn(2, 1); 0.2 * randn];  % 3×1 [x, y, psi]
        step = 1 + randi(20);  % 随机 step (1..21)

        % 随机 u_test (满足动力学的测试点)
        u_test = 0.02 * randn(nu, K);
        nu_test = zeros(3, K+1);
        nu_test(:, 1) = v0;
        for k = 1:K
            nu_test(:, k+1) = nu_test(:, k) + u_test(:, k);
        end

        % 构造 Dense QCQP 和 OCP QCQP
        qp_dense = construct_complete_qp_from_rss(path, step, v0, state, u_anchor, params);
        ocp = construct_ocp_qcqp_from_rss(path, step, v0, state, u_anchor, params);

        % ===== 1. 动力学一致性 =====
        % Dense: x_dense = [u_test(:); nu_test(:, 2:K+1)], 验证 A_eq * x_dense ≈ b_eq
        x_dense_test = [u_test(:); reshape(nu_test(:, 2:K+1), [], 1)];
        dyn_res_dense = norm(qp_dense.A * x_dense_test - qp_dense.b, inf);

        % OCP: 构造 x_ocp 序列, 验证 x_{n+1} = A*x_n + B*u_n + b_n
        x_ocp = zeros(nx, N_stages);
        ref_idx_0 = min(size(path, 2), step);
        e0 = [state(1) - path(1, ref_idx_0); ...
              state(2) - path(2, ref_idx_0); ...
              state(3) - path(3, ref_idx_0)];
        x_ocp(:, 1) = [e0; v0];
        dyn_res_ocp = 0;
        for n = 1:K
            x_next = ocp.A * x_ocp(:, n) + ocp.B * u_test(:, n) + ocp.b{n};
            x_ocp(:, n+1) = x_next;
            % 检查 v 分量与 nu_test 一致 (OCP v_n = paper ν_n)
            v_diff = norm(x_next(4:6) - nu_test(:, n+1), inf);
            dyn_res_ocp = max(dyn_res_ocp, v_diff);
        end

        max_dyn_dense_err = max(max_dyn_dense_err, dyn_res_dense);
        max_dyn_ocp_err = max(max_dyn_ocp_err, dyn_res_ocp);

        % ===== 2. 目标函数一致性 =====
        % Dense: 0.5*x'H*x + g'*x + objective_constant
        obj_dense = 0.5 * x_dense_test' * qp_dense.H * x_dense_test ...
                    + qp_dense.g' * x_dense_test + qp_dense.objective_constant;

        % OCP: const + Σ_n [0.5*xn'Q_n*xn] + Σ_{n<K} [0.5*un'R*un + r_n'un]
        %     (S=0, q=0 by construction)
        obj_ocp = ocp.const;
        for n = 1:N_stages
            xn = x_ocp(:, n);
            Q_n = ocp.Q_eff(:, (n-1)*nx+1:n*nx);
            obj_ocp = obj_ocp + 0.5 * xn' * Q_n * xn;
            if n <= K
                un = u_test(:, n);
                obj_ocp = obj_ocp + 0.5 * un' * ocp.R_eff * un + ocp.r{n}' * un;
            end
        end

        obj_err = abs(obj_dense - obj_ocp);
        max_obj_err = max(max_obj_err, obj_err);

        % ===== 3 & 4. 72 条二次约束一致性 =====
        % Dense constraint value: 0.5*x'Hq*x + gq'*x - uq
        % OCP constraint value:    0.5*x'Qq*x + x'Sq'*u + 0.5*u'Rq*u + qq'*x + rq'*u - uq
        %
        % 约束缩放关系 (经数学推导验证):
        %   wheel:    dense_val = C_wheel,    ocp_val = C_wheel     → dense_val == ocp_val
        %   steering: dense_val = 2*C_steer,  ocp_val = C_steer     → dense_val == 2*ocp_val
        %   (Dense 代码对 steering 的 A/B/L 各项整体 ×2, 但 2*C<=0 ⟺ C<=0)

        % --- 构造 Dense 元数据 (idx → kind, k, wheel, rotation) ---
        % Dense 顺序:
        %   idx 1..24:  wheel, k=1..6, n=1..4   (dense_idx = (k-1)*4 + n)
        %   idx 25..48: steering R1, k=1..6, n=1..4  (dense_idx = 24 + (k-1)*4 + n)
        %   idx 49..72: steering R2, k=1..6, n=1..4  (dense_idx = 48 + (k-1)*4 + n)
        dense_kind = cell(72, 1);
        dense_k = zeros(72, 1);
        dense_wheel = zeros(72, 1);
        dense_rotation = cell(72, 1);
        idx = 0;
        for k = 1:K
            for n = 1:num_wheels
                idx = idx + 1;
                dense_kind{idx} = 'wheel';
                dense_k(idx) = k;
                dense_wheel(idx) = n;
                dense_rotation{idx} = '';
            end
        end
        for k = 1:K
            for n = 1:num_wheels
                idx = idx + 1;
                dense_kind{idx} = 'steering';
                dense_k(idx) = k;
                dense_wheel(idx) = n;
                dense_rotation{idx} = 'R1';
            end
        end
        for k = 1:K
            for n = 1:num_wheels
                idx = idx + 1;
                dense_kind{idx} = 'steering';
                dense_k(idx) = k;
                dense_wheel(idx) = n;
                dense_rotation{idx} = 'R2';
            end
        end
        assert(idx == 72, 'Dense 元数据构造错误: idx=%d != 72', idx);

        % --- 构造 OCP idx 查找表 (key → ocp_idx) ---
        ocp_idx_map = containers.Map();
        for j = 1:72
            key = sprintf('%s_k%d_w%d_%s', ocp.metadata.kind{j}, ocp.metadata.k(j), ...
                ocp.metadata.wheel(j), ocp.metadata.rotation{j});
            assert(~ocp_idx_map.isKey(key), 'OCP 元数据重复: key=%s', key);
            ocp_idx_map(key) = j;
        end

        % --- 逐条比较 ---
        tp_wheel_err = 0;
        tp_steer_err = 0;
        for d_idx = 1:72
            key = sprintf('%s_k%d_w%d_%s', dense_kind{d_idx}, dense_k(d_idx), ...
                dense_wheel(d_idx), dense_rotation{d_idx});
            assert(ocp_idx_map.isKey(key), 'Dense idx %d (key=%s) 在 OCP 中找不到对应约束', d_idx, key);
            o_idx = ocp_idx_map(key);

            % --- Dense value ---
            Hq_d = qp_dense.Hq{d_idx};
            gq_d = qp_dense.gq{d_idx};
            uq_d = qp_dense.uq(d_idx);
            dense_val = 0.5 * x_dense_test' * Hq_d * x_dense_test + gq_d' * x_dense_test - uq_d;

            % --- OCP value ---
            Qq_o = ocp.Qq_stack(:, :, o_idx);
            Sq_o = ocp.Sq_stack(:, :, o_idx);
            Rq_o = ocp.Rq_stack(:, :, o_idx);
            qq_o = ocp.qq_stack(:, o_idx);
            rq_o = ocp.rq_stack(:, o_idx);
            uq_o = ocp.uq_stack(o_idx);

            % 找到该约束所在的 OCP stage (0-indexed)
            s = ocp.metadata.stage(o_idx);
            xn = x_ocp(:, s+1);  % MATLAB 1-indexed
            if s < K
                un = u_test(:, s+1);  % u_s = u_test(:, s+1) (OCP u_n = paper u_{n+1})
            else
                un = zeros(nu, 1);
            end

            ocp_val = 0.5 * xn' * Qq_o * xn + xn' * Sq_o' * un + 0.5 * un' * Rq_o * un ...
                      + qq_o' * xn + rq_o' * un - uq_o;

            % 比较 (wheel: dense==ocp; steering: dense==2*ocp)
            if strcmp(dense_kind{d_idx}, 'wheel')
                err = abs(dense_val - ocp_val);
                tp_wheel_err = max(tp_wheel_err, err);
            else
                err = abs(dense_val - 2 * ocp_val);
                tp_steer_err = max(tp_steer_err, err);
            end
            n_constraints_checked = n_constraints_checked + 1;
        end

        max_wheel_err = max(max_wheel_err, tp_wheel_err);
        max_steer_err = max(max_steer_err, tp_steer_err);

        fprintf('Test point %d (step=%d): dyn_d=%.3e, dyn_o=%.3e, obj=%.3e, wheel=%.3e, steer=%.3e\n', ...
            tp, step, dyn_res_dense, dyn_res_ocp, obj_err, tp_wheel_err, tp_steer_err);
    end

    fprintf('\n========== Summary ==========\n');
    fprintf('Test points:              %d\n', n_test_points);
    fprintf('Constraints checked:      %d (expected %d)\n', n_constraints_checked, n_test_points * 72);
    fprintf('Max dynamics residual (Dense): %.3e (target < 1e-10)\n', max_dyn_dense_err);
    fprintf('Max dynamics residual (OCP):   %.3e (target < 1e-10)\n', max_dyn_ocp_err);
    fprintf('Max objective error:           %.3e (target < 1e-10)\n', max_obj_err);
    fprintf('Max wheel constraint error:    %.3e (target < 1e-10)\n', max_wheel_err);
    fprintf('Max steering constraint error: %.3e (target < 1e-10)\n', max_steer_err);

    pass = max_dyn_dense_err < 1e-10 && max_dyn_ocp_err < 1e-10 && ...
           max_obj_err < 1e-10 && max_wheel_err < 1e-10 && max_steer_err < 1e-10;
    if pass
        fprintf('结论: PASS (构造级逐点恒等)\n');
    else
        fprintf('结论: FAIL\n');
        if max_dyn_dense_err >= 1e-10
            fprintf('  原因: Dense 动力学残差 = %.3e >= 1e-10\n', max_dyn_dense_err);
        end
        if max_dyn_ocp_err >= 1e-10
            fprintf('  原因: OCP 动力学残差 = %.3e >= 1e-10\n', max_dyn_ocp_err);
        end
        if max_obj_err >= 1e-10
            fprintf('  原因: 目标函数误差 = %.3e >= 1e-10\n', max_obj_err);
        end
        if max_wheel_err >= 1e-10
            fprintf('  原因: 轮速约束误差 = %.3e >= 1e-10\n', max_wheel_err);
        end
        if max_steer_err >= 1e-10
            fprintf('  原因: 转向锥约束误差 = %.3e >= 1e-10\n', max_steer_err);
        end
    end
    fprintf('==============================================================\n');
end
