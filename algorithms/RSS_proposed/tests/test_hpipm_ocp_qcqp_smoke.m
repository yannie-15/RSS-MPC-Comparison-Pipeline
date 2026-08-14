function test_hpipm_ocp_qcqp_smoke()
% TEST_HPIPM_OCP_QCQP_SMOKE
% Test B: HPIPM OCP QCQP wrapper smoke tests
%
% 逐步验证 HPIPM OCP QCQP wrapper 的正确性 (从最小问题到全 72 条约束):
%   1. 无 nq (纯 OCP QP, 退化情况)
%   2. 一条 terminal wheel constraint (nu=0 stage 的 Qq/qq/uq)
%   3. 一条 stage steering constraint (含 cross term Sq, 终端前 stage)
%   4. 全部 24 条 wheel constraints (RSS, 无 steering)
%   5. 全部 72 条 constraints (RSS 完整问题)
%
% 每个 test 验收: status==0, 解有限, 残差有限
%
% 用法:
%   addpath('core'); setup_paths;
%   addpath('algorithms/RSS_proposed');
%   addpath('algorithms/RSS_proposed/tests');
%   test_hpipm_ocp_qcqp_smoke

    % 测试文件位于 algorithms/RSS_proposed/tests/, 需添加父目录以访问 hpipm_qp_solver 等
    parent_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(parent_dir);

    % ========== Python 环境设置 ==========
    sys_mod = py.importlib.import_module('sys');
    py.getattr(sys_mod, 'path').insert(0, parent_dir);
    try
        py.importlib.reload(py.importlib.import_module('hpipm_qp_solver'));
    catch
    end
    hpipm_mod = py.importlib.import_module('hpipm_qp_solver');

    fprintf('============== Test B: HPIPM OCP QCQP wrapper smoke ==============\n');
    n_pass = 0;
    n_fail = 0;

    %% ===== Test 1: 无 nq (纯 OCP QP) =====
    fprintf('\n--- Test 1: 无 nq (纯 OCP QP) ---\n');
    try
        % K=2, nx=2, nu=1
        % 动力学: x_{n+1} = x_n + u_n
        % 代价: Σ x_n^2 + u_n^2
        K1 = 2; nx1 = 2; nu1 = 1; N1 = K1 + 1;
        A1 = eye(nx1);
        B1 = [1; 0];
        b_stack1 = zeros(nx1, K1);
        Q_stack1 = 2 * repmat(eye(nx1), 1, N1);  % (2, 2*3) = (2,6), 逐 stage Q=2*I
        R_stack1 = 2 * repmat(eye(nu1), 1, K1);  % (1, 1*2) = (1,2)
        S_stack1 = zeros(nu1, nx1 * K1);          % (1, 2*2) = (1,4)
        q_stack1 = zeros(nx1, N1);
        r_stack1 = zeros(nu1, K1);
        nx_arr1 = int32(repmat(nx1, 1, N1));
        nu_arr1 = int32([repmat(nu1, 1, K1), 0]);
        nq_arr1 = int32(zeros(1, N1));
        nbx_arr1 = int32([nx1, zeros(1, K1)]);
        nq_per_stage1 = int32(zeros(1, N1));
        % 空 stacks (total_nq=0)
        Qq_stack1 = zeros(nx1, nx1, 0);
        Sq_stack1 = zeros(nu1, nx1, 0);
        Rq_stack1 = zeros(nu1, nu1, 0);
        qq_stack1 = zeros(nx1, 0);
        rq_stack1 = zeros(nu1, 0);
        uq_stack1 = zeros(0, 1);
        x0_1 = [1.0; 0.0];
        idxbx1 = int32(0:nx1-1);
        const1 = 0.0;

        result1 = hpipm_mod.solve_ocp_qcqp(...
            py.numpy.array(A1), py.numpy.array(B1), py.numpy.array(b_stack1), ...
            py.numpy.array(Q_stack1), py.numpy.array(S_stack1), py.numpy.array(R_stack1), ...
            py.numpy.array(q_stack1), py.numpy.array(r_stack1), ...
            py.numpy.array(nx_arr1), py.numpy.array(nu_arr1), ...
            py.numpy.array(nq_arr1), py.numpy.array(nbx_arr1), ...
            py.numpy.array(nq_per_stage1), ...
            py.numpy.array(Qq_stack1), py.numpy.array(Sq_stack1), ...
            py.numpy.array(Rq_stack1), py.numpy.array(qq_stack1), ...
            py.numpy.array(rq_stack1), py.numpy.array(uq_stack1), ...
            py.numpy.array(x0_1), py.numpy.array(idxbx1), ...
            py.float(const1), py.bool(false));

        status1 = double(result1{'status'});
        x1 = double(result1{'x'}); x1 = x1(:);
        iters1 = double(result1{'iters'});
        finite1 = all(isfinite(x1));
        fprintf('  status=%d, iters=%d, x=[%.4f, %.4f, ...], finite=%d\n', ...
            status1, iters1, x1(1), x1(2), finite1);
        if status1 == 0 && finite1
            fprintf('  Test 1: PASS\n');
            n_pass = n_pass + 1;
        else
            fprintf('  Test 1: FAIL (status=%d, finite=%d)\n', status1, finite1);
            n_fail = n_fail + 1;
        end
    catch ME
        fprintf('  Test 1: EXCEPTION: %s\n', ME.message);
        n_fail = n_fail + 1;
    end

    %% ===== Test 2: 一条 terminal wheel constraint =====
    fprintf('\n--- Test 2: 一条 terminal wheel constraint ---\n');
    try
        % K=2, nx=2, nu=1
        % 同 Test 1, 但在 terminal stage (N=2) 加一条二次约束: ||x_N||^2 <= 1
        % Qq = 2*I (×2 约定), qq=0, uq=1
        % nq_per_stage = [0, 0, 1]
        nq_per_stage2 = int32([0, 0, 1]);
        nq_arr2 = nq_per_stage2;
        total_nq2 = 1;
        Qq_stack2 = zeros(nx1, nx1, total_nq2);
        Qq_stack2(:, :, 1) = 2 * eye(nx1);  % terminal Qq
        Sq_stack2 = zeros(nu1, nx1, total_nq2);  % terminal nu=0, 不用但需传 (nu,nx,1)
        Rq_stack2 = zeros(nu1, nu1, total_nq2);
        qq_stack2 = zeros(nx1, total_nq2);
        rq_stack2 = zeros(nu1, total_nq2);
        uq_stack2 = ones(total_nq2, 1);  % ||x_N||^2 <= 1

        result2 = hpipm_mod.solve_ocp_qcqp(...
            py.numpy.array(A1), py.numpy.array(B1), py.numpy.array(b_stack1), ...
            py.numpy.array(Q_stack1), py.numpy.array(S_stack1), py.numpy.array(R_stack1), ...
            py.numpy.array(q_stack1), py.numpy.array(r_stack1), ...
            py.numpy.array(nx_arr1), py.numpy.array(nu_arr1), ...
            py.numpy.array(nq_arr2), py.numpy.array(nbx_arr1), ...
            py.numpy.array(nq_per_stage2), ...
            py.numpy.array(Qq_stack2), py.numpy.array(Sq_stack2), ...
            py.numpy.array(Rq_stack2), py.numpy.array(qq_stack2), ...
            py.numpy.array(rq_stack2), py.numpy.array(uq_stack2), ...
            py.numpy.array(x0_1), py.numpy.array(idxbx1), ...
            py.float(const1), py.bool(false));

        status2 = double(result2{'status'});
        x2 = double(result2{'x'}); x2 = x2(:);
        iters2 = double(result2{'iters'});
        finite2 = all(isfinite(x2));
        % 提取 terminal state 检查约束 (用 x_stages 字段, Python list → MATLAB cell)
        x_stages2 = result2{'x_stages'};
        x_N = double(x_stages2{end});  % terminal state x_N (最后 stage)
        x_N = x_N(:);  % 强制列向量
        constraint_val = 0.5 * x_N' * (2 * eye(nx1)) * x_N - 1;  % ||x_N||^2 - 1
        fprintf('  status=%d, iters=%d, x_N=[%.4f, %.4f], constraint_val=%.6e, finite=%d\n', ...
            status2, iters2, x_N(1), x_N(2), constraint_val, finite2);
        if status2 == 0 && finite2 && constraint_val < 1e-6
            fprintf('  Test 2: PASS\n');
            n_pass = n_pass + 1;
        else
            fprintf('  Test 2: FAIL (status=%d, finite=%d, constraint_val=%.6e)\n', ...
                status2, finite2, constraint_val);
            n_fail = n_fail + 1;
        end
    catch ME
        fprintf('  Test 2: EXCEPTION: %s\n', ME.message);
        n_fail = n_fail + 1;
    end

    %% ===== Test 2b: stage 0 纯 Qq 约束 (无 Sq cross term, 区分 status=3 根因) =====
    fprintf('\n--- Test 2b: stage 0 纯 Qq 约束 (无 Sq cross term) ---\n');
    try
        % K=2, nx=2, nu=1
        % 在 stage 0 加一条纯 Qq 约束 (Sq=0, Rq=0):
        %   0.5*x'Qq*x + qq'*x <= uq  (只含 x, 不含 u)
        % Qq = 2*I (×2), Sq=0, Rq=0, qq=0, rq=0, uq=10
        % 约束: x(1)^2 + x(2)^2 <= 10
        nq_per_stage2b = int32([1, 0, 0]);
        nq_arr2b = nq_per_stage2b;
        total_nq2b = 1;
        Qq_stack2b = zeros(nx1, nx1, total_nq2b);
        Qq_stack2b(:, :, 1) = 2 * eye(nx1);
        Sq_stack2b = zeros(nu1, nx1, total_nq2b);  % 全零 (无 cross term)
        Rq_stack2b = zeros(nu1, nu1, total_nq2b);  % 全零
        qq_stack2b = zeros(nx1, total_nq2b);
        rq_stack2b = zeros(nu1, total_nq2b);
        uq_stack2b = 10 * ones(total_nq2b, 1);

        result2b = hpipm_mod.solve_ocp_qcqp(...
            py.numpy.array(A1), py.numpy.array(B1), py.numpy.array(b_stack1), ...
            py.numpy.array(Q_stack1), py.numpy.array(S_stack1), py.numpy.array(R_stack1), ...
            py.numpy.array(q_stack1), py.numpy.array(r_stack1), ...
            py.numpy.array(nx_arr1), py.numpy.array(nu_arr1), ...
            py.numpy.array(nq_arr2b), py.numpy.array(nbx_arr1), ...
            py.numpy.array(nq_per_stage2b), ...
            py.numpy.array(Qq_stack2b), py.numpy.array(Sq_stack2b), ...
            py.numpy.array(Rq_stack2b), py.numpy.array(qq_stack2b), ...
            py.numpy.array(rq_stack2b), py.numpy.array(uq_stack2b), ...
            py.numpy.array(x0_1), py.numpy.array(idxbx1), ...
            py.float(const1), py.bool(false));

        status2b = double(result2b{'status'});
        x2b = double(result2b{'x'}); x2b = x2b(:);
        iters2b = double(result2b{'iters'});
        finite2b = all(isfinite(x2b));
        fprintf('  status=%d, iters=%d, finite=%d\n', status2b, iters2b, finite2b);
        if status2b == 0 && finite2b
            fprintf('  Test 2b: PASS (纯 Qq 约束在 stage 0 正常 → status=3 与 Sq cross term 有关)\n');
            n_pass = n_pass + 1;
        else
            fprintf('  Test 2b: FAIL (status=%d) → status=3 与任何 stage 0 二次约束有关\n', status2b);
            n_fail = n_fail + 1;
        end
    catch ME
        fprintf('  Test 2b: EXCEPTION: %s\n', ME.message);
        n_fail = n_fail + 1;
    end

    %% ===== Test 3: 一条 stage steering constraint (含 cross term) =====
    fprintf('\n--- Test 3: 一条 stage steering constraint (含 cross term) ---\n');
    try
        % K=2, nx=2, nu=1
        % 在 stage 0 加一条含 cross term 的约束:
        %   0.5*x'Qq*x + x'Sq'*u + 0.5*u'Rq*u + qq'*x + rq'*u <= uq
        % Qq = 2*I (×2), Sq = [1, 0; ...] (nu×nx = 1×2), Rq = 1, qq = 0, rq = 0, uq = 10
        % 约束: x(1)^2 + x(2)^2 + x(1)*u + 0.5*u^2 <= 10
        nq_per_stage3 = int32([1, 0, 0]);
        nq_arr3 = nq_per_stage3;
        total_nq3 = 1;
        Qq_stack3 = zeros(nx1, nx1, total_nq3);
        Qq_stack3(:, :, 1) = 2 * eye(nx1);
        Sq_stack3 = zeros(nu1, nx1, total_nq3);
        Sq_stack3(:, :, 1) = [1, 0];  % 1×2, Sq(1,1)=1 → x(1)*u
        Rq_stack3 = zeros(nu1, nu1, total_nq3);
        Rq_stack3(:, :, 1) = 1;  % 0.5*u'*1*u = 0.5*u^2
        qq_stack3 = zeros(nx1, total_nq3);
        rq_stack3 = zeros(nu1, total_nq3);
        uq_stack3 = 10 * ones(total_nq3, 1);

        result3 = hpipm_mod.solve_ocp_qcqp(...
            py.numpy.array(A1), py.numpy.array(B1), py.numpy.array(b_stack1), ...
            py.numpy.array(Q_stack1), py.numpy.array(S_stack1), py.numpy.array(R_stack1), ...
            py.numpy.array(q_stack1), py.numpy.array(r_stack1), ...
            py.numpy.array(nx_arr1), py.numpy.array(nu_arr1), ...
            py.numpy.array(nq_arr3), py.numpy.array(nbx_arr1), ...
            py.numpy.array(nq_per_stage3), ...
            py.numpy.array(Qq_stack3), py.numpy.array(Sq_stack3), ...
            py.numpy.array(Rq_stack3), py.numpy.array(qq_stack3), ...
            py.numpy.array(rq_stack3), py.numpy.array(uq_stack3), ...
            py.numpy.array(x0_1), py.numpy.array(idxbx1), ...
            py.float(const1), py.bool(true));  % verbose=true 诊断 status=3

        status3 = double(result3{'status'});
        x3 = double(result3{'x'}); x3 = x3(:);
        iters3 = double(result3{'iters'});
        finite3 = all(isfinite(x3));
        % 检查约束: x_0 = [1, 0], u_0 = u_stages{1} (Python list → MATLAB cell)
        x_0 = x0_1;
        u_stages3 = result3{'u_stages'};
        u_0 = double(u_stages3{1});  % u at stage 0
        u_0 = u_0(:);
        constraint_val = 0.5 * x_0' * (2 * eye(nx1)) * x_0 + x_0' * [1; 0] * u_0 ...
            + 0.5 * u_0' * 1 * u_0 - 10;
        fprintf('  status=%d, iters=%d, u_0=%.4f, constraint_val=%.6e, finite=%d\n', ...
            status3, iters3, u_0(1), constraint_val, finite3);
        if status3 == 0 && finite3 && constraint_val < 1e-6
            fprintf('  Test 3: PASS\n');
            n_pass = n_pass + 1;
        else
            fprintf('  Test 3: FAIL (status=%d, finite=%d, constraint_val=%.6e)\n', ...
                status3, finite3, constraint_val);
            n_fail = n_fail + 1;
        end
    catch ME
        fprintf('  Test 3: EXCEPTION: %s\n', ME.message);
        n_fail = n_fail + 1;
    end

    %% ===== Test 4 & 5: RSS 问题 (wheel only / full 72) =====
    fprintf('\n--- Test 4 & 5: RSS 问题 ---\n');
    try
        params = config();
        K = 6; nx = 6; nu = 3;
        path = generateReference(params, params.num_path_pts);
        v0 = [0.01; 0.01; 0.01];
        state = [0.05, 0.1, 0.2]';
        step = 1;
        u_anchor = zeros(nu, K);  % u_hat = 0 (与 Dense QCQP 一致)

        ocp = construct_ocp_qcqp_from_rss(path, step, v0, state, u_anchor, params);

        % 数据格式转换
        b_stack = zeros(nx, ocp.K);
        r_stack = zeros(nu, ocp.K);
        for n = 1:ocp.K
            b_stack(:, n) = ocp.b{n};
            r_stack(:, n) = ocp.r{n};
        end

        % --- Test 4: 仅 wheel 约束 (steering 约束放宽为平凡) ---
        fprintf('\n--- Test 4: 仅 wheel 约束 (24 条) ---\n');
        % 策略: 复制 ocp 的 stacks, 但将 steering 约束的 uq 设为大正值 (放宽)
        Qq_stack4 = ocp.Qq_stack;
        Sq_stack4 = ocp.Sq_stack;
        Rq_stack4 = ocp.Rq_stack;
        qq_stack4 = ocp.qq_stack;
        rq_stack4 = ocp.rq_stack;
        uq_stack4 = ocp.uq_stack;
        for j = 1:size(ocp.uq_stack, 1)
            if strcmp(ocp.metadata.kind{j}, 'steering')
                % 放宽 steering 约束: uq = +1e10 (几乎无约束)
                uq_stack4(j) = 1e10;
            end
        end

        result4 = hpipm_mod.solve_ocp_qcqp(...
            py.numpy.array(ocp.A), py.numpy.array(ocp.B), ...
            py.numpy.array(b_stack), ...
            py.numpy.array(ocp.Q_eff), py.numpy.array(ocp.S_eff), ...
            py.numpy.array(ocp.R_eff), py.numpy.array(ocp.q_stack), ...
            py.numpy.array(r_stack), ...
            py.numpy.array(int32(ocp.nx)), py.numpy.array(int32(ocp.nu)), ...
            py.numpy.array(int32(ocp.nq)), py.numpy.array(int32(ocp.nbx)), ...
            py.numpy.array(int32(ocp.nq_per_stage)), ...
            py.numpy.array(Qq_stack4), py.numpy.array(Sq_stack4), ...
            py.numpy.array(Rq_stack4), py.numpy.array(qq_stack4), ...
            py.numpy.array(rq_stack4), py.numpy.array(uq_stack4), ...
            py.numpy.array(ocp.x0), py.numpy.array(int32(ocp.idxbx)), ...
            py.float(ocp.const), py.bool(true));  % verbose=true 诊断 status=3

        status4 = double(result4{'status'});
        x4 = double(result4{'x'}); x4 = x4(:);
        iters4 = double(result4{'iters'});
        finite4 = all(isfinite(x4));
        fprintf('  status=%d, iters=%d, finite=%d\n', status4, iters4, finite4);
        if status4 == 0 && finite4
            fprintf('  Test 4: PASS\n');
            n_pass = n_pass + 1;
        else
            fprintf('  Test 4: FAIL (status=%d)\n', status4);
            n_fail = n_fail + 1;
        end

        % --- Test 5: 全部 72 条约束 ---
        fprintf('\n--- Test 5: 全部 72 条约束 ---\n');
        result5 = hpipm_mod.solve_ocp_qcqp(...
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
            py.float(ocp.const), py.bool(true));  % verbose=true 诊断 status=3

        status5 = double(result5{'status'});
        x5 = double(result5{'x'}); x5 = x5(:);
        iters5 = double(result5{'iters'});
        obj5 = double(result5{'obj_value_full'});
        finite5 = all(isfinite(x5));
        fprintf('  status=%d, iters=%d, obj=%.6f, finite=%d\n', status5, iters5, obj5, finite5);
        if status5 == 0 && finite5
            fprintf('  Test 5: PASS\n');
            n_pass = n_pass + 1;
        else
            fprintf('  Test 5: FAIL (status=%d)\n', status5);
            n_fail = n_fail + 1;
        end
    catch ME
        fprintf('  Test 4/5: EXCEPTION: %s\n', ME.message);
        for si = 1:min(length(ME.stack), 3)
            fprintf('    at %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
        end
        n_fail = n_fail + 2;
    end

    %% ===== Summary =====
    fprintf('\n========== Summary ==========\n');
    fprintf('Pass: %d / %d\n', n_pass, n_pass + n_fail);
    if n_fail == 0
        fprintf('结论: PASS (全部 smoke test 通过)\n');
    else
        fprintf('结论: FAIL (%d 个 test 失败)\n', n_fail);
    end
    fprintf('==============================================================\n');
end
