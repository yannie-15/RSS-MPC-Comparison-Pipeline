function [u, new_state_dot, velocity, diagnostics] = control_RSS(path, step, state_dot, state)
% CONTROL_RSS  RSS proposed 算法控制器 (HPIPM dense QCQP, 3 次 SQP 迭代)
%
% 架构 (与参考仓库一致):
%   1. MATLAB 端调用 construct_complete_qp_from_rss 构造 QP 矩阵
%      (H, g, A, b, Hq, gq, uq), 修复了 Bug1-5
%   2. Python 端 hpipm_qp_solver.solve_qcqp 仅负责求解
%   3. SQP 外层循环 (max_iter=3) 保留在 MATLAB 端
%
% 输入:
%   path      : 3xN 参考轨迹
%   step      : 当前仿真步 (1-based)
%   state_dot : 3x1 上一步车体系速度 (last_vel)
%   state     : 1x3 当前状态 [x, y, psi]
%
% 输出:
%   u            : 3xK 控制增量矩阵 (最后可行解)
%   new_state_dot: 3x1 世界系状态导数
%   velocity     : 3x1 车体系速度
%   diagnostics  : struct 每次子迭代诊断

    params = config();

    % ================= Param Setup =================
    K = 6; rho = 0.01; k1 = 1; epsilon = 0;  % k1=1 硬编码, 匹配原始 CVX 0121 分支
    current_xy = [state(1), state(2)]';
    psi0 = state(3); current_nu = state_dot;

    % ================= 迭代 Setup =================
    max_iter = 3;  % 子问题迭代次数
    u_hat = zeros(3, K);

    global solver_time_array;
    if ~exist('solver_time_array', 'var') || isempty(solver_time_array)
        solver_time_array = [];
    end

    % Phase 6: 诊断结构体
    diagnostics = struct();
    diagnostics.iterations = struct();
    diagnostics.iterations.status = cell(1, max_iter);
    diagnostics.iterations.optval = zeros(1, max_iter);
    diagnostics.iterations.solve_time = NaN(1, max_iter);
    diagnostics.iterations.solver_name = cell(1, max_iter);
    diagnostics.step = step;
    diagnostics.max_iter = max_iter;

    % ================= Python 环境路径设置 =================
    persistent py_path_added py_reloaded;
    if isempty(py_path_added)
        script_path = fileparts(mfilename('fullpath'));
        % hpipm_qp_solver.py 现与 control_RSS.m 同目录 (algorithms/RSS_proposed/)
        if exist(script_path, 'dir')
            sys_mod = py.importlib.import_module('sys');
            py.getattr(sys_mod, 'path').append(script_path);
        end
        py_path_added = true;
    end
    % 首次调用时重新加载 Python 模块, 确保使用最新的 .py 代码
    if isempty(py_reloaded)
        try
            solver_mod = py.importlib.import_module('hpipm_qp_solver');
            py.importlib.reload(solver_mod);
            py_reloaded = true;
        catch
        end
    end

    % ================= SQP 外层循环 =================
    % 与原始 CVX 0121 分支对齐: 无论成功失败都更新 u_hat, 不 break, 不回退
    for m = 1 : max_iter
        try
            % ========== MATLAB 端构造 QP (使用参考仓库验证过的代码) ==========
            qp = construct_complete_qp_from_rss(path, step, current_nu, state, u_hat, params);

            % ========== 把 Hq/gq cell 数组转为 numpy 数组 ==========
            % Hq_list 是 cell array of (n_var x n_var), 堆叠为 3D 数组 (n_var, n_var, n_qcqp)
            n_var = qp.n_var;
            n_qcqp = qp.n_qcqp;
            Hq_3d = zeros(n_var, n_var, n_qcqp);
            gq_2d = zeros(n_var, n_qcqp);
            for i = 1:n_qcqp
                Hq_3d(:, :, i) = qp.Hq{i};
                gq_2d(:, i) = qp.gq{i};
            end

            % ========== 调用 Python HPIPM 求解器 ==========
            result = py.hpipm_qp_solver.solve_qcqp(...
                py.numpy.array(qp.H), ...
                py.numpy.array(qp.g), ...
                py.numpy.array(qp.A), ...
                py.numpy.array(qp.b), ...
                py.numpy.array(Hq_3d), ...
                py.numpy.array(gq_2d), ...
                py.numpy.array(qp.uq'), ...
                py.bool(false) ...
            );

            % ========== 提取结果 ==========
            x = double(result{'x'});           % (n_var,)
            status_code = double(result{'status'});
            optval = double(result{'obj_value'});
            inner_solve_time = double(result{'solve_time'});

            % 从 x 提取 u (前 3*K 个变量)
            u_sol = reshape(x(1:3*K), 3, K);

            % HPIPM status: 0=SUCCESS, 其他=失败
            if status_code == 0
                cvx_status_str = 'Solved';
            else
                cvx_status_str = 'Failed';
            end
            solver_name = 'HPIPM';

        catch ME
            % Python 调用失败
            u_sol = zeros(3, K);
            status_code = -1;
            optval = NaN;
            inner_solve_time = NaN;
            cvx_status_str = 'Failed';
            solver_name = 'HPIPM-Error';
            fprintf('第%d步第%d次迭代 - Python 调用异常: %s\n', step, m, ME.message);
        end

        % ========== 存入全局数组 ==========
        solver_time_array(3*step + m - 3) = inner_solve_time;

        % ========== 记录诊断 ==========
        diagnostics.iterations.status{m} = cvx_status_str;
        diagnostics.iterations.optval(m) = optval;
        diagnostics.iterations.solve_time(m) = inner_solve_time;
        diagnostics.iterations.solver_name{m} = solver_name;

        % 打印结果
        fprintf('第%d步第%d次迭代 - %s内部求解时间：%.6f秒 (status=%d)\n', ...
            step, m, solver_name, inner_solve_time, status_code);
        fprintf('最优代价: %f | 求解状态: %s\n', optval, cvx_status_str);

        % ========== 迭代更新 (与原始 CVX 0121 一致: 无条件更新) ==========
        u = u_sol;
        u_hat = u;
    end

    % 输出状态导数
    new_state_dot =  [cos(state(3)), -sin(state(3)), 0;
                     sin(state(3)),  cos(state(3)), 0;
                         0,              0, 1] * (state_dot + 1.00 * u(:, 1));
    velocity = current_nu + u(:, 1);

    % 汇总诊断
    diagnostics.total_solve_time = sum(diagnostics.iterations.solve_time(~isnan( ...
        diagnostics.iterations.solve_time)));
    if isempty(diagnostics.total_solve_time)
        diagnostics.total_solve_time = NaN;
    end
end
