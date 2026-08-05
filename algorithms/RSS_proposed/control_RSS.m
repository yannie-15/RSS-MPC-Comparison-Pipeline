function [u, new_state_dot, velocity, diagnostics] = control_RSS(path, step, state_dot, state)
% CONTROL_RSS  论文 Algorithm 1 (Trajectory optimizer for SWMRs) 的实现
%
% 论文: RSS26 "Exploit Agile Mobility of Steerable-Wheeled Mobile Robots:
%        A Fast Motion Planning Approach"
%
% 本函数对应论文 Algorithm 1 的全部流程:
%   - 输入: 当前速度 ν_0 (state_dot), 当前位姿 state (ξ_w)
%   - 初始化 u^(0) (论文建议 static init: u^(0)=0)
%   - while m < max_iter (论文 Algorithm 1 line 3-9):
%       1. 构造凸子问题 Q_K(u^(m)) (论文公式 17)
%          → 调用 construct_complete_qp_from_rss 构造 QCQP 矩阵
%       2. 求解 u^(m+1) = S(u^(m)) (论文 Algorithm 1 line 5)
%          → 调用 Python HPIPM dense QCQP 求解器
%       3. 更新 u_hat = u (论文 Algorithm 1 line 9)
%   - 输出: ν_1 = ν_0 + u_1 (论文 Algorithm 1 line 11)
%
% 论文中的符号对照:
%   state_dot   = ν_0 (当前车体系速度, 论文 (8) 式)
%   u           = u ∈ R^{3×K} (优化变量, 控制增量序列, 论文 (9) 式)
%   u_hat       = û (上一次迭代解, 用于凸化, 论文 Prop.1)
%   K = 6       = 预测时域 (论文 IV-A: "prediction horizon K=6")
%   rho = 0.01  = ρ (强凸正则化参数, 论文 (17) 式)
%   k1 = 1      = 输出增益 (论文 Algorithm 1 line 11: ν_1 = ν_0 + u_1)
%
% 求解器: 论文原用 CVX+ECOS, 本实现替换为 HPIPM dense QCQP (Python 接口)

    params = config();

    % ================= Param Setup =================
    % 论文 IV-A: K=6, dt=0.01s, t_end=1s
    K = 6; rho = 0.01; k1 = 1; epsilon = 0;  % k1=1: 论文 Alg.1 line 11 增益
    current_xy = [state(1), state(2)]';       % 当前位置 (世界系)
    psi0 = state(3); current_nu = state_dot;   % 当前航向 / ν_0

    % ================= 迭代 Setup =================
    % 论文 Algorithm 1 line 1: Initialize u^(0) ∈ ri(D(P_K))
    % 采用 static initialization (论文: "set u^(0)=0")
    max_iter = 3;  % 论文 IV-B: "maximum number of iterations is set to 3"
    u_hat = zeros(3, K);  % u^(0) = 0 (static init)

    global solver_time_array;
    if ~exist('solver_time_array', 'var') || isempty(solver_time_array)
        solver_time_array = [];
    end

    % 诊断结构体 (记录每次迭代, 对应论文 Alg.1 的 m 循环)
    diagnostics = struct();
    diagnostics.iterations = struct();
    diagnostics.iterations.status = cell(1, max_iter);
    diagnostics.iterations.optval = zeros(1, max_iter);
    diagnostics.iterations.solve_time = NaN(1, max_iter);
    diagnostics.iterations.solver_name = cell(1, max_iter);
    diagnostics.step = step;
    diagnostics.max_iter = max_iter;

    % ================= Python 环境路径设置 =================
    % persistent: 保证 Python 模块只加载一次, 避免 libhpipm.dll 内存泄漏
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
            solver_mod = py.importlib.import_module('hpipm_qp_solver');
            py.importlib.reload(solver_mod);
            py_reloaded = true;
        catch
        end
    end

    % ================= SQP 外层循环 (论文 Algorithm 1 line 3-9) =================
    % 论文 Theorem 2: 序列 {u^(m)} 单调下降且收敛到 P_K 的驻点
    % 无条件更新 u_hat (论文 Alg.1 无 break-on-failure, descent inequality 保证下降)
    for m = 1 : max_iter
        try
            % ========== 论文 Alg.1 line 4: 构造凸子问题 Q_K(u^(m)) ==========
            % 论文公式 (17): min f(u,û)=g(u)+ρ||u-û||^2  s.t. C_k_{i,n}<=0, u∈∩U_j
            % 这里将 Q_K 转为 dense QCQP 标准形式 (H,g,A,b,Hq,gq,uq)
            qp = construct_complete_qp_from_rss(path, step, current_nu, state, u_hat, params);

            % ========== 数据格式转换: cell → 3D numpy ==========
            % HPIPM Python 接口要求 Hq 为 3D 数组 (n_var, n_var, n_qcqp)
            n_var = qp.n_var;
            n_qcqp = qp.n_qcqp;
            Hq_3d = zeros(n_var, n_var, n_qcqp);
            gq_2d = zeros(n_var, n_qcqp);
            for i = 1:n_qcqp
                Hq_3d(:, :, i) = qp.Hq{i};
                gq_2d(:, i) = qp.gq{i};
            end

            % ========== 论文 Alg.1 line 5: 求解 u^(m+1) = S(u^(m)) ==========
            % 调用 HPIPM dense QCQP 求解器 (替代论文中的 CVX+ECOS)
            % HPIPM 论文: "hpipm: a high-performance quadratic programming framework"
            %   求解: min 0.5*x'Hx + g'x  s.t. Ax=b, 0.5*x'Hq_i*x + gq_i'x <= uq_i
            result = py.hpipm_qp_solver.solve_qcqp(...
                py.numpy.array(qp.H), ...        % H: Hessian (论文 (18)+(19) 展开)
                py.numpy.array(qp.g), ...        % g: 线性项
                py.numpy.array(qp.A), ...        % A: 等式约束 (论文 (20c) 动力学)
                py.numpy.array(qp.b), ...        % b: 等式约束右端
                py.numpy.array(Hq_3d), ...       % Hq: 二次约束 Hessian (论文 (20a)+(20b))
                py.numpy.array(gq_2d), ...       % gq: 二次约束线性项
                py.numpy.array(qp.uq'), ...      % uq: 二次约束右端
                py.bool(false) ...              % verbose
            );

            % ========== 提取结果 ==========
            % x = [u(:); nu(:)] (36维), 论文决策变量 u ∈ R^{3×K}
            x = double(result{'x'});           % (n_var,) = (36,)
            status_code = double(result{'status'});
            optval = double(result{'obj_value'});
            inner_solve_time = double(result{'solve_time'});

            % 从 x 提取 u (前 3*K 个变量, 论文 u={u_1,...,u_K})
            u_sol = reshape(x(1:3*K), 3, K);

            % HPIPM status: 0=SUCCESS (对应论文 S(û) 存在唯一解)
            if status_code == 0
                cvx_status_str = 'Solved';
            else
                cvx_status_str = 'Failed';
            end
            solver_name = 'HPIPM';

        catch ME
            % Python 调用失败 (对应 Q_K 不可解的情况, 但论文 Prop.2 保证可行性)
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

        fprintf('第%d步第%d次迭代 - %s内部求解时间：%.6f秒 (status=%d)\n', ...
            step, m, solver_name, inner_solve_time, status_code);
        fprintf('最优代价: %f | 求解状态: %s\n', optval, cvx_status_str);

        % ========== 论文 Alg.1 line 9: 更新 m ← m+1 ==========
        % 论文: u^(m+1) = S(u^(m)), 即 û ← u_sol 用于下次迭代凸化
        u = u_sol;
        u_hat = u;  % 更新 û, 下次构造 Q_K 时使用
    end

    % ================= 论文 Alg.1 line 11: 输出 ν_1 = ν_0 + u_1 =================
    % 论文 (1) 式: ξ̇_w = [R(ψ_w), 0; 0, 1] * ξ̇_c
    % 这里 new_state_dot = R(ψ) * (ν_0 + k1*u_1), k1=1
    new_state_dot =  [cos(state(3)), -sin(state(3)), 0;
                     sin(state(3)),  cos(state(3)), 0;
                         0,              0, 1] * (state_dot + 1.00 * u(:, 1));
    % 车体系速度 ν_1 = ν_0 + u_1 (论文 Alg.1 line 11)
    velocity = current_nu + u(:, 1);

    % 汇总诊断 (论文 IV-B: computation cost 记录)
    diagnostics.total_solve_time = sum(diagnostics.iterations.solve_time(~isnan( ...
        diagnostics.iterations.solve_time)));
    if isempty(diagnostics.total_solve_time)
        diagnostics.total_solve_time = NaN;
    end
end
