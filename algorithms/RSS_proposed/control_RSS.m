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
%          → 调用 construct_ocp_qp_from_rss 构造 OCP QCQP (误差状态, 逐阶段)
%       2. 求解 u^(m+1) = S(u^(m)) (论文 Algorithm 1 line 5)
%          → 调用 Python HPIPM OCP QCQP 求解器 (ocp_qcqp 接口)
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
% 求解器: 论文原用 CVX+ECOS, 本实现替换为 HPIPM OCP QCQP (Python ocp_qcqp 接口)
%         与 dense QCQP 数学等价, 但利用 OCP 块三对角结构更高效

    params = config();

    % ================= Param Setup =================
    % 论文 IV-A: K=6, dt=0.01s, t_end=1s
    K = 6; rho = 0.01; k1 = 1; epsilon = 0;  % k1=1: 论文 Alg.1 line 11 增益
    current_xy = [state(1), state(2)]';       % 当前位置 (世界系)
    psi0 = state(3); v0 = state_dot;   % 当前航向 / ν_0 (已知量: 当前车体系速度)

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
            % 用临时 .py 文件执行 sys.modules.pop + import (比 reload 更可靠)
            % reload 不会重新解析文件路径, 必须先 pop 再 import
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
    % 论文 Theorem 2: 序列 {u^(m)} 单调下降且收敛到 P_K 的驻点
    % 无条件更新 u_hat (论文 Alg.1 无 break-on-failure, descent inequality 保证下降)
    for m = 1 : max_iter
        try

            % ========== 论文 Alg.1 line 4: 构造凸子问题 Q_K(u^(m)) ==========
            % Dense QCQP (精确凸化二次约束, HPIPM dense_qcqp)
            % 替代 OCP QP + 线性化 (后者在 u_hat=0 处退化)
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
            % x = [u(:); nu(:)] (36维); x 的前 3K 个变量 = paper u
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
            fprintf('  错误位置: %s (line %d) in %s\n', ME.stack(1).name, ME.stack(1).line, ME.stack(1).file);
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
        % 论文: u^(m+1) = S(u^(m)), 即 û ← u_sol 用于下次迭代凸化
        u = u_sol;
        % 求解失败时 u_sol 可能含 NaN/Inf, 直接赋给 u_hat 会导致下一轮
        % construct_complete_qp_from_rss 构造出全 NaN 的 QCQP, 失败级联传播.
        % 仅在解有限时更新 u_hat, 否则保留上一轮 û (论文 Alg.1 假设始终可解,
        % 数值失败时需防止 NaN 污染下一轮凸化)
        if status_code == 0 && all(isfinite(u_sol(:)))
            u_hat = u;
        end
    end

    % ================= 论文 Alg.1 line 11: 输出 ν_1 = ν_0 + u_1 =================
    % 论文 (1) 式: ξ̇_w = [R(ψ_w), 0; 0, 1] * ξ̇_c
    % 这里 new_state_dot = R(ψ) * (ν_0 + k1*u_1), k1=1
    new_state_dot =  [cos(state(3)), -sin(state(3)), 0;
                     sin(state(3)),  cos(state(3)), 0;
                         0,              0, 1] * (state_dot + 1.00 * u(:, 1));
    % 车体系速度 ν_1 = ν_0 + u_1 (论文 Alg.1 line 11)
    velocity = v0 + u(:, 1);

    % ================= [DIAGNOSTIC] 原始约束违反量检查 =================
    % 对最终 u_hat (3 次 SCP 后输出) 计算完整 nu 轨迹, 回代验证:
    %   轮速 SOC: ||H_n * nu_k||_2 <= vimax?  违反量 = max(0, ||Hn*nu_k|| - vimax)
    %   转向锥:   nu_km1' H' R1/R2 H nu_k >= 0? 违反量 = max(0, -term)
    persistent diag_stats;
    if isempty(diag_stats) || step == 1
        diag_stats = struct('max_wheel_viol', 0, 'max_cone_viol', 0, ...
            'cnt_wheel_viol', 0, 'cnt_cone_viol', 0, 'step_wheel_max', 0, 'step_cone_max', 0);
    end
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
                term = a' * Mgg * b;   % 原约束: term >= 0
                viol = -term;          % 违反量 = max(0, -term)
                if viol > mv_cone; mv_cone = viol; end
            end
        end
    end
    % ---- 汇总统计 ----
    if mv_wheel > diag_stats.max_wheel_viol; diag_stats.max_wheel_viol = mv_wheel; diag_stats.step_wheel_max = step; end
    if mv_cone > diag_stats.max_cone_viol; diag_stats.max_cone_viol = mv_cone; diag_stats.step_cone_max = step; end
    if mv_wheel > 1e-6; diag_stats.cnt_wheel_viol = diag_stats.cnt_wheel_viol + 1; end
    if mv_cone > 1e-6; diag_stats.cnt_cone_viol = diag_stats.cnt_cone_viol + 1; end
    % 每 25 步打印一次, 最后一步打印汇总
    if mod(step, 25) == 0 || step == 1
        fprintf('[DIAG step=%d] wheel_viol_cur=%.6f, cone_viol_cur=%.6f\n', step, mv_wheel, mv_cone);
    end
    if step == 100
        fprintf('\n============== [DIAG 汇总 100步约束违反量] ==============\n');
        fprintf('轮速 SOC 最大违反量: %.6f (发生在 step %d) | 违反步数: %d/100\n', ...
            diag_stats.max_wheel_viol, diag_stats.step_wheel_max, diag_stats.cnt_wheel_viol);
        fprintf('转向锥 最大违反量: %.6f (发生在 step %d) | 违反步数: %d/100\n', ...
            diag_stats.max_cone_viol, diag_stats.step_cone_max, diag_stats.cnt_cone_viol);
        if diag_stats.max_wheel_viol < 1e-6 && diag_stats.max_cone_viol < 1e-6
            fprintf('结论: 解在原始约束下可行 (所有违反量 < 1e-6)\n');
        elseif diag_stats.max_wheel_viol < 1e-3 && diag_stats.max_cone_viol < 1e-3
            fprintf('结论: 解在原始约束下基本可行 (违反量 < 1e-3, 数值误差量级)\n');
        else
            fprintf('结论: 解严重违反原始约束 → 线性化问题过度松弛, J_total 偏低源于不可行!\n');
        end
        fprintf('=========================================================\n\n');
    end
    % ================= [END DIAGNOSTIC] =================
    % 汇总诊断 (论文 IV-B: computation cost 记录)
    diagnostics.total_solve_time = sum(diagnostics.iterations.solve_time(~isnan( ...
        diagnostics.iterations.solve_time)));
    if isempty(diagnostics.total_solve_time)
        diagnostics.total_solve_time = NaN;
    end
end
