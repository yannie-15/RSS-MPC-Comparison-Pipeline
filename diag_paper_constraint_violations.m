function diag_paper_constraint_violations()
% DIAG_PAPER_CONSTRAINT_VIOLATIONS
% 复现 paper_fixed 场景 (proposed-3iter), 对每步的完整预测轨迹 (6 步)
% 回代验证原始二次约束 (轮速 SOC / 转向锥) 的违反量。
%
% 用法：
%   cd d:\PROJECT\RSS-MPC-Comparison-Pipeline-rss_hpipm
%   setup_paths();
%   diag_paper_constraint_violations();

    proj_root = 'd:\PROJECT\RSS-MPC-Comparison-Pipeline-rss_hpipm';
    rss_proposed_dir = fullfile(proj_root, 'algorithms', 'RSS_proposed');
    addpath(fullfile(proj_root, 'paper_reproduction'));
    addpath(fullfile(proj_root, 'core'));
    addpath(fullfile(proj_root, 'batch_simulation'));
    addpath(fullfile(proj_root, 'algorithms'));
    addpath(rss_proposed_dir);   % config.m / control_RSS.m 在此目录
    setup_paths();

    % 清理 Python sys.path 中的 RSS_V2 残留
    cleanup_py = fullfile(proj_root, '_cleanup_syspath.py');
    fid = fopen(cleanup_py, 'w');
    fprintf(fid, 'import sys\n');
    fprintf(fid, 'to_remove = [p for p in list(sys.path) if "RSS_V2" in p or "Projects\\\\RSS" in p]\n');
    fprintf(fid, 'for p in to_remove:\n');
    fprintf(fid, '    sys.path.remove(p)\n');
    fprintf(fid, 'target = r"%s"\n', rss_proposed_dir);
    fprintf(fid, 'if target not in sys.path:\n');
    fprintf(fid, '    sys.path.append(target)\n');
    fclose(fid);
    py.runpy.run_path(cleanup_py);
    delete(cleanup_py);
    % 清 pycache
    pycache_dir = fullfile(rss_proposed_dir, '__pycache__');
    if exist(pycache_dir, 'dir')
        pyc_files = dir(fullfile(pycache_dir, 'hpipm_qp_solver*.pyc'));
        for i = 1:length(pyc_files); delete(fullfile(pycache_dir, pyc_files(i).name)); end
    end
    % 重载 Python 模块
    reload_py = fullfile(proj_root, '_reload_solver.py');
    fid = fopen(reload_py, 'w');
    fprintf(fid, 'import sys\n');
    fprintf(fid, 'sys.modules.pop("hpipm_qp_solver", None)\n');
    fprintf(fid, 'import hpipm_qp_solver\n');
    fclose(fid);
    py.runpy.run_path(reload_py);
    delete(reload_py);
    clear functions;

    %% ======== 场景初始化 (paper_fixed 全零初值) ========
    params = config();
    K_pred = 6;
    num_steps = params.num_steps;   % 100
    num_wheels = size(params.wheel_pos, 1);
    vimax = params.vimax;
    delta_theta = params.dt * params.phidotmax;

    % 轮特征矩阵 Hn
    Hn = cell(1, num_wheels);
    for n = 1:num_wheels
        Hn{n} = [1, 0, -params.wheel_pos(n,2); 0, 1, params.wheel_pos(n,1)];
    end
    % 转向锥旋转矩阵 R1, R2
    R1 = [sin(delta_theta), -cos(delta_theta); cos(delta_theta),  sin(delta_theta)];
    R2 = R1';
    R_stack = cat(3, R1, R2);

    state = [0; 0; 0];
    lastBodyVelocity = [0; 0; 0];
    path = generateReference(params, params.num_path_pts);

    % ======== 全局统计变量 ========
    max_wheel_pred = 0;  step_wheel_pred = 0;
    max_cone_pred  = 0;  step_cone_pred  = 0;
    cnt_wheel_viol = 0;  cnt_cone_viol  = 0;
    max_wheel_closed = 0; max_cone_closed = 0;
    cnt_wheel_closed_viol = 0; cnt_cone_closed_viol = 0;

    closed_body_vel = zeros(3, num_steps+1);
    closed_body_vel(:, 1) = lastBodyVelocity;
    executed_u_closed = zeros(3, num_steps);

    fprintf('\n============== 开始 100 步诊断 ==============\n');
    fprintf('step | wheel_pred_viol | cone_pred_viol\n');
    fprintf('-----+-----------------+-----------------\n');

    for k = 1 : num_steps
        [u_full, ~, bodyVelocity, ~] = control_RSS(path, k, lastBodyVelocity, state');
        u_exec = u_full(:, 1);

        % 重建预测完整 nu 轨迹 (ν_0 .. ν_K_pred)
        nu_pred = zeros(3, K_pred+1);
        nu_pred(:, 1) = lastBodyVelocity;
        for kk = 1 : K_pred
            nu_pred(:, kk+1) = nu_pred(:, kk) + u_full(:, kk);
        end

        % 预测 轮速 SOC: ||Hn nu_k|| <= vimax, k = 1..K_pred
        mv_wheel = 0;
        for kk = 1 : K_pred
            nu_k = nu_pred(:, kk+1);
            for n = 1 : num_wheels
                vn = norm(Hn{n} * nu_k, 2);
                viol = vn - vimax;
                if viol > mv_wheel; mv_wheel = viol; end
            end
        end

        % 预测 转向锥: nu_{k-1}' H' Rg H nu_k >= 0, k = 1..K_pred
        mv_cone = 0;
        for kk = 1 : K_pred
            a = nu_pred(:, kk);
            b = nu_pred(:, kk+1);
            for n = 1 : num_wheels
                Hn_i = Hn{n};
                for gg = 1 : 2
                    Rgg = R_stack(:,:,gg);
                    Mgg = Hn_i' * Rgg * Hn_i;
                    term = a' * Mgg * b;
                    viol = -term;
                    if viol > mv_cone; mv_cone = viol; end
                end
            end
        end

        % 更新预测统计
        if mv_wheel > max_wheel_pred
            max_wheel_pred = mv_wheel; step_wheel_pred = k;
        end
        if mv_cone > max_cone_pred
            max_cone_pred = mv_cone; step_cone_pred = k;
        end
        if mv_wheel > 1e-6; cnt_wheel_viol = cnt_wheel_viol + 1; end
        if mv_cone  > 1e-6; cnt_cone_viol  = cnt_cone_viol  + 1; end

        if mod(k, 10) == 0 || k == 1
            fprintf('%4d |     %.4e      |     %.4e\n', k, mv_wheel, mv_cone);
        end

        % 推进闭环 (与 paper_baseline_case 一致)
        executed_u_closed(:, k) = u_exec;
        closed_body_vel(:, k+1) = bodyVelocity;
        worldVelocity = [cos(state(3)), -sin(state(3)), 0;
                         sin(state(3)),  cos(state(3)), 0;
                         0,              0,              1] * bodyVelocity;
        state = propagateState(state, worldVelocity, params);
        lastBodyVelocity = bodyVelocity;
    end

    % 闭环执行轨迹诊断
    for k = 1 : num_steps
        nu_k = closed_body_vel(:, k+1);
        mv_wheel = 0;
        for n = 1 : num_wheels
            vn = norm(Hn{n} * nu_k, 2);
            viol = vn - vimax;
            if viol > mv_wheel; mv_wheel = viol; end
        end
        a = closed_body_vel(:, k);
        b = closed_body_vel(:, k+1);
        mv_cone = 0;
        for n = 1 : num_wheels
            Hn_i = Hn{n};
            for gg = 1 : 2
                Rgg = R_stack(:,:,gg);
                Mgg = Hn_i' * Rgg * Hn_i;
                term = a' * Mgg * b;
                viol = -term;
                if viol > mv_cone; mv_cone = viol; end
            end
        end
        if mv_wheel > max_wheel_closed; max_wheel_closed = mv_wheel; end
        if mv_cone  > max_cone_closed;  max_cone_closed  = mv_cone;  end
        if mv_wheel > 1e-6; cnt_wheel_closed_viol = cnt_wheel_closed_viol + 1; end
        if mv_cone  > 1e-6; cnt_cone_closed_viol  = cnt_cone_closed_viol  + 1; end
    end

    fprintf('\n============== [诊断汇总] paper_fixed (proposed-3iter) 100 步 ==============\n');
    fprintf('--- [预测域约束违反] MPC 内部预测轨迹 (每步 6 步时域) ---\n');
    fprintf('  轮速 SOC:   max viol = %.6f (step %d), 违反步数 = %d/100 (阈值 1e-6)\n', ...
        max_wheel_pred, step_wheel_pred, cnt_wheel_viol);
    fprintf('  转向锥:     max viol = %.6f (step %d), 违反步数 = %d/100\n', ...
        max_cone_pred, step_cone_pred, cnt_cone_viol);
    fprintf('\n--- [闭环执行约束违反] 真实执行轨迹 (100 步) ---\n');
    fprintf('  轮速 SOC:   max viol = %.6f, 违反步数 = %d/100\n', ...
        max_wheel_closed, cnt_wheel_closed_viol);
    fprintf('  转向锥:     max viol = %.6f, 违反步数 = %d/100\n', ...
        max_cone_closed, cnt_cone_closed_viol);

    pred_ok = (max_wheel_pred < 1e-6) && (max_cone_pred < 1e-6);
    closed_ok = (max_wheel_closed < 1e-6) && (max_cone_closed < 1e-6);

    fprintf('\n--- [结论] ---\n');
    if pred_ok && closed_ok
        fprintf('  ALL FEASIBLE: 预测域 + 闭环域都满足原始约束 (viol < 1e-6)\n');
        fprintf('  → J_total=11.7 < 13.4 是合理的：线性化 SCP 找到更优内点 / benchmark CVX 有对偶间隙\n');
    elseif ~pred_ok
        fprintf('  PREDICTION INFEASIBLE! 预测域解了过度松弛的问题\n');
        fprintf('  → 根本原因：SOC/锥线性化是外近似 + SCP 3 次迭代不收敛\n');
        if max_wheel_pred > max_cone_pred
            fprintf('  → 主导：轮速 SOC (max viol %.4f step %d)。需加 SCP 迭代/投影修正\n', max_wheel_pred, step_wheel_pred);
        else
            fprintf('  → 主导：转向锥 (max viol %.4f step %d)。需检查线性化符号/方向\n', max_cone_pred, step_cone_pred);
        end
    else
        fprintf('  PRED FEASIBLE BUT CLOSED LOOP SMALL VIOL: 数值精度/截断效应\n');
    end
    fprintf('================================================================\n');
end
