function report = verify_constraints_RSS()
% VERIFY_CONSTRAINTS_RSS  验证 RSS proposed 算法在 SQP 每次子迭代中是否满足原始约束
%
% 检查两类约束：
%   1. 轮速约束 (凸, CVX直接强制):  ||H_n * nu_k|| <= vimax
%   2. 转向角速率约束 (非凸, SQP线性化):  相邻步轮角变化 <= delta_theta
%
% 关键: 第2类约束在CVX中是被线性化后强制执行的, 线性化近似可能不满足原始约束。
%       本脚本在每次CVX求解后, 用原始非线性公式重新检验。
%
% 用法:
%   cd('D:\PROJECT\RSS_V2\matlab')
%   setup_paths;
%   verify_constraints_RSS

    %% ========== 路径设置 ==========
    % 确保 RSS_proposed 的 config() 和 bezier_path() 可用
    this_dir = fileparts(mfilename('fullpath'));
    rss_proposed_dir = fullfile(this_dir, '..', 'algorithms', 'RSS_proposed');
    if exist(rss_proposed_dir, 'dir')
        addpath(rss_proposed_dir);
    end

    %% ========== 初始化 ==========
    params = config();
    [path] = bezier_path(params.ctrl_pts, params.num_path_pts);

    % 初始状态 (与 main.m 一致)
    q = [0.050, 0.1];
    psi0 = 0.2;
    vx = 0.01; vy = 0.01; omega_b = 0.01;
    last_vel = [vx; vy; omega_b];
    state = [q, psi0];

    % 算法参数 (与 control_RSS.m 一致)
    K = 6;  rho = 0.01; k1 = 1;
    max_iter = 3;

    % H 矩阵
    H = cell(1, 4);
    for n = 1:4
        H{n} = [1, 0, -params.wheel_pos(n, 2); 0, 1, params.wheel_pos(n, 1)];
    end

    % 约束参数
    vimax = params.vimax;
    delta_theta = params.dt * params.phidotmax;  % 单步最大转角变化

    % 结果存储
    num_steps = params.num_steps;
    % 每步每次子迭代的违反记录
    viol.speed_max = zeros(num_steps, max_iter);       % 每步每次迭代最大轮速
    viol.speed_viol = zeros(num_steps, max_iter);      % 每步每次迭代最大轮速超出量
    viol.angle_max = zeros(num_steps, max_iter);       % 每步每次迭代最大转角变化
    viol.angle_viol = zeros(num_steps, max_iter);      % 每步每次迭代最大转角变化超出量
    viol.cvx_status = cell(num_steps, max_iter);       % CVX求解状态
    viol.cvx_optval = zeros(num_steps, max_iter);      % CVX最优值
    viol.solver_failed = false(num_steps, max_iter);   % 求解器是否失败

    fprintf('========== RSS 约束验证开始 ==========\n');
    fprintf('参数: vimax=%.2f, phidotmax=%.4f, delta_theta=%.6f rad, dt=%.3f s\n', ...
        vimax, params.phidotmax, delta_theta, params.dt);
    fprintf('仿真步数: %d, SQP子迭代: %d, 预测时域: %d\n\n', num_steps, max_iter, K);

    %% ========== 闭环仿真主循环 ==========
    for step = 1:num_steps
        % --- control_RSS 内部逻辑 (带约束检查) ---
        current_xy = [state(1), state(2)]';
        psi0_val = state(3);
        current_nu = last_vel;
        R_psi0 = [cos(psi0_val), -sin(psi0_val); sin(psi0_val), cos(psi0_val)];
        u_hat = zeros(3, K);

        for m = 1:max_iter
            % CVX 求解 (与 control_RSS.m 完全一致)
            cvx_begin quiet
                cvx_solver SDPT3;
                variable u(3, K)
                variable nu(3, K)
                expression NU(2, K);
                expression psi(K);
                expression J;

                NU(:, 1) = current_nu(1:2) * params.dt;
                psi(1) = psi0_val + current_nu(3) * params.dt;
                for k = 1:K-1
                    NU(:, k + 1) = NU(:, k) + nu(1:2, k) * params.dt;
                    psi(k + 1) = psi(k) + nu(3, k) * params.dt;
                end

                J = 0;
                for k = 2:K
                    J = J + 30 * sum_square(current_xy - path(1:2, min(size(path, 2), step + k)) + R_psi0 * NU(:, k));
                end
                for k = 1:K
                    J = J + k1 * sum_square(psi(k) - path(3, min(size(path, 2), step + k)));
                end
                minimize(J + 0.3 * sum_square(u(:)) + rho * sum_square(u(:) - u_hat(:)));

                subject to
                    nu(:, 1) == current_nu + u(:, 1);
                    for k = 1:K-1
                        nu(:, k + 1) == nu(:, k) + u(:, k + 1);
                    end
                    for k = 1:K
                        for n = 1:4
                            norm(H{n} * nu(:, k), 2) <= vimax;
                        end
                    end
                    % 非凸约束线性化 (与原代码一致)
                    R_mat = [sin(delta_theta), -cos(delta_theta); cos(delta_theta), sin(delta_theta)];
                    nu_hat = zeros(3, K);
                    nu_hat(:, 1) = current_nu + u_hat(:, 1);
                    for k = 1:K-1
                        nu_hat(:, k + 1) = nu_hat(:, k) + u_hat(:, k + 1);
                    end
                    for k = 1:K
                        for n = 1:4
                            expression LH;
                            LH = 0;
                            if k > 1
                                for l = 1:k-1
                                    LH = LH + 2 * ((eye(2) + R_mat)* H{n} * nu_hat(:, k - 1) + R_mat * H{n} * u_hat(:, k))' * (eye(2) + R_mat) * H{n} * (u(:, l) - u_hat(:, l));
                                end
                                LH = LH + 2 * ((eye(2) + R_mat)* H{n} * nu_hat(:, k - 1) + R_mat * H{n} * u_hat(:, k))' * R_mat * H{n} * (u(:, k) - u_hat(:, k));
                            end
                            if k > 1
                                sum_square(H{n} * nu(:, k - 1)) + sum_square(H{n} * (nu(:, k - 1) + u(:, k)))...
                                    - sum_square((eye(2) + R_mat) * H{n} * nu_hat(:, k-1) + R_mat * H{n} * u_hat(:, k)) - LH <= 0;
                            end
                            if k == 1
                                sum_square(H{n} * current_nu) + sum_square(H{n} * (current_nu + u(:, 1)))...
                                    - sum_square((eye(2) + R_mat) * H{n} * current_nu + R_mat * H{n} * u_hat(:, 1))...
                                    - 2 * ((eye(2) + R_mat) * H{n} * current_nu + R_mat * H{n} * u_hat(:, 1))' * (R_mat * H{n} * (u(:, 1) - u_hat(:, 1))) <= 0;
                            end
                        end
                    end
                    R_mat = [sin(delta_theta), cos(delta_theta); -cos(delta_theta), sin(delta_theta)];
                    for k = 1:K
                        for n = 1:4
                            expression LH;
                            LH = 0;
                            if k > 1
                                for l = 1:k-1
                                    LH = LH + 2 * ((eye(2) + R_mat)* H{n} * nu_hat(:, k - 1) + R_mat * H{n} * u_hat(:, k))' * (eye(2) + R_mat) * H{n} * (u(:, l) - u_hat(:, l));
                                end
                                LH = LH + 2 * ((eye(2) + R_mat)* H{n} * nu_hat(:, k - 1) + R_mat * H{n} * u_hat(:, k))' * R_mat * H{n} * (u(:, k) - u_hat(:, k));
                            end
                            if k > 1
                                sum_square(H{n} * nu(:, k - 1)) + sum_square(H{n} * (nu(:, k - 1) + u(:, k)))...
                                    - sum_square((eye(2) + R_mat) * H{n} * nu_hat(:, k-1) + R_mat * H{n} * u_hat(:, k)) - LH <= 0;
                            end
                            if k == 1
                                sum_square(H{n} * current_nu) + sum_square(H{n} * (current_nu + u(:, 1)))...
                                    - sum_square((eye(2) + R_mat) * H{n} * current_nu + R_mat * H{n} * u_hat(:, 1))...
                                    - 2 * ((eye(2) + R_mat) * H{n} * current_nu + R_mat * H{n} * u_hat(:, 1))' * (R_mat * H{n} * (u(:, 1) - u_hat(:, 1))) <= 0;
                            end
                        end
                    end
            cvx_end

            % ========== 原始约束检验 ==========
            viol.cvx_status{step, m} = cvx_status;
            viol.cvx_optval(step, m) = cvx_optval;
            solver_ok = strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved');
            viol.solver_failed(step, m) = ~solver_ok;

            if ~solver_ok
                fprintf('[Step %d, Iter %d] CVX 状态: %s (未求解成功, 跳过约束检查)\n', ...
                    step, m, cvx_status);
                u_hat = u;
                continue;
            end

            % --- 检查1: 轮速约束 (凸, 应该一定满足) ---
            max_speed = 0;
            max_speed_viol = 0;
            for k = 1:K
                for n = 1:4
                    wheel_vel = H{n} * nu(:, k);
                    speed = norm(wheel_vel, 2);
                    if speed > max_speed
                        max_speed = speed;
                    end
                    if speed > vimax
                        viol_amount = speed - vimax;
                        if viol_amount > max_speed_viol
                            max_speed_viol = viol_amount;
                        end
                        fprintf('[Step %d, Iter %d] 轮速违反: k=%d, n=%d, speed=%.6f > vimax=%.2f (超出 %.6f)\n', ...
                            step, m, k, n, speed, vimax, viol_amount);
                    end
                end
            end
            viol.speed_max(step, m) = max_speed;
            viol.speed_viol(step, m) = max_speed_viol;

            % --- 检查2: 转向角速率约束 (非凸, 可能违反) ---
            % 原始约束: 相邻步的轮角变化 |theta_k - theta_{k-1}| <= delta_theta
            % 其中 theta = atan2(vy, vx) 是轮速度向量的方向角
            max_angle_change = 0;
            max_angle_viol = 0;

            for n = 1:4
                % 计算每步的轮角
                wheel_angles = zeros(1, K);
                wheel_speeds_vec = zeros(1, K);
                for k = 1:K
                    wv = H{n} * nu(:, k);
                    wheel_speeds_vec(k) = norm(wv, 2);
                    wheel_angles(k) = atan2(wv(2), wv(1));
                end

                % 当前状态 (step 0) 的轮角
                wv0 = H{n} * current_nu;
                angle_prev = atan2(wv0(2), wv0(1));
                speed_prev = norm(wv0, 2);

                % 逐步检查角变化
                for k = 1:K
                    angle_curr = wheel_angles(k);
                    speed_curr = wheel_speeds_vec(k);

                    % 只在有实质速度时检查角度变化 (避免零速度附近的数值噪声)
                    if speed_prev > 1e-6 && speed_curr > 1e-6
                        d_angle = angle_curr - angle_prev;
                        % 归一化到 (-pi, pi]
                        d_angle = mod(d_angle + pi, 2*pi) - pi;
                        abs_d_angle = abs(d_angle);

                        if abs_d_angle > max_angle_change
                            max_angle_change = abs_d_angle;
                        end

                        if abs_d_angle > delta_theta * (1 + 1e-6)
                            viol_amount = abs_d_angle - delta_theta;
                            if viol_amount > max_angle_viol
                                max_angle_viol = viol_amount;
                            end
                            fprintf('[Step %d, Iter %d] 转角违反: n=%d, k=%d, |d_angle|=%.6f > delta_theta=%.6f (超出 %.6f rad, %.4f%%)\n', ...
                                step, m, n, k, abs_d_angle, delta_theta, viol_amount, ...
                                viol_amount/delta_theta*100);
                        end
                    end

                    angle_prev = angle_curr;
                    speed_prev = speed_curr;
                end
            end
            viol.angle_max(step, m) = max_angle_change;
            viol.angle_viol(step, m) = max_angle_viol;

            % SQP 热启动更新
            u_hat = u;
        end

        % --- 状态推进 (与 main.m 一致) ---
        u_final = u_hat;  % 最后一次迭代的解
        new_state_dot = [cos(state(3)), -sin(state(3)), 0;
                         sin(state(3)),  cos(state(3)), 0;
                         0,              0, 1] * (last_vel + 1.00 * u_final(:, 1));
        velocity = current_nu + u_final(:, 1);

        vx = new_state_dot(1);
        vy = new_state_dot(2);
        omega_b = new_state_dot(3);

        psi0 = psi0 + omega_b * params.dt;
        q = q + [vx, vy] * params.dt;
        state = [q, psi0];
        last_vel = velocity;

        % 进度报告
        if mod(step, 10) == 0 || step == 1
            fprintf('Step %d/%d 完成\n', step, num_steps);
        end
    end

    %% ========== 汇总报告 ==========
    fprintf('\n========== 约束验证汇总 ==========\n\n');

    % 1. 求解器状态
    failed_count = sum(viol.solver_failed(:));
    total_solves = num_steps * max_iter;
    fprintf('1. 求解器状态:\n');
    fprintf('   总求解次数: %d\n', total_solves);
    fprintf('   失败次数: %d (%.1f%%)\n', failed_count, failed_count/total_solves*100);
    if failed_count > 0
        [fail_steps, fail_iters] = find(viol.solver_failed);
        for i = 1:length(fail_steps)
            fprintf('   -> Step %d, Iter %d: %s\n', fail_steps(i), fail_iters(i), viol.cvx_status{fail_steps(i), fail_iters(i)});
        end
    end
    fprintf('\n');

    % 2. 轮速约束
    fprintf('2. 轮速约束 (凸, CVX直接强制):\n');
    fprintf('   约束: ||H_n * nu_k|| <= vimax = %.2f\n', vimax);
    fprintf('   全局最大轮速: %.6f\n', max(viol.speed_max(:)));
    fprintf('   最大违反量: %.6e\n', max(viol.speed_viol(:)));
    speed_viol_count = sum(viol.speed_viol(:) > 1e-6);
    fprintf('   违反次数: %d / %d (%.2f%%)\n', speed_viol_count, total_solves, speed_viol_count/total_solves*100);
    if speed_viol_count > 0
        [vs, vi] = find(viol.speed_viol > 1e-6);
        for i = 1:length(vs)
            fprintf('   -> Step %d, Iter %d: 最大超出 %.6e\n', vs(i), vi(i), viol.speed_viol(vs(i), vi(i)));
        end
    end
    fprintf('\n');

    % 3. 转向角速率约束
    fprintf('3. 转向角速率约束 (非凸, SQP线性化近似):\n');
    fprintf('   约束: |d_angle| <= delta_theta = %.6f rad (= dt*phidotmax = %.3f*%.4f)\n', ...
        delta_theta, params.dt, params.phidotmax);
    fprintf('   全局最大转角变化: %.6f rad (%.2f%% of limit)\n', ...
        max(viol.angle_max(:)), max(viol.angle_max(:))/delta_theta*100);
    fprintf('   最大违反量: %.6e rad\n', max(viol.angle_viol(:)));
    angle_viol_count = sum(viol.angle_viol(:) > 1e-6);
    fprintf('   违反次数: %d / %d (%.2f%%)\n', angle_viol_count, total_solves, angle_viol_count/total_solves*100);
    if angle_viol_count > 0
        [as, ai] = find(viol.angle_viol > 1e-6);
        % 按违反量排序, 只显示前20个
        viol_values = viol.angle_viol(sub2ind(size(viol.angle_viol), as, ai));
        [sorted_vals, sort_idx] = sort(viol_values, 'descend');
        show_count = min(20, length(sort_idx));
        fprintf('   违反最严重的 %d 个 (共 %d):\n', show_count, length(sort_idx));
        for i = 1:show_count
            s = as(sort_idx(i));
            it = ai(sort_idx(i));
            fprintf('   -> Step %d, Iter %d: |d_angle|=%.6f, 超出=%.6e rad (%.2f%% of limit)\n', ...
                s, it, viol.angle_max(s, it), viol.angle_viol(s, it), viol.angle_viol(s, it)/delta_theta*100);
        end
    end
    fprintf('\n');

    % 4. 按SQP迭代次数分析
    fprintf('4. 按SQP子迭代分析 (转角约束):\n');
    for m = 1:max_iter
        m_viol = sum(viol.angle_viol(:, m) > 1e-6);
        m_max = max(viol.angle_viol(:, m));
        m_max_angle = max(viol.angle_max(:, m));
        fprintf('   Iter %d: 违反 %d/%d 步, 最大超出 %.6e rad, 最大转角变化 %.6f rad (%.1f%% of limit)\n', ...
            m, m_viol, num_steps, m_max, m_max_angle, m_max_angle/delta_theta*100);
    end
    fprintf('\n');

    % 5. 结论
    fprintf('========== 结论 ==========\n');
    if speed_viol_count == 0 && angle_viol_count == 0
        fprintf('✓ 所有 %d 次CVX求解 (100步 x 3次SQP迭代) 均满足原始约束。\n', total_solves);
        fprintf('  轮速约束: 无违反 (凸约束, CVX直接强制)\n');
        fprintf('  转角约束: 无违反 (SQP线性化后的解也满足原始非凸约束)\n');
    else
        if speed_viol_count > 0
            fprintf('✗ 轮速约束有违反: %d 次\n', speed_viol_count);
        else
            fprintf('✓ 轮速约束: 无违反\n');
        end
        if angle_viol_count > 0
            fprintf('✗ 转角约束有违反: %d 次 (%.2f%%)\n', angle_viol_count, angle_viol_count/total_solves*100);
            fprintf('  说明: SQP线性化近似的解不完全满足原始非凸约束。\n');
            fprintf('  这在SQP框架中是可能的, 关键看最终迭代结果是否收敛到可行解。\n');
        else
            fprintf('✓ 转角约束: 无违反 (SQP线性化后的解也满足原始非凸约束)\n');
        end
    end

    % 保存结果
    report = viol;
    report.params = params;
    report.num_steps = num_steps;
    report.max_iter = max_iter;
    report.vimax = vimax;
    report.delta_theta = delta_theta;
    report.total_solves = total_solves;
    report.speed_viol_count = speed_viol_count;
    report.angle_viol_count = angle_viol_count;

    % 保存到 results 目录
    results_dir = fullfile(this_dir, '..', 'results');
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    save(fullfile(results_dir, 'constraint_verification_results.mat'), 'report');
    fprintf('\n结果已保存到 results/constraint_verification_results.mat\n');
end
