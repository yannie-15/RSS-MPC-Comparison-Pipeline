function main()
    % ================= 清零 =================
    clear control
    clear static_counter
    clear all
    params = config();
    % ================= 获取初始轨迹 =================
    [path] = bezier_path(params.ctrl_pts, params.num_path_pts);
    % ================= 初始化=================
    q = [0.050, 0.1];         % 初始位置(m)(1*2)
    psi0 = 0.2;               % 初始朝向(rad)
    vx = 0.01;
    vy = 0.01;
    omega_b = 0.01;
    last_vel = [vx ; vy ; omega_b];    %3*1
    state = [q , psi0];               %1*3
    % ================= 定义历史记录 =================
    q_history = zeros(params.num_steps, 2);
    vi_history = zeros(params.num_steps, 4);
    phidot_history = zeros(params.num_steps, 4);
    psi_history = zeros(params.num_steps, 1);
    t_history = zeros(params.num_steps, 1);
    phi_history = zeros(params.num_steps, 4);
    nmpc_time_history = zeros(params.num_steps, 1);
    iter_num_history = zeros(params.num_steps, 1);
    
    figure('Position', [100, 100, 1500, 300]);
    
    subplot(1, 3, 1);
    plot(path(1,:), path(2,:), '-','Color',slanCL(1148,4),'LineWidth', 2, 'DisplayName', 'Reference Trajectory'); hold on;
    h_q = animatedline('Marker', 'x', 'LineStyle', 'none', 'Color',slanCL(1148,2), 'LineWidth', 1.5, 'DisplayName', 'Simulation Results');
    xlabel('x_w (m)','FontName','Times New Roman');
    ylabel('y_w (m)','FontName','Times New Roman');
    title('Trajectory Tracking Performance','FontName','Times New Roman');
    legend('FontSize', 6,'FontName','Times New Roman');
    grid on; axis equal; set(gca,'FontName','Times New Roman');
    
    subplot(1, 3, 2); hold on;
    yline(params.phidotmax, '--', 'Color', slanCL(1148,5), 'LineWidth', 2, 'DisplayName', 'Constraints');
    yline(-params.phidotmax, '--', 'Color', slanCL(1148,5), 'LineWidth', 2, 'HandleVisibility', 'off');
    h_phidot = gobjects(4,1);
    for idx = 1:4
        h_phidot(idx) = animatedline('Color', slanCL(1148,idx), 'LineWidth', 1.7, 'DisplayName', ['Wheel ' num2str(idx)]);
    end
    xlabel('Time(s)','FontName','Times New Roman'); 
    ylabel('Steering rate (rad/s)','FontName','Times New Roman'); 
    title('Steering Rate of Each Wheel','FontName','Times New Roman');
    legend('FontSize', 6,'FontName','Times New Roman'); 
    grid on; set(gca,'FontName','Times New Roman');
    ylim([-1.5 * params.phidotmax, 1.5 * params.phidotmax]);
    subplot(1, 3, 3); hold on;
    yline(params.vimax, '--', 'Color', slanCL(1148,5), 'LineWidth', 2, 'DisplayName', 'Constraints');
    h_vi = gobjects(4,1);
    for idx = 1:4
        h_vi(idx) = animatedline('Color', slanCL(1148,idx), 'LineWidth', 1.7, 'DisplayName', ['Wheel ' num2str(idx)]);
    end
    xlabel('Time (s)','FontName','Times New Roman'); 
    ylabel('Output Velocity (m/s)','FontName','Times New Roman'); 
    title('Output Velocity','FontName','Times New Roman');
    legend('FontSize', 6,'FontName','Times New Roman');
    grid on; set(gca,'FontName','Times New Roman');
    
    prev_phi = zeros(4, 1);
    
    % ======================================= 仿真循环 =========================================
    for k = 1:params.num_steps
        t = (k - 1) * params.dt;
        t_history(k) = t;
        q_history(k, :) = q;
        psi_history(k) = psi0;
        %================= 耗时统计 ================
        [new_state_dot, velocity, solve_time, iter_num] = control_RSS(path, k, last_vel, state);
        nmpc_time_history(k) = solve_time; % 直接赋值fmincon纯求解耗时
        iter_num_history(k) = iter_num; % 记录当前步迭代次数
        fprintf('【Step %d】NMPC fmincon纯求解耗时: %.4fs, 迭代次数: %d\n', k, solve_time, iter_num); 
        last_vel = velocity;
        vx = new_state_dot(1);
        vy = new_state_dot(2);
        omega_b= new_state_dot(3);
        phi = [0; 0; 0; 0];
        vi = [0; 0; 0; 0];
        %=====================计算轮子角速度========================
        for i = 1:4
            Hj = [1, 0, -params.wheel_pos(i,2);0, 1, params.wheel_pos(i,1)];
            zn = Hj * velocity;
            vxi = zn(1);
            vyi = zn(2);
            vi(i) = sqrt(vxi^2 + vyi^2);
            phi(i) = atan2(vyi, vxi);
        end
        %==========================记录=============================
        vi_history(k,:) = vi;
        phi_history(k,:) = phi;
        
        addpoints(h_q, q(1), q(2));
        for idx = 1:4
            addpoints(h_vi(idx), t, vi(idx));
        end
        if k > 1
            for idx = 1:4
                d_phi = phi(idx) - prev_phi(idx);
                d_phi = mod(d_phi + pi, 2*pi) - pi;
                curr_phidot = d_phi / params.dt;
                addpoints(h_phidot(idx), t, curr_phidot);
            end
        end
        drawnow limitrate;
        prev_phi = phi;
        
        %=========================更新状态========================
        psi0 = psi0 + omega_b * params.dt;
        q = q + [vx , vy]* params.dt;
        state_dot = new_state_dot;
        state = [q , psi0];  
    end%仿真结束
    % 角度归一化：将phidot_history每个元素限制在(-π, π]区间
    phidot_history = diff(phi_history);
    for m = 1:params.num_steps-1
        for n = 1:4
            phi_val = phidot_history(m,n);
            phi_val = mod(phi_val + pi, 2 * pi) - pi;
            phidot_history(m,n) = phi_val / params.dt;
        end
    end
    phidot_history = [zeros(1,4); phidot_history];
    % 绘图
    % plot_results(q_history, vi_history, phidot_history, psi_history, t_history, path);
    % ================= NMPC优化耗时统计 =================
    fprintf('\n=====================================================\n');
    fprintf('=============== NMPC优化耗时及迭代次数统计 ===============\n');
    fprintf('=====================================================\n');
    valid_time = nmpc_time_history(nmpc_time_history > 1e-6); % 过滤无效耗时
    valid_iter = iter_num_history(iter_num_history > 0); % 过滤无效迭代次数
    if ~isempty(valid_time)
        % 打印耗时
        total_time = sum(valid_time);    % 总耗时
        avg_time = mean(valid_time);     % 平均耗时
        max_time = max(valid_time);      % 最长单步耗时
        min_time = min(valid_time);      % 最短单步耗时
        valid_steps = length(valid_time);% 有效优化步数
        fprintf('有效NMPC优化步数: %d 步\n', valid_steps);
        fprintf('所有优化步总耗时: %.4f s\n', total_time);
        fprintf('单步优化平均耗时: %.4f s\n', avg_time);
        fprintf('单步优化最长耗时: %.4f s（%.2f倍平均）\n', max_time, max_time/avg_time);
        fprintf('单步优化最短耗时: %.4f s（%.2f倍平均）\n', min_time, min_time/avg_time);
        % 打印迭代次数
        avg_iter = mean(valid_iter);     % 平均迭代次数
        max_iter = max(valid_iter);      % 最大迭代次数
        min_iter = min(valid_iter);      % 最小迭代次数
        fprintf('单步fmincon平均迭代次数: %.1f 次\n', avg_iter);
        fprintf('单步fmincon最大迭代次数: %d 次\n', max_iter);
        fprintf('单步fmincon最小迭代次数: %d 次\n', min_iter);
    else
        fprintf('无有效NMPC优化步\n');
        fprintf('无有效NMPC迭代次数记录\n');
    end
    fprintf('=====================================================\n');
     sum1=0;
    for i=1:params.num_steps
        sum1 = sum1 + norm(q_history(i) - path(1:2, i))^2;
    end
    fprintf('RMSE是: %f\n', sqrt(sum1/params.num_steps));
end