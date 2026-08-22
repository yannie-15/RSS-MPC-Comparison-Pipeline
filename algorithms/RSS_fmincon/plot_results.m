function plot_results(q_history, vi_history, phidot_history, psi_history, t_history,path)
   params = config();
    % 调整Figure尺寸：1×3布局需更宽的窗口（原1200→1500，避免子图拥挤）
    figure('Position', [100, 100, 1500, 300]);
    path_x=path(1,:);
    path_y=path(2,:);
    
    % Fig5：运动轨迹（1×3的第1个位置，参数正确，保留）
    subplot(1, 3, 1); % 核心：1行3列，第1个图
    plot(path_x, path_y, '-','Color',slanCL(1148,4),'LineWidth', 2, 'DisplayName', 'Reference Trajectory');
    hold on; plot(q_history(:,1), q_history(:,2), 'x','Color',slanCL(1148,2), 'LineWidth', 1.5, 'DisplayName', 'Simulation Results');
    % 【修改1】坐标轴标签+字体
    xlabel('x_w (m)','FontName','Times New Roman');
    ylabel('y_w (m)','FontName','Times New Roman');
    % 【修改2】标题+字体
    title('Trajectory Tracking Performance','FontName','Times New Roman');
    % 【修改3】图例+字体
    legend('FontSize', 6,'FontName','Times New Roman');
    grid on;
    axis equal
    % 【修改4】刻度数字+字体
    set(gca,'FontName','Times New Roman');

    % Fig6：车轮转向速度（核心修改：1×3的第2个位置）
    subplot(1, 3, 2); % 原2,3,2 → 改为1,3,2
    hold on; 
    plot(t_history, phidot_history(:,1), '-', 'Color',slanCL(1148,1), 'LineWidth', 1.7, 'DisplayName', 'Wheel 1');
    plot(t_history, phidot_history(:,2), '-', 'Color',slanCL(1148,2), 'LineWidth', 1.7, 'DisplayName', 'Wheel 2');
    plot(t_history, phidot_history(:,3), '-', 'Color',slanCL(1148,3), 'LineWidth', 1.7, 'DisplayName', 'Wheel 3');
    plot(t_history, phidot_history(:,4), '-', 'Color',slanCL(1148,4), 'LineWidth', 1.7, 'DisplayName', 'Wheel 4');
    plot([0, max(t_history)], [params.phidotmax, params.phidotmax], '--','Color',slanCL(1148,5), 'LineWidth', 2, 'DisplayName', 'Constraints');
    plot([0, max(t_history)], [-params.phidotmax, -params.phidotmax], '--', 'Color',slanCL(1148,5),'LineWidth', 2, 'HandleVisibility', 'off');  % ← 原有修改保留
    hold off; 
    % 【修改5】坐标轴标签+字体
    xlabel('Time(s)','FontName','Times New Roman'); 
    ylabel('Steering rate (rad/s)','FontName','Times New Roman'); 
    % 【修改6】标题+字体
    title('Steering Rate of Each Wheel','FontName','Times New Roman');
    % 【修改7】图例+字体
    legend('FontSize', 6,'FontName','Times New Roman'); 
    grid on;
    % 【修改8】刻度数字+字体
    set(gca,'FontName','Times New Roman');

    % Fig7：车轮驱动速度（核心修改：1×3的第3个位置）
    subplot(1, 3, 3); % 原3,3,3 → 改为1,3,3
    hold on; 
    plot(t_history, vi_history(:,1), '-', 'Color',slanCL(1148,1), 'LineWidth', 1.7, 'DisplayName', 'Wheel 1');
    plot(t_history, vi_history(:,2), '-', 'Color',slanCL(1148,2), 'LineWidth', 1.7, 'DisplayName', 'Wheel 2');
    plot(t_history, vi_history(:,3), '-', 'Color',slanCL(1148,3), 'LineWidth', 1.7, 'DisplayName', 'Wheel 3');
    plot(t_history, vi_history(:,4), '-', 'Color',slanCL(1148,4), 'LineWidth', 1.7, 'DisplayName', 'Wheel 4');
    plot([0, max(t_history)], [params.vimax, params.vimax], '--','Color',slanCL(1148,5), 'LineWidth', 2, 'DisplayName', 'Constraints');
    hold off; 
    % 【修改9】坐标轴标签+字体
    xlabel('Time (s)','FontName','Times New Roman'); 
    ylabel('Output Velocity (m/s)','FontName','Times New Roman'); 
    % 【修改10】标题+字体
    title('Output Velocity','FontName','Times New Roman');
    % 【修改11】图例+字体
    legend('FontSize', 6,'FontName','Times New Roman');
    grid on;
    % 【修改12】刻度数字+字体
    set(gca,'FontName','Times New Roman');
end