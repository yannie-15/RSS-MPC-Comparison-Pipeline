function main()
% MAIN  兼容入口 - 调用 run_closed_loop 完成仿真
%
% 原始 main.m 逻辑已迁移到 run_closed_loop.m。
% 本文件保持向后兼容, 用法不变:
%   >> cd D:\PROJECT\RSS_V2
%   >> main
%
% 如需自定义参数, 直接调用:
%   cfg = config();
%   cfg.live_plot = true;
%   result = run_closed_loop(cfg);

    % 加载默认配置 (从 algorithms/RSS_proposed/config.m)
    cfg = config();

    % 兼容原 main.m 行为: 开实时绘图
    cfg.live_plot = true;
    cfg.save_full_log = false;
    cfg.case_id = 'main_compat';
    cfg.solver = 'sdpt3';

    % 运行闭环仿真
    result = run_closed_loop(cfg);

    % 打印摘要 (与原 main.m 输出格式一致)
    fprintf('\nRMSE是: %f\n', result.position_rmse);
    fprintf('总用时是：%f\n', result.total_solve_time);
    fprintf('代价是：%f\n', result.trajectory_cost);
end
