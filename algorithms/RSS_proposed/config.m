function params = config()
    % ====================== 机器人几何参数 ======================
    params.Lx = 0.655;          % 车身长度(mm)
    params.Ly = 0.335;          % 车身宽度(mm)
    params.a = params.Lx/2;   % 车轮x向偏移(mm)
    params.b = params.Ly/2;   % 车轮y向偏移(mm)
    params.wheel_pos = [params.a, params.b; -params.a, params.b; -params.a, -params.b; params.a, -params.b]; % 4车轮位置
    params.vimax = 5;     % 最大驱动速度(mm/s)
    params.phidotmax = 5 * pi;   % 最大转向速度(rad/s)

    % ====================== 控制器参数 ======================
    params.k1 = 0.15;
    params.k2 = 0.15;
    params.k3 = 0.1;
    params.eps = 0.001;
    % ====================== 仿真参数 ======================
    params.dt = 0.01;         % 时间步长(s)
    params.t_end = 1;       % 仿真时长(s)
    params.num_steps = round(params.t_end/params.dt);

    % ====================== 期望路径参数 ======================
    params.ctrl_pts = [0, 0;
                       0.750, 0.250;
                       0.250, 0.750;
                       1.250, -1.000;
                       1.000, 0]; % 贝塞尔曲线控制点
    params.num_path_pts = params.num_steps; % 路径离散点数
    %===================================================================
    
end