function [wheel_speed, wheel_angle] = computeWheelOutputs(velocity, params)
% COMPUTEWHEELOUTPUTS 计算各轮的输出速度和转向角
%
% 输入：
%   velocity : 3x1 车体坐标系速度 [vx, vy, omega]
%   params   : 配置参数（含 wheel_pos）
%
% 输出：
%   wheel_speed : Nx1 各轮线速度 (N = num_wheels)
%   wheel_angle : Nx1 各轮转向角 phi (论文定义)
%
% 修复记录:
%   - [修复] 轮角定义与子模块 main.m 一致: phi = atan2(vy, vx)
%     其中 vx, vy 是轮速在车体坐标系中的分量
%     与 RSS_active_set/RSS_fmincon/RSS_proposed 子模块 main.m 完全一致
%   - 支持任意轮数，不再写死4

    if nargin < 2
        params = defaultConfig();
    end

    num_wheels = size(params.wheel_pos, 1);
    wheel_speed = zeros(num_wheels, 1);
    wheel_angle = zeros(num_wheels, 1);

    for i = 1:num_wheels
        Hj = [1, 0, -params.wheel_pos(i, 2); 0, 1, params.wheel_pos(i, 1)];
        wheel_velocity = Hj * velocity;
        wheel_speed(i) = sqrt(wheel_velocity(1)^2 + wheel_velocity(2)^2);
        % 与子模块 main.m 一致: phi = atan2(vy, vx)
        wheel_angle(i) = atan2(wheel_velocity(2), wheel_velocity(1));
    end
end
