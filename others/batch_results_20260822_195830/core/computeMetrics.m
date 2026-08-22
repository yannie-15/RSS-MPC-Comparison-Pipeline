function metrics = computeMetrics(path, states, params, solverTimes, wheelSpeeds, wheelAngles, executedU)
% COMPUTEMETRICS 计算仿真结果的完整指标
%
% 输入：
%   path         : 3xN 参考轨迹
%   states       : 3x(M+1) 仿真状态序列 (M步 + 1个最终状态)
%   params       : 配置参数
%   solverTimes  : 1xM 每MPC步求解时间（可选）
%   wheelSpeeds  : NxM 轮速（可选, N=num_wheels）
%   wheelAngles  : NxM 轮转角（可选）
%   executedU    : 3xM 每步执行的控制增量（可选, 用于轨迹代价）
%
% 输出：
%   metrics : 结构体，包含所有评价指标
%
% 修复记录:
%   - [P1] 新增式(21)轨迹代价 trajectoryCost
%   - [P2] 状态/参考索引对齐修正

    if nargin < 3 || isempty(params)
        params = defaultConfig();
    end

    % states(:,1) 是初始状态, states(:,k+1) 是第k步控制后的状态
    % 对应参考 path(:,k) (k从1开始)
    % 所以 states(:,k) 对应 path(:,k) (控制前的状态与参考对齐)
    num_steps = size(states, 2) - 1;  % 最后一列是最终状态
    num_ref = size(path, 2);

    % 处理空仿真 (控制器第一步就失败)
    if num_steps == 0
        metrics.rmse = NaN;
        metrics.maxPositionError = NaN;
        metrics.meanPositionError = NaN;
        metrics.finalPositionError = NaN;
        metrics.rmseOrientation = NaN;
        metrics.maxOrientationError = NaN;
        metrics.meanSolveTime = NaN;
        metrics.maxSolveTime = NaN;
        metrics.totalSolveTime = NaN;
        metrics.wheelSpeedViolationRatio = NaN;
        metrics.maxWheelSpeedViolation = NaN;
        metrics.steeringRateViolationRatio = NaN;
        metrics.maxSteeringRateViolation = NaN;
        metrics.meanWheelSpeedChange = NaN;
        metrics.maxWheelSpeedChange = NaN;
        metrics.pathLength = sum(vecnorm(diff(path(1:2, :), 1, 2), 2, 1));
        metrics.numSteps = 0;
        metrics.trajectoryCost = NaN;
        return;
    end

    %% =====================================================
    % 1. 位置跟踪精度
    % ======================================================

    position_errors = zeros(num_steps, 1);
    for i = 1:num_steps
        ref_idx = min(i, num_ref);
        % states(:,i) 是第i步控制前的状态, 对应参考 path(:,i)
        % [P2修复] 使用 states(:,i) 而不是 states(:,i+1)
        position_errors(i) = norm(states(1:2, i) - path(1:2, ref_idx));
    end

    metrics.rmse = sqrt(mean(position_errors.^2));
    metrics.maxPositionError = max(position_errors);
    metrics.meanPositionError = mean(position_errors);
    metrics.finalPositionError = position_errors(end);

    %% =====================================================
    % 2. 姿态跟踪精度
    % ======================================================

    orientation_errors = zeros(num_steps, 1);
    for i = 1:num_steps
        ref_idx = min(i, num_ref);
        orientation_errors(i) = abs(wrap_angle(states(3, i) - path(3, ref_idx)));
    end

    metrics.rmseOrientation = sqrt(mean(orientation_errors.^2));
    metrics.maxOrientationError = max(orientation_errors);

    %% =====================================================
    % 3. 求解时间
    % ======================================================

    if nargin >= 4 && ~isempty(solverTimes) && any(solverTimes > 0)
        valid_times = solverTimes(solverTimes > 0);
        metrics.meanSolveTime = mean(valid_times);
        metrics.maxSolveTime = max(valid_times);
        metrics.totalSolveTime = sum(valid_times);
    else
        metrics.meanSolveTime = NaN;
        metrics.maxSolveTime = NaN;
        metrics.totalSolveTime = NaN;
    end

    %% =====================================================
    % 4. 约束违反率
    % ======================================================

    if nargin >= 6 && ~isempty(wheelSpeeds)
        vimax = params.vimax;
        speed_violations = sum(wheelSpeeds > vimax * (1 + 1e-6), 'all');
        metrics.wheelSpeedViolationRatio = speed_violations / numel(wheelSpeeds);
        metrics.maxWheelSpeedViolation = max(0, max(wheelSpeeds(:)) - vimax);
    else
        metrics.wheelSpeedViolationRatio = NaN;
        metrics.maxWheelSpeedViolation = NaN;
    end

    if nargin >= 6 && ~isempty(wheelAngles) && size(wheelAngles, 2) > 1
        phidotmax = params.phidotmax;
        d_angles = diff(wheelAngles, 1, 2);
        % 角度归一化到 (-pi, pi]
        d_angles = mod(d_angles + pi, 2*pi) - pi;
        phidot = d_angles / params.dt;
        phidot_violations = sum(abs(phidot) > phidotmax * (1 + 1e-6), 'all');
        metrics.steeringRateViolationRatio = phidot_violations / numel(phidot);
        metrics.maxSteeringRateViolation = max(0, max(abs(phidot(:))) - phidotmax);
    else
        metrics.steeringRateViolationRatio = NaN;
        metrics.maxSteeringRateViolation = NaN;
    end

    %% =====================================================
    % 5. 控制平滑度
    % ======================================================

    if nargin >= 5 && ~isempty(wheelSpeeds)
        if size(wheelSpeeds, 2) > 1
            speed_diff = diff(wheelSpeeds, 1, 2);
            metrics.meanWheelSpeedChange = mean(abs(speed_diff(:)));
            metrics.maxWheelSpeedChange = max(abs(speed_diff(:)));
        else
            metrics.meanWheelSpeedChange = NaN;
            metrics.maxWheelSpeedChange = NaN;
        end
    else
        metrics.meanWheelSpeedChange = NaN;
        metrics.maxWheelSpeedChange = NaN;
    end

    %% =====================================================
    % 6. 路径信息
    % ======================================================

    metrics.pathLength = sum(vecnorm(diff(path(1:2, :), 1, 2), 2, 1));
    metrics.numSteps = num_steps;

    %% =====================================================
    % 7. [P1新增] 式(21)轨迹代价 J_total
    % ======================================================
    % J_total = sum_{k=1}^{M} [ (e^k)'*Q*e^k + (u^k)'*R*u^k ]
    % 其中 e^k = states(:,k) - path(:,k) 为跟踪误差
    %       u^k = executedU(:,k) 为执行的控制增量
    %       Q = diag([w_pos, w_pos, w_psi]) = diag([30, 30, 1])
    %       R = diag([w_control, w_control, w_control]) = diag([0.3, 0.3, 0.3])

    if nargin >= 7 && ~isempty(executedU)
        Q_weights = [30, 30, 1];      % 位置xy权重 + 姿态权重
        R_weights = [0.3, 0.3, 0.3];  % 控制增量权重

        trajectoryCost = 0;
        for k = 1:num_steps
            ref_idx = min(k, num_ref);
            e = states(:, k) - path(:, ref_idx);
            e(3) = wrap_angle(e(3));  % 姿态误差归一化

            u_k = executedU(:, k);

            trajectoryCost = trajectoryCost ...
                + e' * diag(Q_weights) * e ...
                + u_k' * diag(R_weights) * u_k;
        end

        metrics.trajectoryCost = trajectoryCost;
    else
        metrics.trajectoryCost = NaN;
    end

end


function a = wrap_angle(a)
% 将角度归一化到 (-pi, pi]
    a = mod(a + pi, 2*pi) - pi;
end
