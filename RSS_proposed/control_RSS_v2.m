function [u, new_state_dot, velocity] = control_RSS_v2( ...
    path, step, state_dot, state, config ...
)
% CONTROL_RSS_V2
% RSS迭代凸化控制器 - 支持多种求解器
%
% 改进点：
%  - 支持CVX求解器（原始实现）
%  - 支持HPIPM求解器（可选，需要配置）
%  - 性能测量和日志
%  - 求解器选择配置
%
% 输入：
%   path      : 3×N 参考状态
%   step      : 当前MPC步
%   state_dot : 3×1 当前车体坐标系速度
%   state     : 1×3 当前世界坐标系状态
%   config    : 配置结构体，包含:
%               - use_hpipm: 是否尝试使用HPIPM (默认false)
%               - solver_timeout: 求解器超时时间(秒)
%
% 输出：
%   u, new_state_dot, velocity: 同原始control_RSS

    if nargin < 5 || isempty(config)
        config = config_v2();
    end
    
    % 检查是否使用HPIPM
    use_hpipm = false;
    if isfield(config, 'use_hpipm')
        use_hpipm = config.use_hpipm;
    end

    % 检查是否 override 外层迭代次数 (用于 proposed-1iter)
    max_iter = 3;
    if isfield(config, 'max_iter_override') && ~isempty(config.max_iter_override)
        max_iter = config.max_iter_override;
    end
    
    % 全局变量用于记录求解时间
    global solver_time_array;
    if isempty(solver_time_array)
        solver_time_array = [];
    end
    
    %% =====================================================
    % 参数提取
    % ======================================================
    
    K = 6;
    rho = 0.01;
    k1 = 1;
    
    current_xy = [state(1); state(2)];
    psi0 = state(3);
    current_nu = state_dot;
    
    R_psi0 = [
        cos(psi0), -sin(psi0);
        sin(psi0),  cos(psi0)
    ];
    
    % 轮子特征矩阵
    num_wheels = size(config.wheel_pos, 1);
    H = cell(1, num_wheels);
    for n = 1:num_wheels
        H{n} = [
            1, 0, -config.wheel_pos(n, 2);
            0, 1,  config.wheel_pos(n, 1)
        ];
    end
    
    %% =====================================================
    % RSS外层迭代初始化
    % ======================================================

    u_hat = zeros(3, K);
    
    u = zeros(3, K);
    nu = zeros(3, K);
    
    for m = 1:max_iter
        
        fprintf('\n第 %d 个MPC步，第 %d 次RSS迭代\n', step, m);
        
        %% -------------------------------------------------
        % 选择求解器求解
        % --------------------------------------------------
        
        solver_tic = tic;
        
        try
            if use_hpipm
                % 尝试使用HPIPM
                [u, nu, status] = solve_subproblem_hpipm(...
                    path, step, current_xy, psi0, current_nu, R_psi0, ...
                    H, num_wheels, u_hat, config, K ...
                );
                solver_name = 'HPIPM';
            else
                % 使用CVX求解
                [u, nu, status] = solve_subproblem_cvx(...
                    path, step, current_xy, psi0, current_nu, R_psi0, ...
                    H, num_wheels, u_hat, config, K ...
                );
                solver_name = 'CVX';
            end
        catch ME
            % 如果HPIPM失败，回退到CVX
            if use_hpipm
                fprintf('HPIPM求解失败，回退到CVX: %s\n', ME.message);
                fprintf('完整错误: %s\n', getReport(ME));
                [u, nu, status] = solve_subproblem_cvx(...
                    path, step, current_xy, psi0, current_nu, R_psi0, ...
                    H, num_wheels, u_hat, config, K ...
                );
                solver_name = 'CVX (回退)';
            else
                rethrow(ME);
            end
        end
        
        solve_time = toc(solver_tic);
        time_index = max_iter * (step - 1) + m;
        solver_time_array(time_index) = solve_time;
        
        fprintf('求解器: %s\n', solver_name);
        fprintf('求解时间: %.6f 秒\n', solve_time);
        fprintf('状态: %s\n', status);
        
        % 检查求解是否成功 — 使用白名单判断（不能把 Infeasible/Unbounded 当成功）
        solver_success = strcmpi(status, 'Solved') || strcmpi(status, 'Inaccurate/Solved');
        
        if ~solver_success
            error('control_RSS_v2:SolverFailed', ...
                '求解失败在步 %d 迭代 %d', step, m);
        end
        
        if any(~isfinite(u(:))) || any(~isfinite(nu(:)))
            error('control_RSS_v2:InvalidSolution', ...
                '返回了NaN或Inf值');
        end
        
        %% -------------------------------------------------
        % RSS外层更新
        % --------------------------------------------------
        
        u_hat = u;
        
    end
    
    %% =====================================================
    % 输出结果
    % ======================================================
    
    rotation_world_from_body = [
        cos(state(3)), -sin(state(3)), 0;
        sin(state(3)),  cos(state(3)), 0;
        0,              0,             1
    ];
    
    velocity = current_nu + u(:, 1);
    new_state_dot = rotation_world_from_body * velocity;
    
end


%% =========================================================
% CVX求解器实现
% ==========================================================

function [u, nu, status] = solve_subproblem_cvx(...
    path, step, current_xy, psi0, current_nu, R_psi0, ...
    H, num_wheels, u_hat, config, K ...
)
% 使用CVX + SDPT3求解凸子问题
% （这部分是原始control_RSS中的求解代码）

    cvx_begin
        cvx_solver SDPT3

        variable u(3, K)
        variable nu(3, K)

        expression NU(2, K)
        expression psi(K)
        expression J

        % 预测状态序列
        NU(:, 1) = current_nu(1:2) * config.dt;
        psi(1) = psi0 + current_nu(3) * config.dt;

        for k = 1:K-1
            NU(:, k + 1) = NU(:, k) + nu(1:2, k) * config.dt;
            psi(k + 1) = psi(k) + nu(3, k) * config.dt;
        end

        % 代价函数
        J = 0;

        for k = 2:K
            reference_index = min(size(path, 2), step + k);
            position_error = current_xy - path(1:2, reference_index) + R_psi0 * NU(:, k);
            J = J + 30 * sum_square(position_error);
        end

        for k = 1:K
            reference_index = min(size(path, 2), step + k);
            orientation_error = psi(k) - path(3, reference_index);
            J = J + 1 * sum_square(orientation_error);
        end

        minimize(J + 0.3 * sum_square(u(:)) + 0.01 * sum_square(u(:) - u_hat(:)))

        subject to

            % 动力学递推
            nu(:, 1) == current_nu + u(:, 1);
            for k = 1:K-1
                nu(:, k + 1) == nu(:, k) + u(:, k + 1);
            end

            % 最大轮速约束
            for k = 1:K
                for n = 1:num_wheels
                    norm(H{n} * nu(:, k), 2) <= config.vimax;
                end
            end

            % 构造 nu_hat 序列
            delta_theta = config.dt * config.phidotmax;

            nu_hat = zeros(3, K);
            nu_hat(:, 1) = current_nu + u_hat(:, 1);
            for k = 1:K-1
                nu_hat(:, k + 1) = nu_hat(:, k) + u_hat(:, k + 1);
            end

            % 第一组转向锥凸化约束
            R = [
                sin(delta_theta), -cos(delta_theta);
                cos(delta_theta),  sin(delta_theta)
            ];

            for k = 1:K
                for n = 1:num_wheels

                    if k > 1

                        linearization_vector = ...
                            (eye(2) + R) * H{n} * nu_hat(:, k - 1) ...
                            + R * H{n} * u_hat(:, k);

                        LH = 0;
                        for l = 1:k-1
                            LH = LH ...
                                + 2 * linearization_vector' ...
                                * (eye(2) + R) * H{n} ...
                                * (u(:, l) - u_hat(:, l));
                        end

                        LH = LH ...
                            + 2 * linearization_vector' ...
                            * R * H{n} ...
                            * (u(:, k) - u_hat(:, k));

                        sum_square(H{n} * nu(:, k - 1)) ...
                            + sum_square(H{n} * (nu(:, k - 1) + u(:, k))) ...
                            - sum_square( ...
                                (eye(2) + R) * H{n} * nu_hat(:, k - 1) ...
                                + R * H{n} * u_hat(:, k) ...
                            ) ...
                            - LH <= 0;

                    else

                        linearization_vector = ...
                            (eye(2) + R) * H{n} * current_nu ...
                            + R * H{n} * u_hat(:, 1);

                        sum_square(H{n} * current_nu) ...
                            + sum_square(H{n} * (current_nu + u(:, 1))) ...
                            - sum_square(linearization_vector) ...
                            - 2 * linearization_vector' * R * H{n} ...
                            * (u(:, 1) - u_hat(:, 1)) <= 0;

                    end

                end
            end

            % 第二组转向锥凸化约束
            R = [
                sin(delta_theta),  cos(delta_theta);
               -cos(delta_theta),  sin(delta_theta)
            ];

            for k = 1:K
                for n = 1:num_wheels

                    if k > 1

                        linearization_vector = ...
                            (eye(2) + R) * H{n} * nu_hat(:, k - 1) ...
                            + R * H{n} * u_hat(:, k);

                        LH = 0;
                        for l = 1:k-1
                            LH = LH ...
                                + 2 * linearization_vector' ...
                                * (eye(2) + R) * H{n} ...
                                * (u(:, l) - u_hat(:, l));
                        end

                        LH = LH ...
                            + 2 * linearization_vector' ...
                            * R * H{n} ...
                            * (u(:, k) - u_hat(:, k));

                        sum_square(H{n} * nu(:, k - 1)) ...
                            + sum_square(H{n} * (nu(:, k - 1) + u(:, k))) ...
                            - sum_square( ...
                                (eye(2) + R) * H{n} * nu_hat(:, k - 1) ...
                                + R * H{n} * u_hat(:, k) ...
                            ) ...
                            - LH <= 0;

                    else

                        linearization_vector = ...
                            (eye(2) + R) * H{n} * current_nu ...
                            + R * H{n} * u_hat(:, 1);

                        sum_square(H{n} * current_nu) ...
                            + sum_square(H{n} * (current_nu + u(:, 1))) ...
                            - sum_square(linearization_vector) ...
                            - 2 * linearization_vector' * R * H{n} ...
                            * (u(:, 1) - u_hat(:, 1)) <= 0;

                    end

                end
            end

    cvx_end
    
    status = cvx_status;
    
end


%% =========================================================
% HPIPM求解器实现
% ==========================================================

function [u, nu, status] = solve_subproblem_hpipm(...
    path, step, current_xy, psi0, current_nu, R_psi0, ...
    H, num_wheels, u_hat, config, K ...
)
% 使用HPIPM求解凸子问题
% 调用Python接口

    % state 格式: 1x3 [x, y, psi]（与 construct_complete_qp_from_rss 的输入一致）
    state = [current_xy(1), current_xy(2), psi0];

    % 构造QCQP问题
    qp_problem = construct_complete_qp_from_rss(...
        path, step, current_nu, state, u_hat, config ...
    );

    % 调用Python HPIPM求解器（包含QCQP约束）
    try
        solution = solve_qp_with_python_hpipm(...
            qp_problem.H, qp_problem.g, ...
            qp_problem.A, qp_problem.b, ...
            qp_problem.C, qp_problem.d, ...
            qp_problem.lb, qp_problem.ub, ...
            qp_problem.Hq, qp_problem.gq, qp_problem.uq ...
        );
    catch
        % 如果Python失败，回退到CVX
        error('HPIPM:SolverFailed', 'Failed to call Python HPIPM solver');
    end

    % 检查求解状态
    status = solution.status;
    if ~ischar(status)
        status = num2str(status);
    end

    % HPIPM 求解失败时，让外层 catch 回退到 CVX
    if startsWith(status, 'Failed')
        error('HPIPM:SolverFailed', 'HPIPM status: %s', status);
    end

    % 提取解
    x = solution.x;

    % 检查解的有效性
    if any(~isfinite(x))
        error('HPIPM:InvalidSolution', 'HPIPM returned NaN/Inf');
    end

    u = reshape(x(1:3*K), 3, K);
    nu = reshape(x(3*K+1:end), 3, K);

end


%% =========================================================
% 配置函数（v2版本）
% ==========================================================

function params = config_v2()
    % 基础配置
    params = config();  % 调用原始config函数
    
    % v2版本的额外参数
    params.use_hpipm = false;  % 默认不使用HPIPM
    params.solver_timeout = 10;  % 求解器超时时间
    params.enable_solver_logging = true;  % 记录求解器信息
    
end
