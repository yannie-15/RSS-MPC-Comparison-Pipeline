function [max_viol, wheel_viol, cone_viol] = check_qk_violation(u_sol, u_anchor, v0, params)
% CHECK_QK_VIOLATION  检查 u_sol 在固定 Q_K(u_anchor) 下的精确约束违反量
%
% 输入:
%   u_sol     - 3×K 优化解 (控制增量序列)
%   u_anchor  - 3×K RSS 锚点 (固定 B/L 的 û)
%   v0        - 3×1 当前车体系速度 ν_0
%   params    - 配置参数 (dt, phidotmax, vimax, wheel_pos)
%
% 输出:
%   max_viol   - 最大违反量 (max(wheel_viol, cone_viol))
%   wheel_viol - 轮速 SOC 最大违反量
%   cone_viol  - 转向锥凸化约束最大违反量
%
% 约束定义 (固定 Q_K(u_anchor)):
%   1. 轮速 SOC (论文 20b): ||H_n * ν_k||^2 - vimax^2 <= 0
%   2. 转向锥凸化 (论文 15-16): C = A - B(u_anchor) - L(u, u_anchor) <= 0
%      A = 0.5*||H*ν_{k-1}||^2 + 0.5*||H*(ν_{k-1}+u_k)||^2
%      B = 0.5*||(I+R)*H*ν̂_{k-1} + R*H*û_k||^2
%      L = ell_anchor' * (T*(ν_{k-1} - ν̂_{k-1}) + U*(u_k - û_k))

    K = 6;
    dt = params.dt;
    phidotmax = params.phidotmax;
    vimax = params.vimax;
    wheel_pos = params.wheel_pos;
    num_wheels = size(wheel_pos, 1);

    % Hn 矩阵
    Hn = cell(1, num_wheels);
    for n = 1:num_wheels
        Hn{n} = [1, 0, -wheel_pos(n,2); 0, 1, wheel_pos(n,1)];
    end

    % delta_theta
    delta_theta = dt * phidotmax;

    % R1, R2 (论文 12)
    R1 = [sin(delta_theta), -cos(delta_theta); cos(delta_theta),  sin(delta_theta)];
    R2 = [sin(delta_theta),  cos(delta_theta); -cos(delta_theta), sin(delta_theta)];

    % ν_sol 序列 (from u_sol)
    nu_sol = zeros(3, K+1);
    nu_sol(:, 1) = v0;
    for k = 1:K
        nu_sol(:, k+1) = nu_sol(:, k) + u_sol(:, k);
    end

    % ν̂_anchor 序列 (from u_anchor)
    nu_hat_anchor = zeros(3, K+1);
    nu_hat_anchor(:, 1) = v0;
    for k = 1:K
        nu_hat_anchor(:, k+1) = nu_hat_anchor(:, k) + u_anchor(:, k);
    end

    % ---- 1. 轮速 SOC 约束 (k=1..K) ----
    wheel_viol = 0;
    for k = 1:K
        nu_k = nu_sol(:, k+1);
        for n = 1:num_wheels
            Mn = Hn{n}' * Hn{n};
            val = nu_k' * Mn * nu_k - vimax^2;  % ||H*ν_k||^2 - vimax^2
            if val > wheel_viol
                wheel_viol = val;
            end
        end
    end

    % ---- 2. 转向锥凸化约束 C = A - B - L <= 0 (k=1..K) ----
    cone_viol = 0;
    for k = 1:K
        v_km1 = nu_sol(:, k);       % ν_{k-1} (from u_sol)
        u_k = u_sol(:, k);          % u_k (from u_sol)
        v_hat_km1 = nu_hat_anchor(:, k);  % ν̂_{k-1} (from u_anchor)
        u_hat_k = u_anchor(:, k);         % û_k (from u_anchor)

        for n = 1:num_wheels
            Hi = Hn{n};
            Mn = Hi' * Hi;

            for gg = 1:2
                if gg == 1; Rg = R1; else; Rg = R2; end
                T = (eye(2) + Rg) * Hi;   % 2x3
                U_mat = Rg * Hi;           % 2x3

                % B 项 (常数, 固定 u_anchor)
                ell_anchor = T * v_hat_km1 + U_mat * u_hat_k;
                B_const = 0.5 * (ell_anchor' * ell_anchor);

                % A 项 (凸二次, 在 u_sol 处求值)
                A_val = v_km1' * Mn * v_km1 + v_km1' * Mn * u_k + 0.5 * u_k' * Mn * u_k;

                % L 项 (线性, 锚点 u_anchor, 求值点 u_sol)
                L_val = ell_anchor' * (T * (v_km1 - v_hat_km1) + U_mat * (u_k - u_hat_k));

                % C = A - B - L
                C_val = A_val - B_const - L_val;
                if C_val > cone_viol
                    cone_viol = C_val;
                end
            end
        end
    end

    max_viol = max(wheel_viol, cone_viol);
end
