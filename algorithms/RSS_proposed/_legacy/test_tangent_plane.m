function test_tangent_plane()
% TEST_TANGENT_PLANE  验证 construct_ocp_qp_from_rss.m 中转向锥约束的线性化正确性
%
% 验证内容:
%   1. 相切性: 在切点 (xh, uh) 处, 线性化约束值 = 精确约束值
%      C_lin(xh,uh) = grad_C'*[xh;uh] - rhs = C_exact(xh,uh)
%   2. 梯度一致性: 解析梯度 = 数值梯度
%      grad_v C = 2*Mn*v + Mn*u - T'*ell_hat
%      grad_u C = Mn*(v+u) - U'*ell_hat

    addpath(fileparts(mfilename('fullpath')));
    params = config();

    K = 6;
    dt = params.dt;
    phidotmax = params.phidotmax;
    vimax = params.vimax;
    wheel_pos = params.wheel_pos;
    num_wheels = size(wheel_pos, 1);
    delta_theta = dt * phidotmax;

    % 测试数据 (与 paper baseline 一致)
    state = [0.05, 0.1, 0.2]';
    v0 = [0.01; 0.01; 0.01];
    step = 1;
    path = generateReference(params, params.num_path_pts);

    % 随机 u_anchor 和 u_cut (确保不同, 测试 L_cut ≠ 0 情况)
    rng(42);
    u_anchor = 0.01 * (rand(3, K) - 0.5);
    u_cut = 0.01 * (rand(3, K) - 0.5);

    % 构造 OCP QP
    ocp = construct_ocp_qp_from_rss(path, step, v0, state, u_anchor, u_cut, params);

    % 预计算
    Hn = cell(1, num_wheels);
    for n = 1:num_wheels
        Hn{n} = [1, 0, -wheel_pos(n,2); 0, 1, wheel_pos(n,1)];
    end
    R1 = [sin(delta_theta), -cos(delta_theta); cos(delta_theta),  sin(delta_theta)];
    R2 = [sin(delta_theta),  cos(delta_theta); -cos(delta_theta), sin(delta_theta)];

    % nu_hat 序列
    nu_hat_anchor = zeros(3, K+1);
    nu_hat_anchor(:, 1) = v0;
    for k = 1:K
        nu_hat_anchor(:, k+1) = nu_hat_anchor(:, k) + u_anchor(:, k);
    end
    nu_hat_cut = zeros(3, K+1);
    nu_hat_cut(:, 1) = v0;
    for k = 1:K
        nu_hat_cut(:, k+1) = nu_hat_cut(:, k) + u_cut(:, k);
    end

    max_tangent_err = 0;
    max_grad_err = 0;
    n_checked = 0;

    % 遍历所有转向锥约束 (stage 0..K-1, 每stage 2*N 条)
    for n_stage = 0:K-1
        k = n_stage + 1;
        x_anchor = nu_hat_anchor(:, k);
        u_anchor_k = u_anchor(:, k);
        xh = nu_hat_cut(:, k);
        uh = u_cut(:, k);

        for i = 1:num_wheels
            Hi = Hn{i};
            Mn = Hi' * Hi;

            for gg = 1:2
                if gg == 1; Rg = R1; else; Rg = R2; end
                T = (eye(2) + Rg) * Hi;
                U = Rg * Hi;

                ell_hat = T * x_anchor + U * u_anchor_k;
                B_const = 0.5 * (ell_hat' * ell_hat);

                % 精确约束值 C(xh, uh)
                A_xh_uh = xh'*Mn*xh + xh'*Mn*uh + 0.5*uh'*Mn*uh;
                L_cut = ell_hat' * (T*(xh - x_anchor) + U*(uh - u_anchor_k));
                C_xh_uh = A_xh_uh - B_const - L_cut;

                % 解析梯度
                grad_x_C = 2*Mn*xh + Mn*uh - T'*ell_hat;
                grad_u_C = Mn*(xh+uh) - U'*ell_hat;

                % 从 ocp 提取线性化约束
                % 约束顺序 (construct_ocp_qp_from_rss.m):
                %   stage 0:     只有转向锥 (2*N 条), 索引 1..2N
                %   stage 1..K-2: 轮速 (N) + 转向锥 (2*N), 索引 N+1..3N
                %   stage K-1:   轮速 (N) + 终端轮速 (N) + 转向锥 (2*N), 索引 2N+1..4N
                if n_stage == 0
                    cone_start = 0;  % stage 0: 只有转向锥
                elseif n_stage == K-1
                    cone_start = 2 * num_wheels;  % stage K-1: 轮速 + 终端轮速后
                else
                    cone_start = num_wheels;  % stage 1..K-2: 轮速后
                end

                cone_idx = cone_start + (i-1)*2 + gg;
                C_ni = ocp.Cmat{n_stage+1, cone_idx};
                D_ni = ocp.Dmat{n_stage+1, cone_idx};
                ug_ni = ocp.ug{n_stage+1, cone_idx};

                % 检查非空
                if isempty(C_ni)
                    warning('stage=%d, cone_idx=%d 为空, 跳过', n_stage, cone_idx);
                    continue;
                end

                % 提取 grad_x_C 和 grad_u_C (C_ni 后 3 位)
                grad_x_C_code = C_ni(4:6)';
                grad_u_C_code = D_ni';

                % ---- 验证 1: 相切性 ----
                % C_lin(xh,uh) = grad_x_C'*xh + grad_u_C'*uh - rhs
                % 应等于 C_xh_uh
                C_lin_at_tangent = grad_x_C_code'*xh + grad_u_C_code'*uh - ug_ni;
                tangent_err = abs(C_lin_at_tangent - C_xh_uh);
                max_tangent_err = max(max_tangent_err, tangent_err);

                % ---- 验证 2: 梯度一致性 (代码 vs 解析) ----
                grad_err_x = max(abs(grad_x_C_code(:) - grad_x_C(:)));
                grad_err_u = max(abs(grad_u_C_code(:) - grad_u_C(:)));
                max_grad_err = max(max_grad_err, max(grad_err_x, grad_err_u));

                % ---- 验证 3: 数值梯度 ----
                eps_num = 1e-7;
                num_grad_x = zeros(3,1);
                num_grad_u = zeros(3,1);
                for d = 1:3
                    xh_p = xh; xh_p(d) = xh_p(d) + eps_num;
                    xh_m = xh; xh_m(d) = xh_m(d) - eps_num;
                    A_p = xh_p'*Mn*xh_p + xh_p'*Mn*uh + 0.5*uh'*Mn*uh;
                    A_m = xh_m'*Mn*xh_m + xh_m'*Mn*uh + 0.5*uh'*Mn*uh;
                    L_p = ell_hat' * (T*(xh_p - x_anchor) + U*(uh - u_anchor_k));
                    L_m = ell_hat' * (T*(xh_m - x_anchor) + U*(uh - u_anchor_k));
                    C_p = A_p - B_const - L_p;
                    C_m = A_m - B_const - L_m;
                    num_grad_x(d) = (C_p - C_m) / (2*eps_num);

                    uh_p = uh; uh_p(d) = uh_p(d) + eps_num;
                    uh_m = uh; uh_m(d) = uh_m(d) - eps_num;
                    A_p = xh'*Mn*xh + xh'*Mn*uh_p + 0.5*uh_p'*Mn*uh_p;
                    A_m = xh'*Mn*xh + xh'*Mn*uh_m + 0.5*uh_m'*Mn*uh_m;
                    L_p = ell_hat' * (T*(xh - x_anchor) + U*(uh_p - u_anchor_k));
                    L_m = ell_hat' * (T*(xh - x_anchor) + U*(uh_m - u_anchor_k));
                    C_p = A_p - B_const - L_p;
                    C_m = A_m - B_const - L_m;
                    num_grad_u(d) = (C_p - C_m) / (2*eps_num);
                end
                num_grad_err_x = max(abs(num_grad_x(:) - grad_x_C(:)));
                num_grad_err_u = max(abs(num_grad_u(:) - grad_u_C(:)));
                max_grad_err = max(max_grad_err, max(num_grad_err_x, num_grad_err_u));

                n_checked = n_checked + 1;
            end
        end
    end

    max_tangent_err = max(max_tangent_err(:));
    max_grad_err = max(max_grad_err(:));

    fprintf('============== 切平面相切测试 ==============\n');
    fprintf('检查约束数: %d\n', n_checked);
    fprintf('最大相切误差: %.6e (应 < 1e-10)\n', max_tangent_err);
    fprintf('最大梯度误差: %.6e (应 < 1e-6)\n', max_grad_err);
    if max_tangent_err < 1e-10 && max_grad_err < 1e-6
        fprintf('结论: PASS (线性化正确)\n');
    else
        fprintf('结论: FAIL (线性化有误)\n');
    end
    fprintf('===========================================\n');
end
