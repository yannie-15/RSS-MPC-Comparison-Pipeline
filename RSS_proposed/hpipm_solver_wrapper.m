function [solution, info] = hpipm_solver_wrapper(qp_problem, use_hpipm)
% HPIPM_SOLVER_WRAPPER
% MATLAB wrapper for HPIPM QP/QCQP solver with CVX fallback
%
% 输入：
%   qp_problem  : 结构体，包含：
%     - H         : (n×n) Hessian矩阵
%     - g         : (n×1) 一次项系数
%     - A         : (m×n) 等式约束系数矩阵
%     - b         : (m×1) 等式约束右手边
%     - C         : (p×n) 不等式约束系数矩阵 (C*x <= d)
%     - d         : (p×1) 不等式约束右手边
%     - lb        : (n×1) 下界（可选）
%     - ub        : (n×1) 上界（可选）
%   use_hpipm   : 是否使用hpipm求解器 (默认false - 使用CVX)
%
% 输出：
%   solution    : 最优解
%   info        : 求解信息结构体
%     - solver_used   : 使用的求解器名称
%     - status        : 求解状态
%     - obj_value     : 目标函数值
%     - solve_time    : 求解时间

    if nargin < 2
        use_hpipm = false;
    end

    % 从qp_problem提取参数
    H = qp_problem.H;
    g = qp_problem.g;
    A = qp_problem.A;
    b = qp_problem.b;
    C = qp_problem.C;
    d = qp_problem.d;
    
    lb = qp_problem.lb;
    ub = qp_problem.ub;
    
    n = size(H, 1);
    
    % 计时
    solve_tic = tic;
    
    % 尝试使用hpipm（如果可用）
    if use_hpipm
        try
            [solution, status, obj_value] = solve_with_hpipm(H, g, A, b, C, d, lb, ub);
            info.solver_used = 'HPIPM';
            info.status = status;
            info.obj_value = obj_value;
            info.solve_time = toc(solve_tic);
            return;
        catch ME
            fprintf('Warning: HPIPM求解失败，切换到CVX: %s\n', ME.message);
        end
    end
    
    % 使用CVX作为后备求解器
    [solution, status, obj_value] = solve_with_cvx(H, g, A, b, C, d, lb, ub, n);
    
    info.solver_used = 'CVX';
    info.status = status;
    info.obj_value = obj_value;
    info.solve_time = toc(solve_tic);
    
end


%% =========================================================
% 使用 HPIPM 求解
% ==========================================================

function [solution, status, obj_value] = solve_with_hpipm(H, g, A, b, C, d, lb, ub)
% 调用hpipm的MEX接口求解QP问题
% 需要hpipm的MATLAB MEX文件已编译

    % 检查hpipm_qp是否可用
    if ~exist('hpipm_qp', 'file')
        error('HPIPM MEX interface not found. Please compile hpipm first.');
    end

    n = size(H, 1);
    
    % 构造HPIPM的QP结构
    % HPIPM支持形式：
    %   min 0.5*x'*H*x + g'*x
    %   s.t. A*x = b
    %        C*x <= d
    %        lb <= x <= ub
    
    qp.H = H;
    qp.g = g;
    qp.A = A;
    qp.b = b;
    qp.C = C;
    qp.d = d;
    
    if ~isempty(lb)
        qp.lb = lb;
    else
        qp.lb = -inf * ones(n, 1);
    end
    
    if ~isempty(ub)
        qp.ub = ub;
    else
        qp.ub = inf * ones(n, 1);
    end
    
    % 调用hpipm MEX求解
    [solution, ipm_info] = hpipm_qp(qp);
    
    status = ipm_info.status;
    obj_value = 0.5 * solution' * H * solution + g' * solution;
    
end


%% =========================================================
% 使用 CVX 求解（后备方案）
% ==========================================================

function [solution, status, obj_value] = solve_with_cvx(H, g, A, b, C, d, lb, ub, n)

    cvx_begin

        cvx_solver SDPT3
        
        variable x(n)
        
        minimize(0.5 * quad_form(x, H) + g' * x)
        
        subject to
        
            % 等式约束
            if ~isempty(A) && ~isempty(b)
                A * x == b;
            end
            
            % 不等式约束
            if ~isempty(C) && ~isempty(d)
                C * x <= d;
            end
            
            % 边界约束
            if ~isempty(lb)
                x >= lb;
            end
            if ~isempty(ub)
                x <= ub;
            end
    
    cvx_end
    
    solution = x;
    status = cvx_status;
    obj_value = cvx_optval;
    
end
