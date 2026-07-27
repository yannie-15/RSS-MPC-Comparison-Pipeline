function solution = solve_qp_with_python_hpipm( ...
    H, g, A_eq, b_eq, A_ineq, b_ineq, lb, ub, ...
    Hq, gq, uq ...
)
% SOLVE_QP_WITH_PYTHON_HPIPM
% Solve QCQP problem by calling Python hpipm_solver.py
%
% 输入：
%   H, g          : Hessian和线性项
%   A_eq, b_eq    : 等式约束 (可选)
%   A_ineq, b_ineq: 不等式约束 (可选)
%   lb, ub        : 变量上下界 (可选)
%   Hq, gq, uq    : 二次不等式约束列表 (可选)
%                   Hq: cell数组，每个元素为 (n_var, n_var) 矩阵
%                   gq: cell数组，每个元素为 (n_var, 1) 向量
%                   uq: 向量，每个元素为标量上界
%
% 输出：
%   solution: 结构体包含 x (解), status, obj_value, solve_time, solver_used

    % 获取项目根目录（RSS_proposed 的父目录就是项目根目录）
    current_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(current_dir);
    python_dir = fullfile(project_root, 'python');

    % 临时文件
    input_json = fullfile(tempdir(), 'qp_input.json');
    output_json = fullfile(tempdir(), 'qp_output.json');

    % 构造输入数据
    qp_data = struct();
    qp_data.H = H;
    qp_data.g = g;

    if nargin >= 3 && ~isempty(A_eq)
        qp_data.A = A_eq;
        qp_data.b = b_eq;
    else
        qp_data.A = [];
        qp_data.b = [];
    end

    if nargin >= 5 && ~isempty(A_ineq)
        qp_data.C = A_ineq;
        qp_data.d = b_ineq;
    else
        qp_data.C = [];
        qp_data.d = [];
    end

    if nargin >= 7 && ~isempty(lb)
        qp_data.lb = lb;
    else
        qp_data.lb = [];
    end

    if nargin >= 8 && ~isempty(ub)
        qp_data.ub = ub;
    else
        qp_data.ub = [];
    end

    % QCQP 二次约束
    if nargin >= 11 && ~isempty(Hq)
        % 将 cell 数组转为普通数组列表（JSON序列化需要）
        Hq_serialized = cell2mat_cell(Hq);
        gq_serialized = cell2mat_cell(gq);
        qp_data.Hq = Hq_serialized;
        qp_data.gq = gq_serialized;
        qp_data.uq = uq(:)';
    else
        qp_data.Hq = [];
        qp_data.gq = [];
        qp_data.uq = [];
    end

    % 写入JSON文件
    write_qp_json(input_json, qp_data);

    % 调用Python脚本（使用 venv Python 确保 hpipm_python 可用）
    python_exe = 'C:\Users\18702\sjtu-agent\sjtu-agent\.venv\Scripts\python.exe';
    if ~isfile(python_exe)
        python_exe = 'python';  % fallback to system python
    end
    python_script = fullfile(python_dir, 'hpipm_solver.py');
    cmd = sprintf('"%s" "%s" --input "%s" --output "%s"', ...
        python_exe, python_script, input_json, output_json);

    fprintf('执行命令: %s\n', cmd);
    [status, output] = system(cmd);

    if status ~= 0
        fprintf('Python执行失败:\n%s\n', output);
        error('Failed to solve QCQP with Python hpipm');
    end

    % 读取输出JSON
    if isfile(output_json)
        raw = fileread(output_json);
        solution = jsondecode(raw);
        solution.x = solution.x(:);  % 确保列向量
    else
        error('Output JSON file not found');
    end

    % 清理临时文件
    if isfile(input_json), delete(input_json); end
    if isfile(output_json), delete(output_json); end

end


%% =========================================================
% 辅助函数：将QP数据写入JSON
% ==========================================================

function write_qp_json(filename, qp_data)
    % 将cell数组中的矩阵转为嵌套列表以便JSON序列化
    data = struct();

    data.H = qp_data.H;
    data.g = qp_data.g;

    if isfield(qp_data, 'A') && ~isempty(qp_data.A)
        data.A = qp_data.A;
        data.b = qp_data.b;
    end

    if isfield(qp_data, 'C') && ~isempty(qp_data.C)
        data.C = qp_data.C;
        data.d = qp_data.d;
    end

    if isfield(qp_data, 'lb') && ~isempty(qp_data.lb)
        data.lb = qp_data.lb;
    end

    if isfield(qp_data, 'ub') && ~isempty(qp_data.ub)
        data.ub = qp_data.ub;
    end

    % QCQP 约束：Hq, gq 是 cell 数组，需要转为嵌套列表
    if isfield(qp_data, 'Hq') && ~isempty(qp_data.Hq)
        data.Hq = qp_data.Hq;
        data.gq = qp_data.gq;
        data.uq = qp_data.uq;
    end

    json_str = jsonencode(data);
    fid = fopen(filename, 'w');
    fprintf(fid, '%s', json_str);
    fclose(fid);
end


%% =========================================================
% 辅助函数：cell数组转嵌套cell（用于JSON序列化）
% ==========================================================

function result = cell2mat_cell(c)
% 将 {M1, M2, ...} 转为 {M1, M2, ...}（保持原样）
% jsonencode会自动将cell数组的每个矩阵转为嵌套列表
    result = c;
end
