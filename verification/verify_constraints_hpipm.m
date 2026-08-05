function report = verify_constraints_hpipm()
% VERIFY_CONSTRAINTS_HPIPM
% 验证 RSS+HPIPM 控制器在每个 MPC 时刻最终输出的预测解是否满足原始约束。
%
% 检查内容：
%   1. 轮速约束
%        ||H_n * nu_k||_2 <= vimax
%
%   2. RSS 原始双线性转向锥约束
%        (R_i * z_{n,k-1})' * z_{n,k} >= 0,  i = 1,2
%        z_{n,k} = H_n * nu_k
%        R_1 = R(pi/2-delta_theta), R_2 = R_1'
%
%   3. 方向角变化辅助检查
%        |wrapToPi(theta_k-theta_{k-1})| <= delta_theta
%      仅当相邻轮速均大于低速阈值时检查；低速情况单独记为 undefined。
%
%   4. 控制器接口一致性
%        bodyVelocity == current_nu + u(:,1)
%        worldVelocity == T_world_body * bodyVelocity
%
% 说明：
%   - 本脚本验证每个 MPC 时刻 control_RSS 最终返回的 u_full。
%   - 若需验证每次 RSS 内部迭代，control_RSS 还需在 diag 中返回每次迭代的 u。
%   - 未执行或未检查的条目使用 NaN，不会被误统计为“零违反”。
%   - 求解失败后不使用失败解推进状态。
%
% 用法：
%   cd('D:\PROJECT\RSS_V2\matlab')
%   setup_paths
%   report = verify_constraints_hpipm();

%% 路径设置
this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(this_dir);

% 添加 core/ (共享工具) 和 algorithms/RSS_proposed/ (控制器 + 求解器)
core_dir = fullfile(project_root, 'core');
if exist(core_dir, 'dir')
    addpath(core_dir);
end
rss_proposed_dir = fullfile(project_root, 'algorithms', 'RSS_proposed');
if exist(rss_proposed_dir, 'dir')
    addpath(rss_proposed_dir);
end

% Force reload Python HPIPM solver module
% (MATLAB Engine caches Python modules; must reload after .py changes;
%  hpipm_qp_solver.py add_dll_directory code must re-execute to set DLL path)
try
    % hpipm_qp_solver.py 现与 control_RSS.m 同目录 (algorithms/RSS_proposed/)
    sys_mod = py.importlib.import_module('sys');
    py.getattr(sys_mod, 'path').append(rss_proposed_dir);
    solver_mod = py.importlib.import_module('hpipm_qp_solver');
    py.importlib.reload(solver_mod);
    hpipm_ok = py.getattr(solver_mod, '_HPIPM_OK');
    fprintf('[verify] Python HPIPM module reloaded, _HPIPM_OK=%d\n', ...
        double(hpipm_ok));
    if ~double(hpipm_ok)
        hpipm_err = py.getattr(solver_mod, '_HPIPM_ERR');
        error('HPIPM Python module load failed: %s', char(hpipm_err));
    end
catch err
    fprintf('[verify] Warning: Python module reload failed: %s\n', err.message);
end

%% 初始化
params = config();
path = generateReference(params, params.num_path_pts);

assert(isnumeric(path) && size(path,1) >= 3, ...
    'generateReference 的输出应至少为 3xN: [x;y;psi]。');
assert(all(isfinite(path(:))), ...
    '参考轨迹包含 NaN 或 Inf。');

state = [0.05; 0.1; 0.2];
last_vel = [0.01; 0.01; 0.01];

num_steps = params.num_steps;
num_wheels = size(params.wheel_pos,1);

assert(num_wheels >= 1, ...
    'params.wheel_pos 至少应包含一个轮位。');
assert(size(params.wheel_pos,2) == 2, ...
    'params.wheel_pos 应为 N x 2。');

vimax = params.vimax;
delta_theta = params.dt * params.phidotmax;

assert(params.dt > 0, 'params.dt 必须为正。');
assert(params.phidotmax > 0, ...
    'params.phidotmax 必须为正，单位应为 rad/s。');
assert(delta_theta > 0 && delta_theta < pi/2, ...
    '必须满足 0 < dt*phidotmax < pi/2。');

% 轮子特征矩阵
H_cell = cell(1,num_wheels);
for n = 1:num_wheels
    dx = params.wheel_pos(n,1);
    dy = params.wheel_pos(n,2);
    H_cell{n} = [1,0,-dy; 0,1,dx];
end

% RSS 原始锥约束的两个旋转矩阵
R1 = [sin(delta_theta), -cos(delta_theta); ...
      cos(delta_theta),  sin(delta_theta)];
R2 = R1';
R_set = {R1,R2};

%% 数值容差
tol.speed_abs = 1e-7;
tol.speed_rel = 1e-6;
tol.cone_abs = 1e-9;
tol.cone_rel = 1e-7;
tol.angle_abs = 1e-7;
tol.angle_rel = 1e-6;
tol.direction_speed_eps = 1e-6;
tol.interface = 1e-8;

%% 结果存储
viol.attempted = false(num_steps,1);
viol.checked = false(num_steps,1);
viol.solve_failed = false(num_steps,1);
viol.hpipm_status = repmat({'NotRun'},num_steps,1);
viol.exception = repmat({''},num_steps,1);

viol.horizon_K = nan(num_steps,1);

viol.speed_max = nan(num_steps,1);
viol.speed_max_violation = nan(num_steps,1);
viol.speed_fail_count = nan(num_steps,1);
viol.speed_worst_k = nan(num_steps,1);
viol.speed_worst_n = nan(num_steps,1);

viol.cone_min_margin = nan(num_steps,1);
viol.cone_max_violation = nan(num_steps,1);
viol.cone_fail_count = nan(num_steps,1);
viol.cone_worst_k = nan(num_steps,1);
viol.cone_worst_n = nan(num_steps,1);
viol.cone_worst_i = nan(num_steps,1);

viol.angle_max = nan(num_steps,1);
viol.angle_max_violation = nan(num_steps,1);
viol.angle_fail_count = nan(num_steps,1);
viol.angle_undefined_count = nan(num_steps,1);
viol.angle_worst_k = nan(num_steps,1);
viol.angle_worst_n = nan(num_steps,1);

viol.body_velocity_consistency_error = nan(num_steps,1);
viol.world_velocity_consistency_error = nan(num_steps,1);

fprintf('============================================================\n');
fprintf('              HPIPM RSS 原始约束验证开始\n');
fprintf('============================================================\n');
fprintf('vimax       = %.9f m/s\n',vimax);
fprintf('phidotmax   = %.9f rad/s\n',params.phidotmax);
fprintf('dt          = %.9f s\n',params.dt);
fprintf('delta_theta = %.9f rad (%.6f deg)\n', ...
    delta_theta,rad2deg(delta_theta));
fprintf('steps       = %d\n',num_steps);
fprintf('wheels      = %d\n\n',num_wheels);

%% 闭环验证
for step = 1:num_steps
    viol.attempted(step) = true;
    current_nu = last_vel;

    try
        [u_full,worldVelocity,bodyVelocity,diag] = ...
            control_RSS(path,step,last_vel,state);

        [solver_ok,status_text] = extract_solver_status(diag);
        viol.hpipm_status{step} = status_text;
        viol.solve_failed(step) = ~solver_ok;

        if ~solver_ok
            fprintf('[Step %d] HPIPM 求解失败: %s\n', ...
                step,status_text);
            break;
        end

        %% 输出数据检查
        assert(isnumeric(u_full) && size(u_full,1) == 3, ...
            'u_full 应为 3xK 数值矩阵。');
        assert(size(u_full,2) >= 1, ...
            'u_full 的预测时域 K 必须大于等于 1。');
        assert(all(isfinite(u_full(:))), ...
            'u_full 包含 NaN 或 Inf。');

        bodyVelocity = bodyVelocity(:);
        worldVelocity = worldVelocity(:);

        assert(numel(bodyVelocity) == 3, ...
            'bodyVelocity 应包含 3 个元素。');
        assert(numel(worldVelocity) == 3, ...
            'worldVelocity 应包含 3 个元素。');
        assert(all(isfinite(bodyVelocity)), ...
            'bodyVelocity 包含 NaN 或 Inf。');
        assert(all(isfinite(worldVelocity)), ...
            'worldVelocity 包含 NaN 或 Inf。');

        K = size(u_full,2);
        viol.horizon_K(step) = K;

        %% 重构 nu 序列
        nu_seq = zeros(3,K);
        nu_seq(:,1) = current_nu + u_full(:,1);

        for k = 1:K-1
            nu_seq(:,k+1) = ...
                nu_seq(:,k) + u_full(:,k+1);
        end

        %% 控制器接口一致性
        body_error = norm( ...
            bodyVelocity-nu_seq(:,1),inf);

        psi = state(3);
        T_world_body = [ ...
            cos(psi),-sin(psi),0; ...
            sin(psi), cos(psi),0; ...
            0,        0,       1];

        expected_world_velocity = ...
            T_world_body * bodyVelocity;

        world_error = norm( ...
            worldVelocity-expected_world_velocity,inf);

        viol.body_velocity_consistency_error(step) = body_error;
        viol.world_velocity_consistency_error(step) = world_error;

        if body_error > tol.interface
            error('HPIPMVerification:BodyVelocityMismatch', ...
                ['bodyVelocity 与 current_nu+u(:,1) 不一致，', ...
                 '误差为 %.6e。'],body_error);
        end

        if world_error > tol.interface
            error('HPIPMVerification:WorldVelocityMismatch', ...
                ['worldVelocity 与坐标变换结果不一致，', ...
                 '误差为 %.6e。'],world_error);
        end

        %% 原始约束检查
        check = check_original_constraints( ...
            current_nu,nu_seq,H_cell,R_set, ...
            vimax,delta_theta,tol);

        viol.speed_max(step) = check.speed_max;
        viol.speed_max_violation(step) = ...
            check.speed_max_violation;
        viol.speed_fail_count(step) = ...
            check.speed_fail_count;
        viol.speed_worst_k(step) = check.speed_worst_k;
        viol.speed_worst_n(step) = check.speed_worst_n;

        viol.cone_min_margin(step) = check.cone_min_margin;
        viol.cone_max_violation(step) = ...
            check.cone_max_violation;
        viol.cone_fail_count(step) = ...
            check.cone_fail_count;
        viol.cone_worst_k(step) = check.cone_worst_k;
        viol.cone_worst_n(step) = check.cone_worst_n;
        viol.cone_worst_i(step) = check.cone_worst_i;

        viol.angle_max(step) = check.angle_max;
        viol.angle_max_violation(step) = ...
            check.angle_max_violation;
        viol.angle_fail_count(step) = ...
            check.angle_fail_count;
        viol.angle_undefined_count(step) = ...
            check.angle_undefined_count;
        viol.angle_worst_k(step) = check.angle_worst_k;
        viol.angle_worst_n(step) = check.angle_worst_n;

        viol.checked(step) = true;

        if check.speed_fail_count > 0
            fprintf([ ...
                '[Step %d] 轮速约束违反: count=%d, ', ...
                'max_speed=%.9f, max_violation=%.6e, ', ...
                'worst=(k=%d,n=%d)\n'], ...
                step,check.speed_fail_count, ...
                check.speed_max,check.speed_max_violation, ...
                check.speed_worst_k,check.speed_worst_n);
        end

        if check.cone_fail_count > 0
            fprintf([ ...
                '[Step %d] 原始锥约束违反: count=%d, ', ...
                'min_margin=%.6e, max_violation=%.6e, ', ...
                'worst=(k=%d,n=%d,i=%d)\n'], ...
                step,check.cone_fail_count, ...
                check.cone_min_margin, ...
                check.cone_max_violation, ...
                check.cone_worst_k,check.cone_worst_n, ...
                check.cone_worst_i);
        end

        if check.angle_fail_count > 0
            fprintf([ ...
                '[Step %d] 方向角变化违反: count=%d, ', ...
                'max_angle=%.9f, max_violation=%.6e, ', ...
                'worst=(k=%d,n=%d)\n'], ...
                step,check.angle_fail_count, ...
                check.angle_max,check.angle_max_violation, ...
                check.angle_worst_k,check.angle_worst_n);
        end

        %% 使用已检查的成功解推进状态
        state = propagateState(state,worldVelocity,params);
        last_vel = bodyVelocity;

    catch ME
        viol.solve_failed(step) = true;
        viol.hpipm_status{step} = ...
            ['Exception: ',ME.identifier];
        viol.exception{step} = getReport( ...
            ME,'extended','hyperlinks','off');

        fprintf('[Step %d] 异常: %s\n',step,ME.message);
        break;
    end

    if mod(step,10) == 0 || step == 1
        fprintf('Step %d/%d 完成\n',step,num_steps);
    end
end

%% 汇总
attempted_count = nnz(viol.attempted);
checked_count = nnz(viol.checked);
failed_count = nnz(viol.solve_failed & viol.attempted);
valid_mask = viol.checked & ~viol.solve_failed;

speed_fail_step_count = nnz( ...
    valid_mask & viol.speed_fail_count > 0);
cone_fail_step_count = nnz( ...
    valid_mask & viol.cone_fail_count > 0);
angle_fail_step_count = nnz( ...
    valid_mask & viol.angle_fail_count > 0);

undefined_angle_transition_count = safe_sum( ...
    viol.angle_undefined_count(valid_mask));

global_max_speed = safe_max(viol.speed_max(valid_mask));
global_max_speed_violation = safe_max( ...
    viol.speed_max_violation(valid_mask));

global_min_cone_margin = safe_min( ...
    viol.cone_min_margin(valid_mask));
global_max_cone_violation = safe_max( ...
    viol.cone_max_violation(valid_mask));

global_max_angle = safe_max(viol.angle_max(valid_mask));
global_max_angle_violation = safe_max( ...
    viol.angle_max_violation(valid_mask));

all_steps_checked = checked_count == num_steps;
no_solver_failure = failed_count == 0;
no_speed_failure = speed_fail_step_count == 0;
no_cone_failure = cone_fail_step_count == 0;
no_angle_failure = angle_fail_step_count == 0;

verification_passed = ...
    all_steps_checked && ...
    no_solver_failure && ...
    no_speed_failure && ...
    no_cone_failure && ...
    no_angle_failure;

fprintf('\n');
fprintf('============================================================\n');
fprintf('              HPIPM RSS 原始约束验证汇总\n');
fprintf('============================================================\n');
fprintf('计划仿真步数       : %d\n',num_steps);
fprintf('尝试步数           : %d\n',attempted_count);
fprintf('完成约束检查步数   : %d\n',checked_count);
fprintf('求解/运行失败步数  : %d\n\n',failed_count);

fprintf('轮速约束\n');
fprintf('  全局最大轮速     : %.9f m/s\n',global_max_speed);
fprintf('  最大超出量       : %.9e m/s\n', ...
    global_max_speed_violation);
fprintf('  违反步数         : %d/%d\n\n', ...
    speed_fail_step_count,checked_count);

fprintf('RSS 原始双线性锥约束\n');
fprintf('  全局最小 margin  : %.9e\n',global_min_cone_margin);
fprintf('  最大负向违反量   : %.9e\n', ...
    global_max_cone_violation);
fprintf('  违反步数         : %d/%d\n\n', ...
    cone_fail_step_count,checked_count);

fprintf('方向角变化辅助检查\n');
fprintf('  角度上限         : %.9f rad (%.6f deg)\n', ...
    delta_theta,rad2deg(delta_theta));
fprintf('  最大可定义角变化 : %.9f rad\n',global_max_angle);
fprintf('  最大超出量       : %.9e rad\n', ...
    global_max_angle_violation);
fprintf('  违反步数         : %d/%d\n', ...
    angle_fail_step_count,checked_count);
fprintf('  低速未定义转移数 : %d\n\n', ...
    undefined_angle_transition_count);

if failed_count > 0
    fprintf('失败详情\n');
    failed_steps = find(viol.solve_failed & viol.attempted);
    for idx = 1:numel(failed_steps)
        s = failed_steps(idx);
        fprintf('  Step %d: %s\n', ...
            s,viol.hpipm_status{s});
        if ~isempty(viol.exception{s})
            fprintf('%s\n',viol.exception{s});
        end
    end
    fprintf('\n');
end

fprintf('============================================================\n');
fprintf('                         结论\n');
fprintf('============================================================\n');

if verification_passed
    fprintf('PASS: HPIPM 最终输出通过原始约束验证。\n');
    fprintf('  - 所有 %d 个 MPC 步均成功并完成检查。\n', ...
        num_steps);
    fprintf('  - 轮速约束无显著违反。\n');
    fprintf('  - RSS 原始双线性锥约束无显著违反。\n');
    fprintf('  - 可定义的方向角变化无显著违反。\n');

    if undefined_angle_transition_count > 0
        fprintf([ ...
            '  - 注意：%d 个低速转移无法稳定定义方向角，', ...
            '其可行性以原始锥约束结果为准。\n'], ...
            undefined_angle_transition_count);
    end
else
    fprintf('FAIL: HPIPM 原始约束验证未通过。\n');

    if ~all_steps_checked
        fprintf('  - 只完成 %d/%d 步约束检查。\n', ...
            checked_count,num_steps);
    end
    if failed_count > 0
        fprintf('  - 求解或运行失败 %d 步。\n',failed_count);
    end
    if speed_fail_step_count > 0
        fprintf('  - 轮速约束违反 %d 步。\n', ...
            speed_fail_step_count);
    end
    if cone_fail_step_count > 0
        fprintf('  - 原始双线性锥约束违反 %d 步。\n', ...
            cone_fail_step_count);
    end
    if angle_fail_step_count > 0
        fprintf('  - 方向角变化辅助检查违反 %d 步。\n', ...
            angle_fail_step_count);
    end
end

%% 保存报告
report = viol;
report.params = params;
report.tolerance = tol;
report.num_steps = num_steps;
report.attempted_count = attempted_count;
report.completed_steps = checked_count;
report.failed_count = failed_count;
report.vimax = vimax;
report.delta_theta = delta_theta;
report.speed_fail_step_count = speed_fail_step_count;
report.cone_fail_step_count = cone_fail_step_count;
report.angle_fail_step_count = angle_fail_step_count;
report.undefined_angle_transition_count = ...
    undefined_angle_transition_count;
report.verification_passed = verification_passed;

results_dir = fullfile( ...
    project_root,'verification','results');

if ~exist(results_dir,'dir')
    mkdir(results_dir);
end

result_file = fullfile( ...
    results_dir,'hpipm_constraint_verification.mat');

save(result_file,'report');
fprintf('\n结果已保存到：%s\n',result_file);
end

%% 检查原始约束
function check = check_original_constraints( ...
    current_nu,nu_seq,H_cell,R_set, ...
    vimax,delta_theta,tol)

num_wheels = numel(H_cell);
K = size(nu_seq,2);

%% 轮速约束
speed_max = -inf;
speed_max_violation = 0;
speed_fail_count = 0;
speed_worst_k = NaN;
speed_worst_n = NaN;

speed_tol = tol.speed_abs + ...
    tol.speed_rel * max(1,vimax);

for k = 1:K
    for n = 1:num_wheels
        z_curr = H_cell{n} * nu_seq(:,k);
        speed_curr = norm(z_curr,2);

        if speed_curr > speed_max
            speed_max = speed_curr;
            speed_worst_k = k;
            speed_worst_n = n;
        end

        violation = max(0,speed_curr-vimax);
        speed_max_violation = max( ...
            speed_max_violation,violation);

        if speed_curr > vimax + speed_tol
            speed_fail_count = speed_fail_count + 1;
        end
    end
end

%% RSS 原始双线性锥约束
cone_min_margin = inf;
cone_max_violation = 0;
cone_fail_count = 0;
cone_worst_k = NaN;
cone_worst_n = NaN;
cone_worst_i = NaN;

for n = 1:num_wheels
    z_prev = H_cell{n} * current_nu;

    for k = 1:K
        z_curr = H_cell{n} * nu_seq(:,k);

        for i = 1:numel(R_set)
            margin = (R_set{i} * z_prev)' * z_curr;

            if margin < cone_min_margin
                cone_min_margin = margin;
                cone_worst_k = k;
                cone_worst_n = n;
                cone_worst_i = i;
            end

            violation = max(0,-margin);
            cone_max_violation = max( ...
                cone_max_violation,violation);

            scale = max(1,norm(z_prev)*norm(z_curr));
            cone_tol = tol.cone_abs + tol.cone_rel*scale;

            if margin < -cone_tol
                cone_fail_count = cone_fail_count + 1;
            end
        end

        z_prev = z_curr;
    end
end

%% 方向角变化辅助检查
angle_max = 0;
angle_max_violation = 0;
angle_fail_count = 0;
angle_undefined_count = 0;
angle_worst_k = NaN;
angle_worst_n = NaN;

angle_tol = tol.angle_abs + ...
    tol.angle_rel * max(1,delta_theta);

for n = 1:num_wheels
    z_prev = H_cell{n} * current_nu;

    for k = 1:K
        z_curr = H_cell{n} * nu_seq(:,k);

        speed_prev = norm(z_prev,2);
        speed_curr = norm(z_curr,2);

        if speed_prev <= tol.direction_speed_eps || ...
           speed_curr <= tol.direction_speed_eps
            angle_undefined_count = ...
                angle_undefined_count + 1;
        else
            angle_prev = atan2(z_prev(2),z_prev(1));
            angle_curr = atan2(z_curr(2),z_curr(1));

            d_angle = wrap_to_pi_local( ...
                angle_curr-angle_prev);
            abs_d_angle = abs(d_angle);

            if abs_d_angle > angle_max
                angle_max = abs_d_angle;
                angle_worst_k = k;
                angle_worst_n = n;
            end

            violation = max( ...
                0,abs_d_angle-delta_theta);
            angle_max_violation = max( ...
                angle_max_violation,violation);

            if abs_d_angle > delta_theta + angle_tol
                angle_fail_count = angle_fail_count + 1;
            end
        end

        z_prev = z_curr;
    end
end

check.speed_max = speed_max;
check.speed_max_violation = speed_max_violation;
check.speed_fail_count = speed_fail_count;
check.speed_worst_k = speed_worst_k;
check.speed_worst_n = speed_worst_n;

check.cone_min_margin = cone_min_margin;
check.cone_max_violation = cone_max_violation;
check.cone_fail_count = cone_fail_count;
check.cone_worst_k = cone_worst_k;
check.cone_worst_n = cone_worst_n;
check.cone_worst_i = cone_worst_i;

check.angle_max = angle_max;
check.angle_max_violation = angle_max_violation;
check.angle_fail_count = angle_fail_count;
check.angle_undefined_count = angle_undefined_count;
check.angle_worst_k = angle_worst_k;
check.angle_worst_n = angle_worst_n;
end

%% 提取求解状态 (适配当前 control_RSS.m 的 diagnostics 结构)
function [solver_ok,status_text] = extract_solver_status(diag)
solver_ok = false;
status_text = 'Unknown';

if ~isstruct(diag)
    return;
end

% 当前 control_RSS.m 的 diagnostics 结构:
%   diag.iterations.status = cell(1, max_iter)  % 字符串 cell 数组
%   'Solved' / 'Inaccurate/Solved' = 成功, 其他 = 失败
if isfield(diag,'iterations') && isstruct(diag.iterations) ...
        && isfield(diag.iterations,'status') ...
        && ~isempty(diag.iterations.status)
    status_cells = diag.iterations.status;
    if iscell(status_cells)
        last_status = status_cells{end};
    else
        last_status = status_cells(end);
    end
    status_text = value_to_text(last_status);

    % 'Solved' 或 'Inaccurate/Solved' 视为成功
    if ischar(last_status) || isstring(last_status)
        s = lower(char(last_status));
        if strcmp(s,'solved') || ...
           strcmp(s,'inaccurate/solved') || ...
           contains(s,'solved')
            solver_ok = true;
        end
    end
end

% 兼容: 若有 success 字段直接用
if isfield(diag,'success') && ~isempty(diag.success)
    solver_ok = logical(diag.success(1));
end

% 兼容: 若有 finalStatus 字段
if isfield(diag,'finalStatus') && ~isempty(diag.finalStatus)
    status_text = value_to_text(diag.finalStatus);
end
end

%% 值转文本
function s = value_to_text(v)
if ischar(v)
    s = v;
elseif isstring(v)
    s = char(v);
elseif isnumeric(v)
    if isscalar(v)
        s = num2str(v);
    else
        s = mat2str(v);
    end
else
    s = 'Unknown';
end
end

%% 角度归一化到 (-pi, pi]
function a = wrap_to_pi_local(a)
a = mod(a + pi, 2*pi) - pi;
end

%% 安全求和 (忽略 NaN)
function s = safe_sum(x)
if isempty(x)
    s = 0;
else
    s = sum(x,'omitnan');
    if isnan(s), s = 0; end
end
end

%% 安全最大值 (空则 NaN)
function m = safe_max(x)
if isempty(x)
    m = NaN;
else
    m = max(x,[],'omitnan');
end
end

%% 安全最小值 (空则 NaN)
function m = safe_min(x)
if isempty(x)
    m = NaN;
else
    m = min(x,[],'omitnan');
end
end
