% CONFIG_HPIPM
% HPIPM求解器的配置文件
% 在config.m中调用此文件以启用HPIPM支持

function hpipm_config = config_hpipm()
    % 返回HPIPM相关的配置参数
    
    hpipm_config = struct();
    
    %% =====================================================
    % 求解器选择
    % ======================================================
    
    % 是否使用HPIPM求解器
    % 可选值：
    %   false   : 使用CVX + SDPT3（默认，最兼容）
    %   'python': 使用Python HPIPM接口（推荐）
    %   'mex'   : 使用HPIPM MEX接口（需要编译）
    hpipm_config.solver_type = 'python';  % 或 false, 'python', 'mex'
    
    % 如果主求解器失败是否回退到CVX
    hpipm_config.fallback_to_cvx = true;
    
    %% =====================================================
    % 算法参数
    % ======================================================
    
    % IPM最大迭代数
    hpipm_config.ipm_max_iter = 50;
    
    % 收敛容差
    hpipm_config.ipm_tolerance = 1e-6;
    
    % 使用热启动（上次迭代的解作为初值）
    hpipm_config.warm_start = true;
    
    % 线性求解器类型
    % 可选：'choleski', 'qr'
    hpipm_config.lin_alg_type = 'choleski';
    
    %% =====================================================
    % 性能参数
    % ======================================================
    
    % 求解器超时时间（秒）
    hpipm_config.solver_timeout = 10;
    
    % 是否记录求解器详细日志
    hpipm_config.verbose = false;
    
    % 是否记录求解时间
    hpipm_config.measure_time = true;
    
    %% =====================================================
    % Python特定参数
    % ======================================================
    
    % Python解释器路径（如果为空则使用系统默认）
    hpipm_config.python_executable = '';  % 例如：'C:\Python39\python.exe'
    
    % Python脚本目录
    hpipm_config.python_script_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'python');
    
    %% =====================================================
    % MEX特定参数
    % ======================================================
    
    % MEX文件目录
    hpipm_config.mex_lib_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'third_party', 'hpipm', 'interfaces', 'matlab_octave');
    
    %% =====================================================
    % 调试参数
    % ======================================================
    
    % 保存QP问题到文件用于调试
    hpipm_config.save_qp_to_file = false;
    hpipm_config.qp_debug_dir = fullfile(tempdir(), 'qp_debug');
    
    % 比较CVX和HPIPM的结果（用于验证）
    hpipm_config.compare_solvers = false;
    
end


%% =========================================================
% 使用示例
% ==========================================================

% 在main.m中使用：
%
% % 加载基础配置
% params = config();
%
% % 加载HPIPM配置
% hpipm_cfg = config_hpipm();
%
% % 将HPIPM配置加入params
% params.hpipm = hpipm_cfg;
%
% % 使用control_RSS_v2
% [u, vel, new_state_dot] = control_RSS_v2(path, step, state_dot, state, params);


%% =========================================================
% 预设配置
% ==========================================================

function cfg = config_hpipm_fast()
    % 快速求解配置（精度降低，速度快）
    cfg = config_hpipm();
    cfg.ipm_tolerance = 1e-4;
    cfg.ipm_max_iter = 30;
    cfg.verbose = false;
end

function cfg = config_hpipm_accurate()
    % 精确求解配置（精度高，速度慢）
    cfg = config_hpipm();
    cfg.ipm_tolerance = 1e-8;
    cfg.ipm_max_iter = 100;
    cfg.verbose = true;
end

function cfg = config_hpipm_balanced()
    % 均衡配置（中等精度和速度）
    cfg = config_hpipm();
    cfg.ipm_tolerance = 1e-6;
    cfg.ipm_max_iter = 50;
    cfg.verbose = false;
end
