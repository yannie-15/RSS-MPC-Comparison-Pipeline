% fix_paths_and_run.m
% 一次性修复脚本：清理 MATLAB path + Python sys.path 中的 RSS_V2 残留，
% 删除 __pycache__，重载模块，验证签名，运行 paper_reproduction。
%
% 用法：在 MATLAB 中执行：fix_paths_and_run

function fix_paths_and_run()
    fprintf('========== 开始一次性修复 ==========\n');

    proj_root = 'd:\PROJECT\RSS-MPC-Comparison-Pipeline-rss_hpipm';
    rss_proposed_dir = fullfile(proj_root, 'algorithms', 'RSS_proposed');

    %% 1. 清理 MATLAB path 中的 RSS_V2 / Projects\RSS 残留
    fprintf('\n[1/6] 清理 MATLAB path...\n');
    all_dirs = strsplit(path, pathsep);
    removed_matlab = 0;
    for i = 1:length(all_dirs)
        d = all_dirs{i};
        if isempty(d), continue; end
        if contains(d, 'RSS_V2') || contains(d, 'Projects\RSS')
            rmpath(d);
            fprintf('  MATLAB rmpath: %s\n', d);
            removed_matlab = removed_matlab + 1;
        end
    end
    fprintf('  MATLAB path 清理完成 (移除 %d 条)\n', removed_matlab);

    %% 2. 添加当前项目路径到 MATLAB path
    fprintf('\n[2/6] 添加当前项目路径到 MATLAB path...\n');
    addpath(fullfile(proj_root, 'paper_reproduction'));
    addpath(fullfile(proj_root, 'core'));
    addpath(fullfile(proj_root, 'batch_simulation'));
    addpath(fullfile(proj_root, 'algorithms'));
    setup_paths();

    %% 3. 清理 Python sys.path 中的 RSS_V2 残留 (用临时 .py 文件执行)
    fprintf('\n[3/6] 清理 Python sys.path...\n');
    cleanup_py = fullfile(proj_root, '_cleanup_syspath.py');
    fid = fopen(cleanup_py, 'w');
    fprintf(fid, 'import sys\n');
    fprintf(fid, 'to_remove = [p for p in list(sys.path) if "RSS_V2" in p or "Projects\\\\RSS" in p]\n');
    fprintf(fid, 'for p in to_remove:\n');
    fprintf(fid, '    sys.path.remove(p)\n');
    fprintf(fid, '    print("  Python sys.path remove:", p)\n');
    fprintf(fid, 'target = r"%s"\n', rss_proposed_dir);
    fprintf(fid, 'if target not in sys.path:\n');
    fprintf(fid, '    sys.path.append(target)\n');
    fprintf(fid, '    print("  Python sys.path append:", target)\n');
    fprintf(fid, 'print("  Python sys.path count:", len(sys.path))\n');
    fclose(fid);
    py.runpy.run_path(cleanup_py);
    delete(cleanup_py);

    %% 4. 删除 __pycache__ 中的旧 .pyc 文件
    fprintf('\n[4/6] 清理 __pycache__...\n');
    pycache_dir = fullfile(rss_proposed_dir, '__pycache__');
    if exist(pycache_dir, 'dir')
        pyc_files = dir(fullfile(pycache_dir, 'hpipm_qp_solver*.pyc'));
        for i = 1:length(pyc_files)
            fpath = fullfile(pycache_dir, pyc_files(i).name);
            delete(fpath);
            fprintf('  删除: %s\n', pyc_files(i).name);
        end
        fprintf('  __pycache__ 清理完成 (删除 %d 个 .pyc)\n', length(pyc_files));
    else
        fprintf('  __pycache__ 目录不存在，跳过\n');
    end

    %% 5. 彻底重新加载 Python 模块并验证签名 (用临时 .py 文件执行)
    fprintf('\n[5/6] 重载 Python 模块并验证签名...\n');
    reload_py = fullfile(proj_root, '_reload_solver.py');
    fid = fopen(reload_py, 'w');
    fprintf(fid, 'import sys\n');
    fprintf(fid, 'sys.modules.pop("hpipm_qp_solver", None)\n');
    fprintf(fid, 'import hpipm_qp_solver\n');
    fprintf(fid, 'print("  模块文件路径:", hpipm_qp_solver.__file__)\n');
    fprintf(fid, 'import inspect\n');
    fprintf(fid, 'sig = inspect.signature(hpipm_qp_solver.solve_ocp_qp)\n');
    fprintf(fid, 'print("  solve_ocp_qp 签名:", sig)\n');
    fprintf(fid, 'n_params = len(sig.parameters)\n');
    fprintf(fid, 'print("  参数数量:", n_params)\n');
    fprintf(fid, 'if n_params >= 19:\n');
    fprintf(fid, '    print("  签名正确")\n');
    fprintf(fid, 'else:\n');
    fprintf(fid, '    print("  签名错误: 旧版本")\n');
    fprintf(fid, '    sys.exit(1)\n');
    fclose(fid);
    py.runpy.run_path(reload_py);
    delete(reload_py);

    % 强制 MATLAB 侧重新导入 (clear python 清除 MATLAB 的 Python 模块缓存)
    try
        clear python;
    catch
    end
    % 重新触发 import, 让 py.hpipm_qp_solver 指向新模块
    % 用 py.getattr 访问 __file__ 属性
    % 用 char([95 95 102 105 108 101 95 95]) 构造 '__file__' 字符串, 避免双下划线编码问题
    file_attr = char([95 95 102 105 108 101 95 95]);  % '__file__'
    tmp_mod = py.importlib.import_module('hpipm_qp_solver');
    tmp_file = char(py.getattr(tmp_mod, file_attr));
    fprintf('  MATLAB 侧模块路径: %s\n', tmp_file);
    if contains(tmp_file, 'RSS_V2') || contains(tmp_file, 'Projects\RSS')
        fprintf('  错误: MATLAB 侧仍加载旧版本\n');
        return;
    end
    fprintf('  MATLAB 侧模块引用已更新\n');

    %% 6. 清除 MATLAB 函数缓存并运行
    fprintf('\n[6/6] 运行 paper_reproduction...\n');
    clear functions;
    cd(proj_root);

    % 验证 MATLAB 加载的文件路径
    fprintf('  MATLAB which paper_reproduction: %s\n', which('paper_reproduction'));

    paper_reproduction({'proposed-3iter'});

    fprintf('\n========== 修复脚本完成 ==========\n');
end
