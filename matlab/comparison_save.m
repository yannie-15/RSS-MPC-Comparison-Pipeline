function comparison_save(comparison, checkpoint_file)
%COMPARISON_SAVE  增量保存 checkpoint (覆盖写)
%
% 保持与原 compare_algorithms.m 的 save_checkpoint 完全等价的行为,
% 仅做名字上的拆分以便 main.m / run_batch_simulation.m 共用。
%
% 输入:
%   comparison      : comparison 结构体
%   checkpoint_file : checkpoint 文件全路径

    save(checkpoint_file, 'comparison', '-v7.3');
end
