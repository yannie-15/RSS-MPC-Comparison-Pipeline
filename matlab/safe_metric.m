function val = safe_metric(r, field_name)
%SAFE_METRIC  安全地从 result.metrics 中读取一个数值字段
%
% 输入:
%   r          : 单个 summary 结构体 (含 .metrics 字段)
%   field_name : 字段名 (字符串)
%
% 输出:
%   val : 数值; 字段不存在 / 非数值时返回 NaN

    if isfield(r, 'metrics') && isfield(r.metrics, field_name)
        val = r.metrics.(field_name);
        if ~isnumeric(val), val = NaN; end
    else
        val = NaN;
    end
end
