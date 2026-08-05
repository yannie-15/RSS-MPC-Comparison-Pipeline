function path = generateReference(params, num_points)
    if nargin < 2
        num_points = params.num_path_pts;
    end

    ctrl_pts = params.ctrl_pts;
    t = linspace(0, 1, num_points);
    n = size(ctrl_pts, 1) - 1;
    path_x = zeros(size(t));
    path_y = zeros(size(t));

    for i = 0:n
        coeff = nchoosek(n, i) .* (t .^ i) .* ((1 - t) .^ (n - i));
        path_x = path_x + coeff .* ctrl_pts(i + 1, 1);
        path_y = path_y + coeff .* ctrl_pts(i + 1, 2);
    end

    dx = diff(path_x);
    dy = diff(path_y);
    s_list = [0, cumsum(sqrt(dx .^ 2 + dy .^ 2))];
    Lp = s_list(end);
    path_theta = (2 * pi / Lp^2) * s_list.^2;
    path = [path_x; path_y; path_theta];
end
