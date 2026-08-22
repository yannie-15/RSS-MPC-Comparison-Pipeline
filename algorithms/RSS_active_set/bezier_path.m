% 辅助子函数：生成贝塞尔路径以及期望朝向路径，path为3*N矩阵
function [path] = bezier_path(ctrl_pts, num_pts)
    t = linspace(0, 1, 100);
    n = size(ctrl_pts,1)-1;
    path_x = 0; path_y = 0;
    for i = 0:n
        coeff = nchoosek(n,i).*t.^i.*(1-t).^(n-i);
        path_x = path_x + coeff*ctrl_pts(i+1,1);
        path_y = path_y + coeff*ctrl_pts(i+1,2);
    end
     dx = diff(path_x); dy = diff(path_y);
    s_list =  [0,cumsum(sqrt(dx.^2 + dy.^2))];
    Lp = s_list( end );
    path_theta = (2*pi / Lp^2) * s_list.^2;
    path =[path_x;path_y;path_theta];    
end