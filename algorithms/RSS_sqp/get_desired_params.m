    function [Pd, psit, Cc, thetad, dthetad_ds, d2thetad_ds2] = get_desired_params(s, path_x, path_y, s_list, Lp)
    s = min(s, Lp);
  
    % 插值期望位置
    x = interp1(s_list, path_x, s, 'linear');
    y = interp1(s_list, path_y, s, 'linear');
    Pd = [x, y];
   
    % 先算完整路径的导数向量，再插值到当前s点
    dx_ds_vec = gradient(path_x, s_list);
    d2x_ds2_vec = gradient(dx_ds_vec, s_list);
    dy_ds_vec = gradient(path_y, s_list);
    d2y_ds2_vec = gradient(dy_ds_vec, s_list);
    
    dx_ds = interp1(s_list, dx_ds_vec, s, 'linear');
    dy_ds = interp1(s_list, dy_ds_vec, s, 'linear');
    d2x_ds2 = interp1(s_list, d2x_ds2_vec, s, 'linear');
    d2y_ds2 = interp1(s_list, d2y_ds2_vec, s, 'linear');
    
    % 切向角和曲率
    psit = atan2(dy_ds, dx_ds);
    denom = (dx_ds^2 + dy_ds^2)^(3/2) + 1e-6;
    Cc = (dx_ds * d2y_ds2 - dy_ds * d2x_ds2) / denom;
    
    % 期望朝向
    thetad = (2*pi / Lp^2) * s^2;
    dthetad_ds = (4*pi / Lp^2) * s;
    d2thetad_ds2 = 4*pi / Lp^2;
end