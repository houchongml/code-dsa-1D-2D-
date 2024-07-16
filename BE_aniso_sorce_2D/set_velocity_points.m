function [vel_points, vel_weights] = set_velocity_points(Nv)
    % 使用均匀分布在[0, 2*pi]上的角度
    theta = linspace(0, 2 * pi, Nv + 1);
    theta(end) = []; % 去掉最后一个点，避免重复
    
    % 速度点
    vel_points = [cos(theta); sin(theta)];
    
    % 权重
    vel_weights = ones(1, Nv) * (2 * pi / Nv);
end