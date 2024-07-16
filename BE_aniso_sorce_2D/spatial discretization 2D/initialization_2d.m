function f = initialization_2d(fcn_ic, ...
                               reference_mass, ...
                               quad_point_x, quad_weight_x, ...
                               quad_point_y, quad_weight_y, ...
                               v, ...
                               x_center, y_center, dx, dy)
    n_basis = 3; % 基函数数量，对应于 1, x, y
    local_rhs = zeros(n_basis, 1);
    x_quad = x_center + 0.5 * quad_point_x * dx; % 实际网格点 (x 方向)
    y_quad = y_center + 0.5 * quad_point_y * dy; % 实际网格点 (y 方向)
    
    for k = 0:(n_basis - 1)
        record = 0.0;
        for qx = 1:length(quad_point_x)
            for qy = 1:length(quad_point_y)
                phi = legendre_poly_2d(k, quad_point_x(qx), quad_point_y(qy));
                record = record + quad_weight_x(qx) * quad_weight_y(qy) * fcn_ic(x_quad(qx), y_quad(qy), v) * phi;
            end
        end
        local_rhs(k + 1) = record;
    end
    
    f = reference_mass \ local_rhs;
end

% 二维勒让德多项式函数
function val = legendre_poly_2d(k, x, y)
    if k == 0
        val = 1;
    elseif k == 1
        val = x;
    elseif k == 2
        val = y;
    else
        error('Unsupported polynomial order');
    end
end