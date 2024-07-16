function rhs_vector = impose_source_2d(x, y, dx, dy, velocity_point, ...
                                       quad_point_x, quad_weight_x, ...
                                       quad_point_y, quad_weight_y, ...
                                       N_x, N_y, ...
                                       source_fun)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 生成右端项向量 (source,phi_i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 分配内存
n_basis = 3; % 基函数数量，对应于 1, x, y
N_tot_dof = N_x * N_y * n_basis;
rhs_vector = zeros(N_tot_dof, 1);

% 各向量
for i_row = 1:N_x
    for j_row = 1:N_y
        xx = x(i_row) + 0.5 * dx(i_row) * quad_point_x;
        yy = y(j_row) + 0.5 * dy(j_row) * quad_point_y;
        for k_row = 0:(n_basis - 1)
            ind = index_map_space_2d(i_row, j_row, k_row, N_y);
            record = 0.0;
            for qx = 1:length(quad_point_x)
                for qy = 1:length(quad_point_y)
                    phi = legendre_poly_2d(k_row, quad_point_x(qx), quad_point_y(qy));
                    record = record + quad_weight_x(qx) * quad_weight_y(qy) * source_fun(xx(qx), yy(qy), velocity_point) * phi;
                end
            end
            rhs_vector(ind) = 0.25 * dx(i_row) * dy(j_row) * record;
        end
    end
end

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

% 二维索引映射函数
function ind = index_map_space_2d(i, j, k, N_y)
n_basis = 3; % 基函数数量
ind = (i - 1) * n_basis * N_y + (j - 1) * n_basis + k + 1;
end