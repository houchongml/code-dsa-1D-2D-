function mat_cross_section = impose_cross_section_2d(x, y, dx, dy, ...
                                                     quad_point_x, quad_weight_x, ...
                                                     quad_point_y, quad_weight_y, ...
                                                     N_x, N_y, cross_section_fun)

% 基函数数量为3，对应基函数为1, x, y
n_basis = 3;
N_tot_dof = N_x * N_y * n_basis;

% 分配内存
N_quadrature_x = length(quad_point_x);
N_quadrature_y = length(quad_point_y);
mat_cross_section = sparse(N_tot_dof, N_tot_dof);

% 分配值
for i = 1:N_x
    for j = 1:N_y
        for k_row = 0:(n_basis - 1) % 0代表基函数1，1代表基函数x，2代表基函数y
            ind_row = index_map_space_2d(i, j, k_row, N_y);
            for k_col = 0:(n_basis - 1)
                ind_col = index_map_space_2d(i, j, k_col, N_y);
                % 使用高斯求积计算 (\sigma_s phi_{k_col}, phi_{k_row})
                record = 0.0;
                for quad_count_x = 1:N_quadrature_x
                    for quad_count_y = 1:N_quadrature_y
                        xx = x(i) + 0.5 * dx(i) * quad_point_x(quad_count_x);
                        yy = y(j) + 0.5 * dy(j) * quad_point_y(quad_count_y);
                        phi_row = legendre_poly_2d(k_row, quad_point_x(quad_count_x), quad_point_y(quad_count_y));
                        phi_col = legendre_poly_2d(k_col, quad_point_x(quad_count_x), quad_point_y(quad_count_y));
                        record = record + 0.25 * dx(i) * dy(j) * quad_weight_x(quad_count_x) * quad_weight_y(quad_count_y) ...
                                          * cross_section_fun(xx, yy) ...
                                          * phi_row ...
                                          * phi_col;
                    end
                end
                mat_cross_section(ind_row, ind_col) = record;
            end % k_col
        end % k_row
    end % j
end % i

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