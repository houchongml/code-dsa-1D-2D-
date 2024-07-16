function mat_D_plus = assemble_D_plus_2d(N_x, N_y, reference_gradient_x, reference_gradient_y)

% 基函数数量为3，对应基函数为1, x, y
n_basis = 3;
N_tot = N_x * N_y * n_basis;

% 初始化稀疏矩阵
mat_D_plus = sparse(N_tot, N_tot);

% 生成流量矩阵
[flux_right_edge_plus, flux_right_edge_minus, ...
 flux_left_edge_plus, flux_left_edge_minus, ...
 flux_top_edge_plus, flux_top_edge_minus, ...
 flux_bottom_edge_plus, flux_bottom_edge_minus] = generate_flux();

for i = 1:N_x
    for j = 1:N_y
        for k_row = 0:(n_basis - 1) % 0代表基函数1，1代表基函数x，2代表基函数y
            ind_row = index_map_space_2d(i, j, k_row, N_y);
            for k_col = 0:(n_basis - 1)
                ind_col = index_map_space_2d(i, j, k_col, N_y);
                % 计算 (phi_i, -\partial_x phi_j)
                mat_D_plus(ind_row, ind_col) = -reference_gradient_x(k_row + 1, k_col + 1) - reference_gradient_y(k_row + 1, k_col + 1);
                
                % 左边界上的通量
                mat_D_plus(ind_row, ind_col) = mat_D_plus(ind_row, ind_col) - flux_left_edge_plus(k_row + 1, k_col + 1);
                
                % 右边界上的通量
                if i < N_x
                    ind_col_right = index_map_space_2d(i + 1, j, k_col, N_y);
                    mat_D_plus(ind_row, ind_col_right) = mat_D_plus(ind_row, ind_col_right) + flux_right_edge_plus(k_row + 1, k_col + 1);
                end
                
                % 上边界上的通量
                if j < N_y
                    ind_col_top = index_map_space_2d(i, j + 1, k_col, N_y);
                    mat_D_plus(ind_row, ind_col_top) = mat_D_plus(ind_row, ind_col_top) + flux_top_edge_plus(k_row + 1, k_col + 1);
                end
                
                % 下边界上的通量
                
                ind_col_bottom = index_map_space_2d(i, j , k_col, N_y);
                mat_D_plus(ind_row, ind_col_bottom) = mat_D_plus(ind_row, ind_col_bottom) + flux_bottom_edge_plus(k_row + 1, k_col + 1);
                
            end % k_col
        end % k_row
    end % j
end % i

end

% 二维索引映射函数
function ind = index_map_space_2d(i, j, k, N_y)
n_basis = 3; % 基函数数量
ind = (i - 1) * n_basis * N_y + (j - 1) * n_basis + k + 1;
end