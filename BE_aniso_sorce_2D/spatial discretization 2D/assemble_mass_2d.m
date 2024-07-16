function mat_mass = assemble_mass_2d(N_x, N_y, dx, dy, reference_mass)

% 基函数数量为3，对应基函数为1, x, y
n_basis = 3;
N_tot = N_x * N_y * n_basis;

% 初始化稀疏矩阵
mat_mass = sparse(N_tot, N_tot);

for i_row = 1:N_x
    for j_row = 1:N_y
        for k_row = 0:(n_basis - 1) % 0代表基函数1，1代表基函数x，2代表基函数y
            ind_row = index_map_space_2d(i_row, j_row, k_row, N_y);
            for k_col = 0:(n_basis - 1)
                ind_col = index_map_space_2d(i_row, j_row, k_col, N_y);
                % 计算 (phi_i, phi_j)
                mat_mass(ind_row, ind_col) = 0.25 * dx(i_row) * dy(j_row) * reference_mass(k_row + 1, k_col + 1);
            end % k_col
        end % k_row
    end % j_row
end % i_row

end

% 二维索引映射函数
function ind = index_map_space_2d(i, j, k, N_y)
n_basis = 3; % 基函数数量
ind = (i - 1) * n_basis * N_y + (j - 1) * n_basis + k + 1;
end