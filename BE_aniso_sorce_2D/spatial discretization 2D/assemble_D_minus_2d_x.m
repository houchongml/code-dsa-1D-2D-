function mat_D_minus_x = assemble_D_minus_2d_x(N_x, N_y, reference_gradient_x)
    n_basis = 3; % 基函数数量为3，对应基函数为1, x, y
    N_tot = N_x * N_y * n_basis;
    mat_D_minus_x = sparse(N_tot, N_tot);

    [flux_right_edge_plus, flux_right_edge_minus, ...
     flux_left_edge_plus, flux_left_edge_minus, ...
     flux_top_edge_plus, flux_top_edge_minus, ...
     flux_bottom_edge_plus, flux_bottom_edge_minus] = generate_flux_2d();

    for i_row = 1:N_x
        for j_row = 1:N_y
            for k_row = 0:(n_basis - 1)
                ind_row = index_map_space_2d(i_row, j_row, k_row, N_y);
                for k_col = 0:(n_basis - 1)
                    ind_col = index_map_space_2d(i_row, j_row, k_col, N_y);
                    % 计算 (phi_i, -\partial_x phi_j)
                    mat_D_minus_x(ind_row, ind_col) = -reference_gradient_x(k_row + 1, k_col + 1);
                    % flux for the left edge
                    if i_row > 1
                        ind_col = index_map_space_2d(i_row - 1, j_row, k_col, N_y);
                        mat_D_minus_x(ind_row, ind_col) = mat_D_minus_x(ind_row, ind_col) ...
                                                          - flux_left_edge_minus(k_row + 1, k_col + 1);
                    end
                    % flux for the right edge
                    ind_col = index_map_space_2d(i_row, j_row, k_col, N_y);
                    mat_D_minus_x(ind_row, ind_col) = mat_D_minus_x(ind_row, ind_col) ...
                                                      + flux_right_edge_minus(k_row + 1, k_col + 1);
                end % k_col
            end % k_row
        end % j_row
    end % i_row
end