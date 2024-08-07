function mat_D_plus_y = assemble_D_plus_2d_y(N_x, N_y, reference_gradient_y)
    n_basis = 3; % 基函数数量为3，对应基函数为1, x, y
    N_tot = N_x * N_y * n_basis;
    mat_D_plus_y = sparse(N_tot, N_tot);

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
                    % 计算 (phi_i, -\partial_y phi_j)
                    mat_D_plus_y(ind_row, ind_col) = -reference_gradient_y(k_row + 1, k_col + 1);
                    % flux for the bottom edge
                    mat_D_plus_y(ind_row, ind_col) = mat_D_plus_y(ind_row, ind_col) ...
                                                     - flux_bottom_edge_plus(k_row + 1, k_col + 1);
                    % flux for the top edge
                    if j_row < N_y
                        ind_col = index_map_space_2d(i_row, j_row + 1, k_col, N_y);
                        mat_D_plus_y(ind_row, ind_col) = mat_D_plus_y(ind_row, ind_col) ...
                                                         + flux_top_edge_plus(k_row + 1, k_col + 1);
                    end
                end % k_col
            end % k_row
        end % j_row
    end % i_row
end