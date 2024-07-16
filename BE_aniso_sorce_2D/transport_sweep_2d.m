function f = transport_sweep_2d(f, f_last_time_step, ...
                    local_rhs, ...
                    mat_scattering, mat_total, ...
                    rho, ...
                    source_vec, ...
                    bc_conditions)
global dx dy

global Nv
global N_x N_y
global local_dof
global spatial_dof
global vel_points
global epsilon
global reference_mass
global Dminus_left_cell_x Dplus_right_cell_x
global Dminus_current_cell_x Dplus_current_cell_x
global Dminus_left_cell_y Dplus_right_cell_y
global Dminus_current_cell_y Dplus_current_cell_y
global left_edge_basis_value right_edge_basis_value
global top_edge_basis_value bottom_edge_basis_value

for j = 1:Nv
    ind_start = (j-1) * spatial_dof;
    mu = vel_points(1, j);
    eta = vel_points(2, j);
    
    if (mu > 0 && eta > 0)
        for i = 1:N_x
            for k = 1:N_y
                local_ind = index_map_space_2d(i, k, 0:2, N_y);
                ind = ind_start + local_ind;
                if (i == 1)
                    mflux_left_x = -left_edge_basis_value * bc_conditions.f_l(j);
                else
                    mflux_left_x = Dminus_left_cell_x * f(ind - 3 * N_y);
                end
                if (k == 1)
                    mflux_bottom_y = -bottom_edge_basis_value * bc_conditions.f_b(j);
                else
                    mflux_bottom_y = Dminus_left_cell_y * f(ind - local_dof);
                end
                local_rhs =  ( - mu * mflux_left_x - eta * mflux_bottom_y ...
                                   + mat_scattering(local_ind, local_ind) * rho(local_ind) ...
                                   + source_vec(ind) );
                f(ind) =  (mat_total(local_ind, local_ind) ...
                                 + mu * Dminus_current_cell_x + eta * Dminus_current_cell_y) ...
                           \ local_rhs;
            end
        end
        
    elseif (mu > 0 && eta < 0)
        for i = 1:N_x
            for k = N_y:-1:1
                local_ind = index_map_space_2d(i, k, 0:2, N_y);
                ind = ind_start + local_ind;
                if (i == 1)
                    mflux_left_x = -left_edge_basis_value * bc_conditions.f_l(j);
                else
                    mflux_left_x = Dminus_left_cell_x * f(ind - 3 * N_y);
                end
                if (k == N_y)
                    mflux_top_y = top_edge_basis_value * bc_conditions.f_t(j);
                else
                    mflux_top_y = Dplus_right_cell_y * f(ind + local_dof);
                end
                local_rhs =  ( - mu * mflux_left_x - eta * mflux_top_y ...
                                   + mat_scattering(local_ind, local_ind) * rho(local_ind) ...
                                   + source_vec(ind) );
                f(ind) = (mat_total(local_ind, local_ind) ...
                                 + mu * Dminus_current_cell_x + eta * Dplus_current_cell_y) ...
                           \ local_rhs;
            end
        end
        
    elseif (mu < 0 && eta > 0)
        for i = N_x:-1:1
            for k = 1:N_y
                local_ind = index_map_space_2d(i, k, 0:2, N_y);
                ind = ind_start + local_ind;
                if (i == N_x)
                    mflux_right_x = right_edge_basis_value * bc_conditions.f_r(j);
                else
                    mflux_right_x = Dplus_right_cell_x * f(ind + 3 * N_y);
                end
                if (k == 1)
                    mflux_bottom_y = -bottom_edge_basis_value * bc_conditions.f_b(j);
                else
                    mflux_bottom_y = Dminus_left_cell_y * f(ind - local_dof);
                end
                local_rhs =  ( - mu * mflux_right_x - eta * mflux_bottom_y ...
                                   + mat_scattering(local_ind, local_ind) * rho(local_ind) ...
                                   + source_vec(ind) );
                f(ind) =  (mat_total(local_ind, local_ind) + mu * Dplus_current_cell_x + eta * Dminus_current_cell_y) ...
                            \ local_rhs;
            end
        end
        
    elseif (mu < 0 && eta < 0)
        for i = N_x:-1:1
            for k = N_y:-1:1
                local_ind = index_map_space_2d(i, k, 0:2, N_y);
                ind = ind_start + local_ind;
                if (i == N_x)
                    mflux_right_x = right_edge_basis_value * bc_conditions.f_r(j);
                else
                    mflux_right_x = Dplus_right_cell_x * f(ind + 3 * N_y);
                end
                if (k == N_y)
                    mflux_top_y = top_edge_basis_value * bc_conditions.f_t(j);
                else
                    mflux_top_y = Dplus_right_cell_y * f(ind + local_dof);
                end
                local_rhs =  ( - mu * mflux_right_x - eta * mflux_top_y ...
                                   + mat_scattering(local_ind, local_ind) * rho(local_ind) ...
                                   + source_vec(ind) );
                f(ind) =  (mat_total(local_ind, local_ind) ...
                                 + mu * Dplus_current_cell_x + eta * Dplus_current_cell_y) ...
                           \ local_rhs;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end for the function