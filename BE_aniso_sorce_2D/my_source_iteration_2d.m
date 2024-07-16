function [f, rho, iter, f_difference] = my_source_iteration_2d(...
                        max_iter, ...
                        f_difference_tol, ...
                        f, f_last_step, ...
                        local_rhs, ...
                        mat_scattering, mat_total, ...
                        rho, ...
                        source_vec, ...
                        bc_conditions, ...
                        diffusion_operator, ...
                        apply_dsa)

global dx dy
global dt
global Nv
global N_x N_y
global local_dof
global spatial_dof
global vel_points
global vel_weights
global epsilon
global reference_mass
global Dminus_left_cell_x Dplus_right_cell_x
global Dminus_current_cell_x Dplus_current_cell_x
global Dminus_left_cell_y Dplus_right_cell_y
global Dminus_current_cell_y Dplus_current_cell_y
global left_edge_basis_value right_edge_basis_value
global top_edge_basis_value bottom_edge_basis_value

% Step 1: source iteration for one step
local_rhs = zeros(local_dof, 1);
f = transport_sweep_2d(f, f_last_step, ...
                    local_rhs, ...
                    mat_scattering, mat_total, ...
                    rho, ...
                    source_vec, ...
                    bc_conditions);

% reconstruct rho
rho_old = rho;
rho(:) = 0.0;
for j = 1:Nv
    ind = (j-1) * spatial_dof + 1:j * spatial_dof;
    rho = rho + (1/2*pi) * vel_weights(j) * f(ind);
end

% Diffusion synthetic acceleration (DSA)
rho_difference = rho - rho_old;
if (apply_dsa)
    rho_correction = diffusion_operator \ (mat_scattering * rho_difference);
    rho = rho + rho_correction;
end

% Step 2: source iteration
f_difference = f_difference_tol + 1;
iter = 0;
for i = 1:max_iter
    if (f_difference < f_difference_tol)
        break;
    end

    rho_old = rho; 
    f_old = f;
    
    % Transport sweep
    f = transport_sweep_2d(f, f_last_step, ...
                    local_rhs, ...
                    mat_scattering, mat_total, ...
                    rho, ...
                    source_vec, ...
                    bc_conditions);

    % reconstruct rho
    rho(:) = 0.0;
    for j = 1:Nv
        ind = (j-1) * spatial_dof + 1:j * spatial_dof;
        rho = rho + (1/2*pi) * vel_weights(j) * f(ind);
    end

    % Diffusion synthetic acceleration (DSA)
    rho_difference = rho - rho_old;
    if (apply_dsa)
        rho_correction = diffusion_operator \ (mat_scattering * rho_difference);
        rho = rho + rho_correction;
    end

    % f difference
    f_difference = max(max(abs(f - f_old)));
    iter = iter + 1;
end

end % end of function