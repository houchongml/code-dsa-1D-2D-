function f = transport_sweep(f,f_last_time_step,...
                    local_rhs,...
                    mat_scattering,mat_total,...
                    rho,...
                    source_vec,...
                    bc_conditions)
global dx
global dt
global Nv
global N
global local_dof
global spatial_dof
global vel_point 
global epsilon
global reference_mass
global Dminus_left_cell
global Dplus_right_cell
global Dminus_current_cell
global Dplus_current_cell
global left_edge_basis_value
global right_edge_basis_value

for j = 1:Nv
    ind_start = (j-1)*spatial_dof;
    if (vel_point(j)>0)
        for i = 1:N
            local_ind = (i-1)*local_dof+1:i*local_dof;
            ind = ind_start+local_ind;
            if (i==1)
                mflux_left = -left_edge_basis_value*bc_conditions.f_l(j);
            else
                mflux_left = Dminus_left_cell*f(ind-local_dof);
            end
            local_rhs=0.5*dx(i)*reference_mass*f_last_time_step(ind)...
                      +dt*( -vel_point(j)*mflux_left...
                            +mat_scattering(local_ind,local_ind)*rho(local_ind)...
                            +source_vec(ind) );
            f(ind) = ( 0.5*dx(i)*reference_mass...
                      +dt*(mat_total(local_ind,local_ind)...
                          +vel_point(j)*Dminus_current_cell)...
                      )\local_rhs;
        end
    else
        for i = N:-1:1
            local_ind = (i-1)*local_dof+1:i*local_dof;
            ind = ind_start+local_ind;
            if (i==N)
                flux_right = right_edge_basis_value*bc_conditions.f_r(j);
            else
                flux_right = Dplus_right_cell*f(ind+local_dof);
            end
            local_rhs=0.5*dx(i)*reference_mass*f_last_time_step(ind)...
                         +dt*( -vel_point(j)*flux_right...
                               +mat_scattering(local_ind,local_ind)*rho(local_ind)...
                               +source_vec(ind) );
            f(ind) = ( 0.5*dx(i)*reference_mass...
                       +dt*(mat_total(local_ind,local_ind)....
                           +vel_point(j)*Dplus_current_cell)...
                           )\local_rhs;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end for the function