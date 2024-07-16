function rhs_vector = impose_source(x,dx,velocity_point,...
									quad_point_x,quad_weight_x,...
									N_x_tot,poly_order,...
                                    time,...
									source_fun)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate righthand side vetor (source,phi_i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocate memory
N_tot_dof = N_x_tot*(poly_order+1);
rhs_vector = zeros(N_tot_dof,1);
N_quadrature = size(quad_point_x);

% vector
for i_row = 1:N_x_tot
    xx = x(i_row) + 0.5*dx(i_row)*quad_point_x;
	for k_row = 0:poly_order
		ind = index_map_space(i_row,k_row,poly_order);
		record = 0.;        
        
		rhs_vector(ind) = 0.5*dx(i_row)*dot(quad_weight_x,...
                          source_fun(xx,velocity_point,time).*legendre_poly(k_row,quad_point_x));     
    end
end

end
