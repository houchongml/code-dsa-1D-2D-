function rhs_vector = impose_local_inflow_bc_vector(rhs_vector,velocity_point,inflow_bc,N_x_tot,poly_order)

addpath('..');

% impose inflow bc

if(velocity_point>0)
	% left bc
	i_row = 1;
	for k_row = 0:poly_order 
		ind = index_map_space(i_row,k_row,poly_order);
		rhs_vector(ind) = rhs_vector(ind)...
						 +velocity_point...
						 *inflow_bc...
						 *legendre_poly(k_row,-1.);
	end
else	
	% right bc
	i_row = N_x_tot;
	for k_row = 0:poly_order
		ind = index_map_space(i_row,k_row,poly_order);
		rhs_vector(ind) = rhs_vector(ind)...
						 -velocity_point...
						 *inflow_bc...
						 *legendre_poly(k_row,1.);
	end
end


end
