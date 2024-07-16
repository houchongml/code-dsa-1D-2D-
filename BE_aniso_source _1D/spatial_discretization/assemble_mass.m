function mat_mass=assemble_mass(N_x,poly_order,dx,reference_mass)

	N_tot = N_x*(poly_order+1);
	mat_mass = sparse(N_tot,N_tot);

	for i_row = 1:N_x
		for k_row = 0:poly_order
			ind_row = index_map_space(i_row,k_row,poly_order);
			for k_col = 0:poly_order
				ind_col = index_map_space(i_row,k_col,poly_order);
				% calcualte (phi_i,-\partial_x phi_j)
				mat_mass(ind_row,ind_col) = 0.5*dx(i_row)*reference_mass(k_row+1,k_col+1);
			end % k_col
		end % k_row
	end % i_row


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

