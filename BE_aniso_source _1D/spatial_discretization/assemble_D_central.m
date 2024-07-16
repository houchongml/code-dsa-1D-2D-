function mat_D_central=assemble_D_central(N_x,poly_order,reference_gradient)

	N_tot = N_x*(poly_order+1);
	mat_D_central = sparse(N_tot,N_tot);

	[flux_right_edge_plus...
	 flux_right_edge_minus... 
	 flux_left_edge_plus... 
	 flux_left_edge_minus] = generate_flux(poly_order);

	for i_row = 1:N_x
		for k_row = 0:poly_order
			ind_row = index_map_space(i_row,k_row,poly_order);
			for k_col = 0:poly_order
				ind_col = index_map_space(i_row,k_col,poly_order);
				% calcualte (phi_i,-\partial_x phi_j)
				mat_D_central(ind_row,ind_col) = -reference_gradient(k_row+1,k_col+1);
				% flux for the left edge, use reflection to impose zero boundary condition
				if(i_row==1)

				else
					ind_col = index_map_space(i_row-1,k_col,poly_order);
					mat_D_central(ind_row,ind_col) = mat_D_central(ind_row,ind_col)...
											     -0.5*flux_left_edge_minus(k_row+1,k_col+1);
				end
				ind_col = index_map_space(i_row,k_col,poly_order);
				mat_D_central(ind_row,ind_col) = mat_D_central(ind_row,ind_col)...
										     -0.5*flux_left_edge_plus(k_row+1,k_col+1);
				% flux for the right edge
				if(i_row == N_x)

				else
					ind_col = index_map_space(i_row+1,k_col,poly_order);
					mat_D_central(ind_row,ind_col) = mat_D_central(ind_row,ind_col)...
											     +0.5*flux_right_edge_plus(k_row+1,k_col+1);
				end
				ind_col = index_map_space(i_row,k_col,poly_order);
				mat_D_central(ind_row,ind_col) = mat_D_central(ind_row,ind_col)...
										     +0.5*flux_right_edge_minus(k_row+1,k_col+1);
			end % k_col
		end % k_row
	end % i_row


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

