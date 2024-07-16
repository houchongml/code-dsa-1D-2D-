function [flux_right_edge_plus flux_right_edge_minus flux_left_edge_plus flux_left_edge_minus]=generate_flux(poly_order)
	for k_row=0:poly_order
		for k_col=0:poly_order
			flux_right_edge_plus(k_row+1,k_col+1) = legendre_poly(k_row,1.)*legendre_poly(k_col,-1.);
			flux_right_edge_minus(k_row+1,k_col+1) = legendre_poly(k_row,1)*legendre_poly(k_col,1.);
			flux_left_edge_plus(k_row+1,k_col+1) = legendre_poly(k_row,-1)*legendre_poly(k_col,-1);
			flux_left_edge_minus(k_row+1,k_col+1) = legendre_poly(k_row,-1)*legendre_poly(k_col,1);
		end % k_col
	end % k_row
end
