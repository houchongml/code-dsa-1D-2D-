function mat_upwind = assemble_local_upwind_matrix(N_x,poly_order,velocity_point,reference_gradient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemble upwind matrix corresponding to a particular velocity point 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_tot = N_x*(poly_order+1);
mat_upwind = zeros(N_tot,N_tot);

% generate matrices storing flux terms 
	[flux_right_edge_plus...
	 flux_right_edge_minus... 
	 flux_left_edge_plus... 
	 flux_left_edge_minus] = generate_flux(poly_order);
% loop
for i_row = 1:N_x
	for k_row = 0:poly_order
		ind_row = index_map_space(i_row,k_row,poly_order);
		for k_col = 0:poly_order
			ind_col = index_map_space(i_row,k_col,poly_order);
		%----------part 1: calculate -(phi_j,v\dot\phi_i)-------------------
			mat_upwind(ind_row,ind_col) = -velocity_point...
										  *reference_gradient(k_row+1,k_col+1); %precomputed
		%----------part 2: flux term---------------------------------------
			% positive velocity point
			if(velocity_point>0)
				mat_upwind(ind_row,ind_col) = mat_upwind(ind_row,ind_col)...
											+velocity_point*flux_right_edge_minus(k_row+1,k_col+1);
				if(i_row==1) % left boudnary, do nothing

				else % interior
					ind_col = index_map_space(i_row-1,k_col,poly_order);
					mat_upwind(ind_row,ind_col) = mat_upwind(ind_row,ind_col) ...
									-velocity_point*flux_left_edge_minus(k_row+1,k_col+1);
				end
			% negative velocity point
			else
				mat_upwind(ind_row,ind_col) = mat_upwind(ind_row,ind_col)...
											-velocity_point*flux_left_edge_plus(k_row+1,k_col+1);
				if(i_row == N_x) % right boudnary, do nothing

				else
					ind_col = index_map_space(i_row+1,k_col,poly_order);
					mat_upwind(ind_row,ind_col) = mat_upwind(ind_row,ind_col) ...
									+velocity_point*flux_right_edge_plus(k_row+1,k_col+1);
				end
			end % if
		end % k_col
	end % k_row
end % i_row

% return f^+(1)*phi(1), f^-*phi(1), f^+*phi(-1), f^-*phi(-1)
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




end
