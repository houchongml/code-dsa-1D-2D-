function mat_cross_section = impose_cross_section(x,dx,...
												quad_point_x,quad_weight_x,...
												N_x,poly_order,...
												cross_section_fun)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the matrix for the discretization of (\sigma_s\lgl f \rgl, \phi)
%
% Consider isotropic scattering for now
% 
% input: 
%	x: center of each element
%	dx: length of each element 
%	quad_point_x: quadrature points on the reference element
%	quad_weight_x: quadrauter weights on the reference element
%	N_x_tot: total number of element
%	poly_order: polynomial order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..');

% allocate memory
N_quadrature = length(quad_point_x);
N_tot_dof = N_x*(poly_order+1);
mat_cross_section = sparse(N_tot_dof,N_tot_dof);

% assign values
for i_row = 1:N_x
	for k_row = 0:poly_order
		ind_row = index_map_space(i_row,k_row,poly_order);
		for k_col = 0:poly_order
			ind_col = index_map_space(i_row,k_col,poly_order);
			% use quadrature to calcualte (\sigma_s phi_{k_col},phi_{k_row})
			record = 0.;
			for quad_count = 1:N_quadrature
				xx = x(i_row) + 0.5*dx(i_row) * quad_point_x(quad_count);%对应网格上的高斯点
				record = record + 0.5*dx(i_row)*quad_weight_x(quad_count)...
								  *cross_section_fun(xx)...
								  *legendre_poly( k_row,quad_point_x(quad_count) )...
								  *legendre_poly( k_col,quad_point_x(quad_count) );
			end
			mat_cross_section(ind_row,ind_col) = record;
		end % k_col
	end % k_row
end % i_row
