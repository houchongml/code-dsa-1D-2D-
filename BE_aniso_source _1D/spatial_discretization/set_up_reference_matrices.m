% generate mass and stiff matrices on the reference element
function [reference_mass reference_gradient] = set_up_reference_matrices(poly_order)

addpath('..');

reference_mass = zeros(poly_order+1,poly_order+1);
reference_gradient = zeros(poly_order+1,poly_order+1);

%%%%%%%%%%%%%%%%%%%%%%
% mass matrix on the reference element
for k=1:poly_order+1
	reference_mass(k,k) = 2./(2*k-1);
end

%%%%%%%%%%%%%%%%%%%%%%
% stiff matrices (\phi_j,\grad phi_i) in 1d
if(poly_order>=1)
	reference_gradient(2,1) = 2.;	
end

if(poly_order>=2)
	reference_gradient(3,2) = 2.;
end


end
