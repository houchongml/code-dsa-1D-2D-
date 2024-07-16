% generate quadrature rule on the reference element,
% with N_quadrauter points
function [quad_point quad_weight] = set_up_reference_quadrature_rule(N_quadrature)

quad_point = zeros(N_quadrature,1);
quad_weight = zeros(N_quadrature,1);
if(N_quadrature==5)

	quad_point(1) =-0.906179845938664; 
	quad_point(2) =-0.538469310105683;                
	quad_point(3) = 0.;
	quad_point(4) = 0.538469310105683;
	quad_point(5) = 0.906179845938664;

	quad_weight(1) = 0.236926885056189; 
	quad_weight(2) = 0.478628670499366;                
	quad_weight(3) = 0.568888888888889; 
	quad_weight(4) = 0.478628670499366;
	quad_weight(5) = 0.236926885056189; 

else

end

end
