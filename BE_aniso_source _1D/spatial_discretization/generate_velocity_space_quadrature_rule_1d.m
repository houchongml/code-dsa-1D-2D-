% Given points \pm v_1, \pm v_2, ..., \pm v_k, where v_j \neq 0, forall j
% generate a quadrauture rule accurate enough for \int_{-1}^1 v^j dv, where j = 1,2,...,2k-1
% 
% input: 
%	velocity: descending order, velocity(1) > velocity(2)> ... > velocity(2k) = -velocity(1)
% output:
%	weight: corresponding weights in quadrautre rule
% Note: 
%	the result is normalized quadrature rule in [-1,1], \sum_i weight(i) =1 instead of 2
%

function [weight] = generate_velocity_space_quadrature_rule_1d(velocity)

addpath('..');

% set_up 
	N_v = length(velocity);
	k = N_v/2;
% assemble matrix
	Vandermonde_mat =zeros(k,k);
	for i=1:k
		for j=1:k
			Vandermonde_mat(i,j) = velocity(j)^(2*i-2);
		end
	end

% assemble vector
	rhs_vec = zeros(k,1);
	for i=1:k
		rhs_vec(i) = 1/(2*i-1);
	end

% obtain quadrature rule
	half_weight = zeros(k,1);
	half_weight = Vandermonde_mat\rhs_vec;
	weight = zeros(2*k,1);
	for i=1:k	
		weight(i)   = half_weight(i)/2;
		weight(2*k+1-i) = half_weight(i)/2;
	end


end
