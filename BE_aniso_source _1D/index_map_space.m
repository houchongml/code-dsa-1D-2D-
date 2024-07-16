function ind = index_map_space(i,k,poly_order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suppose f(x)|_{I_i} = \sum f(i,k) phi_k(x)
% Ordering (f(1,0),..., f(1,k),..., f(N_x,0), ..., f(N_x,k) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = (i-1)*(poly_order+1) + k+1;

end
