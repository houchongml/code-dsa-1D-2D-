function ind = index_map_space_2d(i, j, k,N_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suppose f(x,y)|_{I_{i,j}} = \sum f(i,j,k) phi_k(x,y)
% Ordering (f(1,1,0), f(1,1,1), f(1,1,2), ..., f(N_x,N_y,2))
% i: x方向上的单元索引
% j: y方向上的单元索引
% k: 基函数索引（0对应phi_0 = 1, 1对应phi_1 = x, 2对应phi_2 = y）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_basis = 3; % 基函数数量
ind = (i - 1) * n_basis * N_y + (j - 1) * n_basis + k + 1;
end