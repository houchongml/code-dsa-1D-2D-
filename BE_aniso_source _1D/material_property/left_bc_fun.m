%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = left_bc_fun(v,time)
    global x_l
	res = f_exact(x_l,v,time);
end

