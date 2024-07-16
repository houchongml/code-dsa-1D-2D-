function res = right_bc_fun(v,time)
    global x_r
	res = f_exact(x_r,v,time);
end