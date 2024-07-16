function res = f_initial_condition(x,v)
res = zeros(length(x),1);
res(:) = f_exact(x,v,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end