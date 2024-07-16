function res = source_function(x,v,time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% source_funciton(x) return the value of a source funciton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test 1
res = (-1/3*(sin(x)-v*cos(x))... % f_t
       +v^2*sin(x)... % f_x-(rho-f)/eps
      )*exp(-time/3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
