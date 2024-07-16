function res = f_exact(x,v,t)
    res = ( sin(x)-v*cos(x) )*exp(-1/3*t);
end