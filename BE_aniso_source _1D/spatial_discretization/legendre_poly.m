%% Legendre polynomial on reference element [-1,1]
function res=legendre_poly(k,x)
    if(k==0)
        res=ones(size(x));
    elseif (k==1)
        res=x;
    elseif (k==2)
        res=0.5*(3*x.^2-1);
    end

end