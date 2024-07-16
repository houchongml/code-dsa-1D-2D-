function val = legendre_poly_2d(k, x, y)
    if k == 0
        val = 1;
    elseif k == 1
        val = x;
    elseif k == 2
        val = y;
    else
        error('Unsupported polynomial order');
    end
end