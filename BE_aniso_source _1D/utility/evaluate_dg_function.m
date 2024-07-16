function [res,x_evaluate] = evaluate_dg_function(u,points_per_element)

global N
global poly_order
global x_center
global dx
global local_dof

res = zeros(N*points_per_element,1);
x_evaluate = zeros(N*points_per_element,1);

x_point = linspace(-1,1,points_per_element+2);
evaluate_point = x_point(2:end-1);

ind = 1;
for i = 1:N
    for j = 1:points_per_element
        xx = evaluate_point(j); 
        x_evaluate(ind) = x_center(i)+0.5*dx(i)*xx;
        for k = 0:poly_order
            u_ind = index_map_space(i,k,poly_order);
            %u_ind = (i-1)*local_dof+k+1;
            res(ind) = res(ind)+u(u_ind)*legendre_poly(k,xx);
        end
        ind = ind+1;
    end
end

end