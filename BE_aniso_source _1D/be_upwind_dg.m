%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/Users/houchong/Documents/python文件/输运方程/code_dsa_peng/1d_code/unscaled/BE_aniso_source copy/spatial_discretization');
addpath('/Users/houchong/Documents/python文件/输运方程/code_dsa_peng/1d_code/unscaled/BE_aniso_source copy/material_property')
addpath('/Users/houchong/Documents/python文件/输运方程/code_dsa_peng/1d_code/unscaled/BE_aniso_source copy/utility')

clear;clc;
global N
N = 80;
time_final = 0.1;
Nv = 8;

fprintf("N = %d\n",N);

global poly_order
poly_order = 1;
apply_dsa =  true;
%apply_dsa = false;

test_accuracy = false;

global x_center
global dx
global dt
global Nv


global local_dof %单个单元多项式维数
global spatial_dof %总的多项式维数
global vel_point 
global vel_weight
global reference_mass %基函数的求积积分
global reference_gradient
global Dminus_left_cell
global Dplus_right_cell
global Dminus_current_cell
global Dplus_current_cell
global left_edge_basis_value
global right_edge_basis_value

global x_r
global x_l

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: set up mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f_difference_tol = 1e-10;%1e-10;
max_iter = 100;

local_dof = poly_order+1;
spatial_dof = local_dof*N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_l = 0.;
x_r = 2*pi;

bc_conditions.f_l = zeros(Nv,1);
bc_conditions.f_r = zeros(Nv,1);


x_edge = linspace(x_l,x_r,N+1);
dx_min = (x_r-x_l)/N;
dx =zeros(N,1);
for i=1:N
 	dx(i) = x_edge(i+1)-x_edge(i);
end
x_center = 0.5*(x_edge(1:N)+x_edge(2:N+1));

% generate quadrature points and quadrature weights on the reference
% element [-1,1], sum(quad_weight)=1
N_quad = 5;
[quad_point,quad_weight] = legpts(N_quad);
fprintf("sum=%f",sum(quad_weight));
%quad_weight = 0.5 .* quad_weight ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: angular discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vel_point,vel_weight] = legpts(Nv);
fprintf("sum=%d", sum(vel_weight));
%vel_weight = 0.5 .* vel_weight ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: assemble DG spatial discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass matrix 
[reference_mass reference_gradient] = set_up_reference_matrices(poly_order); 
 reference_mass_inv = inv(reference_mass);

% Upwind term, matrices for fluxes
[flux_right_edge_plus...
 flux_right_edge_minus... 
 flux_left_edge_plus... 
 flux_left_edge_minus] = generate_flux(poly_order);

% Local matrix for D_minus 
Dminus_current_cell = flux_right_edge_minus-reference_gradient;
Dminus_left_cell = -flux_left_edge_minus;

% Local matrix for D_plus 
Dplus_current_cell = -flux_left_edge_plus-reference_gradient;
Dplus_right_cell = flux_right_edge_plus;


% Matrix for the total cross section 
 mat_total = impose_cross_section( ...
                        x_center,dx, ...
                        quad_point,quad_weight, ...
                        N,poly_order,...
                        @sigma_t);

 % Matrix for the scattering cross section
 mat_scattering = impose_cross_section( ...
                        x_center,dx, ...
                        quad_point,quad_weight, ...
                        N,poly_order,...
                        @sigma_s);

 % Left edge basis value
 left_edge_basis_value = zeros(local_dof,1);
 for k = 1:local_dof
    left_edge_basis_value(k) = legendre_poly(k-1,-1);
 end

 % right edge basis value
 right_edge_basis_value = zeros(local_dof,1);
 for k = 1:local_dof
    right_edge_basis_value(k) = legendre_poly(k-1,1);
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_dof = local_dof*N;
f = zeros(Nv*spatial_dof,1);%要求函数
rho = zeros(spatial_dof,1);%速度积分\phi

for j = 1:Nv
    ind_start = (j-1)*spatial_dof;
    for i = 1:N%N网格剖分
        f(ind_start+(i-1)*local_dof+1:ind_start+i*local_dof) = initialization(...
                            @f_initial_condition, ...
                            reference_mass, ...
                            poly_order,...
                            quad_point,quad_weight,...
                            vel_point(j),...
                            x_center(i),dx(i)...
                            );
    end
    ind = ind_start+(1:spatial_dof);
    rho = rho+0.5*vel_weight(j)*f(ind);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Assemble RHS due to the time indendent source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble righthand side due to the source
source_vec = zeros(Nv*spatial_dof,1);
f_old = f;
rho_old = rho;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 6: set up timestep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (test_accuracy)
    dt = dx_min^2;
else
    dt = dx_min;
end
%dt = 1e-3;
dt_eps_ratio = dt
N_time_step = ceil(time_final/dt);
%dt = time_final/N_time_step;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 7: set up preconditioner and allocate memory
% for SISA algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
local_rhs = zeros(local_dof,1);
% construct diffusion operator for the DSA
if (apply_dsa)
    % mass matrix
    mat_mass = assemble_mass(N,poly_order,dx,reference_mass);
    diffusion_operator = DSA_diffusion_operator(...
                               dt,...
							   x_center,dx,...
							   N,poly_order,...
                               mat_mass,...
                               mat_total,...
                               mat_scattering...
							   );
else
    diffusion_operator = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 7: time marching with backward Euler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = 0.0;
for time_step = 1:3
%for time_step = 1:N_time_step
    f_old = f;
    rho_old = rho;

    % impose boundary conditions
    time = time_step*dt;

    for j = 1:Nv
        if (vel_point(j)>0)
            bc_conditions.f_l(j) = left_bc_fun(vel_point(j),time);
        else
            bc_conditions.f_r(j) = right_bc_fun(vel_point(j),time);
        end
    end

    % impose source
    for j = 1:Nv
        ind = (j-1)*spatial_dof+(1:spatial_dof);
        source_vec(ind) = impose_source(x_center,dx,vel_point(j),...
									quad_point,quad_weight,...
									N,poly_order,...
                                    time,...
									@source_function);
    end

    [f,rho,iter,f_difference] = my_source_iteration(...
                        max_iter,...
                        f_difference_tol,...
                        f,f_old,...
                        local_rhs,...
                        mat_scattering,mat_total,...
                        rho,...
                        source_vec,...
                        bc_conditions,...
                        diffusion_operator,...
                        apply_dsa);

    %if (mod(time_step,10)==0)
    fprintf("time=%f, iter=%d, f_diff=%e \n",time,iter,f_difference);
    
    % If you only want to run the loop once, add a break statement here
    %break;
   
    %if(test_case == 10)
    %    fprintf("time=%f, iter=%d, f_diff=%e \n",time,iter,f_difference);
    %end
    %{
    clf
    figure(10)
    plot(x_center,f(1:local_dof:spatial_dof),'bo','Linewidth',1.5);
    hold on 
    plot(x_center,f_initial_condition(x_center)*exp(-1/3*time),'r-','Linewidth',1.5);
    %plot(x_center,f_initial_condition(x_center-vel_point(1)*time),'Linewidth',1.5)
    xlim([x_l,x_r])
    %vylim([-1,1])
    grid on
    pause(0.5)
    %}
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 8: present results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
% points_per_element = 5;
% f_error_inf = 0.0;
% f_error_l2  = 0.0;
% v_max = 0;
% for j = 1:Nv
%     ind_start = (j-1)*spatial_dof;
%     ind = ind_start+(1:spatial_dof);
%     [f_evaluate,x_evaluate] = evaluate_dg_function(f(ind),points_per_element);
%     f_ex = f_exact(x_evaluate,vel_point(j),time);
% 
%     figure(10*j)
%     plot(x_evaluate,f_evaluate,'bo')
%     xlim([x_l,x_r]);
%     set(gca,'FontSize',15);
%     hold on
%     plot(x_evaluate,f_ex,'r','LineWidth',1.5);
% 
%     f_error = abs(f_evaluate-f_ex);
%     f_error_l2_local = norm(f_error)*sqrt(dx_min/points_per_element);
%     fprintf("Error = %e at %d with weight %e\n",f_error_l2_local,j,vel_weight(j));
%     if (max(f_error)>f_error_inf)
%         f_error_inf = max(max(f_error),f_error_inf);
%         v_max = j;
%     end
%     f_error_l2 = max(f_error_l2,f_error_l2_local);
% end
% 
% [rho_evaluate,x_evaluate] = evaluate_dg_function(rho,points_per_element);
% rho_ex = sin(x_evaluate)*exp(-1/3*time);
% 
% rho_ee = zeros(size(rho_evaluate));
% for j = 1:Nv
%     ind = (j-1)*spatial_dof+1:j*spatial_dof;
%     [f_evaluate,x_evaluate] = evaluate_dg_function(f(ind),points_per_element);
%     rho_ee = rho_ee+0.5*vel_weight(j)*f_evaluate;
% end
% 
% 
% 
% figure(1)
% plot(x_evaluate,rho_evaluate,'bo');
% xlim([x_l,x_r]);
% set(gca,'FontSize',15);
% hold on
% plot(x_evaluate,rho_ex,'r','LineWidth',1.5);
% 
% rho_error = abs(rho_evaluate-rho_ex);
% rho_error_inf = max(rho_error);
% rho_error_l2 = norm(rho_error)*sqrt(dx_min/points_per_element);
% fprintf("dx = %e, dt = %e\n",dx_min,dt);
% fprintf("rho max error = %e\n",rho_error_inf);
% fprintf("rho l2 error = %e\n",rho_error_l2);
% fprintf("f max error = %e at %d\n",f_error_inf,v_max);
% fprintf("f l2 error = %e at %d\n",f_error_l2,v_max);
% fprintf("\n")
% 
% figure(10086)
% plot(x_evaluate,rho_ee-rho_evaluate,'o');
%     
