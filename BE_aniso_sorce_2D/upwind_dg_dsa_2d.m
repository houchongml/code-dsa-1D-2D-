clear;clc;
addpath('/Users/houchong/Documents/python文件/输运方程/code_dsa_peng/1d_code/unscaled/BE_aniso_sorce_2D/spatial discretization 2D')
addpath('/Users/houchong/Documents/python文件/输运方程/code_dsa_peng/1d_code/unscaled/BE_aniso_sorce_2D/material_property')
global N_x N_y
N_x = 8;
N_y = 8;

Nv = 6;

fprintf("N_x = %d, N_y = %d\n", N_x, N_y);

global poly_order
poly_order = 1;
%apply_dsa = true;
apply_dsa = false;

test_accuracy = false;

global x_center y_center
global dx dy
global dt
global Nv

global local_dof %单个单元多项式维数
global spatial_dof %总的多项式维数
global vel_points 
global vel_weights
global reference_mass %基函数的求积积分
global reference_gradient_x reference_gradient_y

global Dminus_left_cell_x Dplus_right_cell_x
global Dminus_current_cell_x Dplus_current_cell_x
global Dminus_left_cell_y Dplus_right_cell_y
global Dminus_current_cell_y Dplus_current_cell_y
global left_edge_basis_value right_edge_basis_value
global top_edge_basis_value bottom_edge_basis_value




global x_r x_l
global y_r y_l
global bc_conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: set up mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_difference_tol = 1e-10;
max_iter = 100;

local_dof = 3; % 基函数为1, x, y
spatial_dof = local_dof * N_x * N_y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_l = 0.0;
x_r = 2 * pi;
y_l = 0.0;
y_r = 2 * pi;


bc_conditions.f_l = zeros(Nv, 1);
bc_conditions.f_r = zeros(Nv, 1);
bc_conditions.f_b = zeros(Nv, 1);
bc_conditions.f_t = zeros(Nv, 1);
% 设置x方向的网格
x_edge = linspace(x_l, x_r, N_x + 1);
dx_min = (x_r - x_l) / N_x;
dx = zeros(N_x, 1);
for i = 1:N_x
    dx(i) = x_edge(i + 1) - x_edge(i);
end
x_center = 0.5 * (x_edge(1:N_x) + x_edge(2:N_x + 1));

% 设置y方向的网格
y_edge = linspace(y_l, y_r, N_y + 1);
dy_min = (y_r - y_l) / N_y;
dy = zeros(N_y, 1);
for j = 1:N_y
    dy(j) = y_edge(j + 1) - y_edge(j);
end
y_center = 0.5 * (y_edge(1:N_y) + y_edge(2:N_y + 1));

% generate quadrature points and quadrature weights on the reference element [-1,1], sum(quad_weight)=1
N_quad = 5;
[quad_point_x, quad_weight_x] = legpts(N_quad);
[quad_point_y, quad_weight_y] = legpts(N_quad);
%fprintf("sum=%f", sum(quad_weight));

% 设置速度点和权重
[vel_points, vel_weights] = set_velocity_points(Nv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: angular discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [vel_points, vel_weights] = set_velocity_points(Nv)
%     % 使用均匀分布在[0, 2*pi]上的角度
%     theta = linspace(0, 2 * pi, Nv + 1);
%     theta(end) = []; % 去掉最后一个点，避免重复
%     
%     % 速度点
%     vel_points = [cos(theta); sin(theta)];
%     
%     % 权重
%     vel_weights = ones(1, Nv) * (2 * pi / Nv);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: assemble DG spatial discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass matrix
[reference_mass, reference_gradient_x, reference_gradient_y] = set_up_reference_matrices_2d(); 
reference_mass_inv = inv(reference_mass);

% Upwind term, matrices for fluxes
[flux_right_edge_plus, flux_right_edge_minus, ...
 flux_left_edge_plus, flux_left_edge_minus, ...
 flux_top_edge_plus, flux_top_edge_minus, ...
 flux_bottom_edge_plus, flux_bottom_edge_minus] = generate_flux_2d();

% Local matrix for D_minus in x direction
Dminus_current_cell_x = flux_right_edge_minus - reference_gradient_x;
Dminus_left_cell_x = -flux_left_edge_minus;

% Local matrix for D_plus in x direction
Dplus_current_cell_x = -flux_left_edge_plus - reference_gradient_x;
Dplus_right_cell_x = flux_right_edge_plus;

% Local matrix for D_minus in y direction
Dminus_current_cell_y = flux_top_edge_minus - reference_gradient_y;
Dminus_left_cell_y = -flux_bottom_edge_minus;

% Local matrix for D_plus in y direction
Dplus_current_cell_y = -flux_bottom_edge_plus - reference_gradient_y;
Dplus_right_cell_y = flux_top_edge_plus;
% Matrix for the total cross section
mat_total = impose_cross_section_2d( ...
                        x_center, y_center, dx, dy, ...
                        quad_point_x, quad_weight_x, ...
                        quad_point_y, quad_weight_y, ...
                        N_x, N_y, ...
                        @sigma_t_2d);

% Matrix for the scattering cross section
mat_scattering = impose_cross_section_2d( ...
                        x_center, y_center, dx, dy, ...
                        quad_point_x, quad_weight_x, ...
                        quad_point_y, quad_weight_y, ...
                        N_x, N_y, ...
                        @sigma_s_2d);

% Left edge basis value
    left_edge_basis_value = zeros(local_dof, 1);
    for k = 1:local_dof
        left_edge_basis_value(k) = legendre_poly_2d(k - 1, -1, 0);
    end

    % Right edge basis value
    right_edge_basis_value = zeros(local_dof, 1);
    for k = 1:local_dof
        right_edge_basis_value(k) = legendre_poly_2d(k - 1, 1, 0);
    end

    % Top edge basis value
    top_edge_basis_value = zeros(local_dof, 1);
    for k = 1:local_dof
        top_edge_basis_value(k) = legendre_poly_2d(k - 1, 0, 1);
    end

    % Bottom edge basis value
    bottom_edge_basis_value = zeros(local_dof, 1);
    for k = 1:local_dof
        bottom_edge_basis_value(k) = legendre_poly_2d(k - 1, 0, -1);
    end
% Step 4: Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_dof = local_dof * N_x * N_y;
f = zeros(Nv * spatial_dof, 1); % 要求函数
rho = zeros(spatial_dof, 1); % 速度积分\phi

for j = 1:Nv
    ind_start = (j-1) * spatial_dof;
    for i = 1:N_x % N_x网格剖分
        for k = 1:N_y % N_y网格剖分
            local_ind = (i-1) * N_y + k; % 计算二维网格索引
            f(ind_start + (local_ind-1) * local_dof + 1 : ind_start + local_ind * local_dof) = initialization_2d(...
                @f_initial_condition, ...
                reference_mass, ...
                quad_point_x, quad_weight_x, ...
                quad_point_y, quad_weight_y, ...
                vel_points(:, j), ...
                x_center(i), dx(i), ...
                y_center(k), dy(k));
        end
    end
    ind = ind_start + (1:spatial_dof);
    rho = rho + (1/2*pi) * vel_weights(j) * f(ind);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Assemble RHS due to the time-independent source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble righthand side due to the source
source_vec = zeros(Nv * spatial_dof, 1);
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
% dt = 1e-3;
% dt_eps_ratio = dt;
% N_time_step = ceil(time_final / dt);
% % dt = time_final / N_time_step;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 7: set up preconditioner and allocate memory
% for SISA algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
local_rhs = zeros(local_dof, 1);
% construct diffusion operator for the DSA
if (apply_dsa)
    % mass matrix
    %mat_mass = assemble_mass_2d(N_x, N_y, dx, dy, reference_mass);
    diffusion_operator = DSA_diffusion_operator_2d(...
                               dt, ...
                               x_center, y_center, dx, dy, ...
                               N_x, N_y, ...
                               mat_total, ...
                               mat_scattering ...
                               );
else
    diffusion_operator = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 8: Source Iteration (no time dependency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_old = f;
rho_old = rho;
% impose boundary conditions and source term
for j = 1:Nv
    if vel_points(1, j) > 0
        bc_conditions.f_l(j) = left_bc_fun(vel_points(:, j));
    else
        bc_conditions.f_r(j) = right_bc_fun(vel_points(:, j));
    end
    
    if vel_points(2, j) > 0
        bc_conditions.f_b(j) = bottom_bc_fun(vel_points(:, j));
    else
        bc_conditions.f_t(j) = top_bc_fun(vel_points(:, j));
    end
end

for j = 1:Nv
    ind = (j-1) * spatial_dof + (1:spatial_dof);
    source_vec(ind) = impose_source_2d(x_center, y_center, dx, dy, vel_points(:, j), ...
                                       quad_point_x, quad_weight_x, ...
                                       quad_point_y, quad_weight_y, ...
                                       N_x, N_y, ...
                                       @source_function_2d);
end

[f, rho, iter, f_difference] = my_source_iteration_2d(...
                    max_iter, ...
                    f_difference_tol, ...
                    f, f_old, ...
                    local_rhs, ...
                    mat_scattering, mat_total, ...
                    rho, ...
                    source_vec, ...
                    bc_conditions, ...
                    diffusion_operator, ...
                    apply_dsa);

fprintf("iter=%d, f_diff=%e\n", iter, f_difference);

% function [vel_points, vel_weights] = set_velocity_points(Nv)
%     % 使用均匀分布在[0, 2*pi]上的角度
%     theta = linspace(0, 2 * pi, Nv + 1);
%     theta(end) = []; % 去掉最后一个点，避免重复
%     
%     % 速度点
%     vel_points = [cos(theta); sin(theta)];
%     
%     % 权重
%     vel_weights = ones(1, Nv) * (2 * pi / Nv);
% end
function [vel_points, vel_weights] = set_velocity_points(Nv)
    % 使用 [-1, 1] 区间上的高斯点
    [theta, w] = legpts(Nv);
    
    % 将 [-1, 1] 区间上的高斯点映射到 [0, 2*pi)
    theta = (theta + 1) * pi; % 从 [-1, 1] 映射到 [0, 2*pi)
    w = w * pi; % 权重也需要进行相应调整

    % 速度点
    vel_points = [cos(theta)'; sin(theta)'];
    
    % 权重
    vel_weights = w';
end

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