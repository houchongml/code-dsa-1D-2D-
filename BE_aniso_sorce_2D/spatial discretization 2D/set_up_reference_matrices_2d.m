function [reference_mass, reference_gradient_x, reference_gradient_y] = set_up_reference_matrices_2d()

% 基函数数量为3，对应基函数为1, x, y
n_basis = 3;

% 初始化质量矩阵和梯度矩阵
reference_mass = zeros(n_basis, n_basis);
reference_gradient_x = zeros(n_basis, n_basis);
reference_gradient_y = zeros(n_basis, n_basis);

% 质量矩阵的计算
% \int_{-1}^{1} \int_{-1}^{1} phi_i * phi_j dx dy
reference_mass(1, 1) = 4;     % \int 1 * 1 dx dy
reference_mass(2, 2) = 4/3;   % \int x * x dx dy
reference_mass(3, 3) = 4/3;   % \int y * y dx dy

% 梯度矩阵的计算
% x方向的梯度矩阵
% \int_{-1}^{1} \int_{-1}^{1} phi_i * d(phi_j)/dx dx dy
reference_gradient_x(1, 2) = 4;   % \int 1 * d(x)/dx dx dy


% y方向的梯度矩阵
% \int_{-1}^{1} \int_{-1}^{1} phi_i * d(phi_j)/dy dx dy
reference_gradient_y(1, 3) = 4;   % \int 1 * d(y)/dy dx dy


end