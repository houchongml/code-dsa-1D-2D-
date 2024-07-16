% function [flux_right_edge_plus, flux_right_edge_minus, flux_left_edge_plus, flux_left_edge_minus, ...
%           flux_top_edge_plus, flux_top_edge_minus, flux_bottom_edge_plus, flux_bottom_edge_minus] = generate_flux_2d()
%     % 基函数数量
%     num_basis = 3; % 对应 1, x, y
%     
%     % 初始化flux矩阵
%     flux_right_edge_plus = cell(num_basis, num_basis);
%     flux_right_edge_minus = cell(num_basis, num_basis);
%     flux_left_edge_plus = cell(num_basis, num_basis);
%     flux_left_edge_minus = cell(num_basis, num_basis);
%     flux_top_edge_plus = cell(num_basis, num_basis);
%     flux_top_edge_minus = cell(num_basis, num_basis);
%     flux_bottom_edge_plus = cell(num_basis, num_basis);
%     flux_bottom_edge_minus = cell(num_basis, num_basis);
%     
%     % 右边界的 x 坐标为 1
%     % 左边界的 x 坐标为 -1
%     x_right = 1;
%     x_left = -1;
%     
%     % 上边界的 y 坐标为 1
%     % 下边界的 y 坐标为 -1
%     y_top = 1;
%     y_bottom = -1;
%     
%     % 遍历所有基函数对，计算左右边界flux
%     for k_row = 1:num_basis
%         for k_col = 1:num_basis
%             flux_right_edge_plus(k_row, k_col) = basis_function_value(k_row, x_right, 'x') * basis_function_value(k_col, x_left, 'x');
%             flux_right_edge_minus(k_row, k_col) = basis_function_value(k_row, x_right, 'x') * basis_function_value(k_col, x_right, 'x');
%             flux_left_edge_plus(k_row, k_col) = basis_function_value(k_row, x_left, 'x') * basis_function_value(k_col, x_left, 'x');
%             flux_left_edge_minus(k_row, k_col) = basis_function_value(k_row, x_left, 'x') * basis_function_value(k_col, x_right, 'x');
%             
%             flux_top_edge_plus(k_row, k_col) = basis_function_value(k_row, y_top, 'y') * basis_function_value(k_col, y_bottom, 'y');
%             flux_top_edge_minus(k_row, k_col) = basis_function_value(k_row, y_top, 'y') * basis_function_value(k_col, y_top, 'y');
%             flux_bottom_edge_plus(k_row, k_col) = basis_function_value(k_row, y_bottom, 'y') * basis_function_value(k_col, y_bottom, 'y');
%             flux_bottom_edge_minus(k_row, k_col) = basis_function_value(k_row, y_bottom, 'y') * basis_function_value(k_col, y_top, 'y');
%         end
%     end
% end
% 
% function val = basis_function_value(index, coord, direction)
%     % 根据索引和方向返回基函数值
%     switch index
%         case 1
%             val = 1;  % 基函数 1
%         case 2
%             if strcmp(direction, 'x')
%                 val = coord;  % 基函数 x 在 x 方向
%             else
%                 val = 'x';  % 基函数 x 在 y 方向，恒为 1
%             end
%         case 3
%             if strcmp(direction, 'y')
%                 val = coord;  % 基函数 y 在 y 方向
%             else
%                 val = 'y';  % 基函数 y 在 x 方向，恒为 1
%             end
%         otherwise
%             error('Invalid index');
%     end
% end
function [flux_right_edge_plus, flux_right_edge_minus, flux_left_edge_plus, flux_left_edge_minus, ...
          flux_top_edge_plus, flux_top_edge_minus, flux_bottom_edge_plus, flux_bottom_edge_minus] = generate_flux_2d()
    % 基函数数量
    num_basis = 3; % 对应 1, x, y
    
    % 初始化flux矩阵
    flux_right_edge_plus = zeros(num_basis, num_basis);
    flux_right_edge_minus = zeros(num_basis, num_basis);
    flux_left_edge_plus = zeros(num_basis, num_basis);
    flux_left_edge_minus = zeros(num_basis, num_basis);
    flux_top_edge_plus = zeros(num_basis, num_basis);
    flux_top_edge_minus = zeros(num_basis, num_basis);
    flux_bottom_edge_plus = zeros(num_basis, num_basis);
    flux_bottom_edge_minus = zeros(num_basis, num_basis);
    
    % 符号变量
    syms x y;
    
    % 右边界的 x 坐标为 1
    % 左边界的 x 坐标为 -1
    x_right = 1;
    x_left = -1;
    
    % 上边界的 y 坐标为 1
    % 下边界的 y 坐标为 -1
    y_top = 1;
    y_bottom = -1;
    
    % 遍历所有基函数对，计算左右边界flux
    for k_row = 1:num_basis
        for k_col = 1:num_basis
            % 右边界：右侧为外部单元
            flux_right_edge_plus(k_row, k_col) = double(int(basis_function_value(k_row, x_right, 'x') * basis_function_value(k_col, x_left, 'x'), y, -1, 1));
            flux_right_edge_minus(k_row, k_col) = double(int(basis_function_value(k_row, x_right, 'x') * basis_function_value(k_col, x_right, 'x'), y, -1, 1));
            % 左边界：左侧为外部单元
            flux_left_edge_plus(k_row, k_col) = double(int(basis_function_value(k_row, x_left, 'x') * basis_function_value(k_col, x_right, 'x'), y, -1, 1));
            flux_left_edge_minus(k_row, k_col) = double(int(basis_function_value(k_row, x_left, 'x') * basis_function_value(k_col, x_left, 'x'), y, -1, 1));
            
            % 上边界：上侧为外部单元
            flux_top_edge_plus(k_row, k_col) = double(int(basis_function_value(k_row, y_top, 'y') * basis_function_value(k_col, y_bottom, 'y'), x, -1, 1));
            flux_top_edge_minus(k_row, k_col) = double(int(basis_function_value(k_row, y_top, 'y') * basis_function_value(k_col, y_top, 'y'), x, -1, 1));
            % 下边界：下侧为外部单元
            flux_bottom_edge_plus(k_row, k_col) = double(int(basis_function_value(k_row, y_bottom, 'y') * basis_function_value(k_col, y_top, 'y'), x, -1, 1));
            flux_bottom_edge_minus(k_row, k_col) = double(int(basis_function_value(k_row, y_bottom, 'y') * basis_function_value(k_col, y_bottom, 'y'), x, -1, 1));
        end
    end
end

function val = basis_function_value(index, coord, direction)
    % 根据索引和方向返回基函数值
    syms x y;
    switch index
        case 1
            val = 1;  % 基函数 1
        case 2
            if strcmp(direction, 'x')
                val = coord;  % 基函数 x 在 x 方向
            else
                val = x;  % 基函数 x 在 y 方向，还是x
            end
        case 3
            if strcmp(direction, 'y')
                val = coord;  % 基函数 y 在 y 方向
            else
                val = y;  % 基函数 y 在 x 方向，还是y
            end
        otherwise
            error('Invalid index');
    end
end