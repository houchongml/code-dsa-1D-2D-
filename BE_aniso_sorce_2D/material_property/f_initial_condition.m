function res = f_initial_condition(x, y, v)
    % 初始条件设为一个常数
    initial_value = 1.0;
    res = initial_value * ones(length(x), 1);
end