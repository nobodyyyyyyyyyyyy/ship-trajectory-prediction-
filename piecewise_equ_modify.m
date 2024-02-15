%Author: Zhang Ziqing
%University: Wuhan University of Technology
%Contact Email: 322980@whut.edu.cn
%Feedback and Support
%If you have any questions, suggestions, or would like to contribute, feel free to contact me 
% or raise an issue on the GitHub Issues page.
function [reserve_x, reserve_y] = piecewise_equ_modify(flag, intersection, initial_point , exclusion , flag_1)
%     Clear command window and workspace
%     clc,clear
    % Define symbolic variables
    syms t x y 
    syms x_opposite
%     flag = 'D_1';
%     equation0 = (2^(1/2)*((11*2^(1/2))/2 + 10)*(t + x + 2))/11 - (3^(1/2)*(4*t - 1797))/2 - t/2 - 1;
%     time = 100;
%     intersection = (t/2 + (3^(1/2)*(4*t - 1797))/2 - (2^(1/2)*((11*2^(1/2))/2 + 10)*(t + 2))/11 + 1)/((11*2^(1/2))/(2*((11*2^(1/2))/2 + 10)) + (2^(1/2)*((11*2^(1/2))/2 + 10))/11);
%     initial_point = 2*t - (3^(1/2)*(t + 2))/2 - 1797/2;


%     Constants and calculations

global  r_fort r_aft r_starb r_port 
global equation_orthographic equation_combine
global  x_position_11 x_position_12 x_position_21 x_position_22

if flag == 'D_0'
    if ~ismember(exclusion, 1)
        formula_ = matlabFunction(equation_combine);
        parameter_11 = @(x) abs(power(2 * x ./ ((1 + 1) * r_starb - (1 - 1) * r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 + 1) * r_fort - (1 - 1) * r_aft), 2) - 1);
        lb = 0;
        ub = x_position_11;
        options = optimset( 'MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(1), feval(1)] = fminbnd(parameter_11, lb, ub, options);
    end


if ~ismember(exclusion , 2)
        formula_ = matlabFunction(equation_combine);
        parameter_11 = @(x) abs(power(2 * x ./ ((1 - 1) * r_starb - (1 + 1) * r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 + 1) * r_fort - (1 - 1) * r_aft), 2) - 1);
        lb = x_position_12;
        ub = 0;
        options = optimset( 'MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(2), feval(2)] = fminbnd(parameter_11, lb, ub, options);
end
if ~ismember(exclusion , 3)
        formula_ = matlabFunction(equation_combine);
        parameter_11 = @(x) abs(power(2 * x ./ ((1 - 1) * r_starb - (1 + 1) * r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 - 1) * r_fort - (1 + 1) * r_aft), 2) - 1);
        lb = x_position_12;
        ub = 0;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(3), feval(3)] = fminbnd(parameter_11, lb, ub, options);
end
if ~ismember(exclusion , 4)
        formula_ = matlabFunction(equation_combine);
        parameter_11 = @(x) abs(power(2 * x ./ ((1 + 1) * r_starb - (1 - 1) * r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 - 1) * r_fort - (1 + 1) * r_aft), 2) - 1);
        lb = 0;
        ub = x_position_11;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(4), feval(4)] = fminbnd(parameter_11, lb, ub, options);
end
non_zero_indices = find(feval ~= 0);
feval = feval(non_zero_indices);
non_zero_indices = find(x_min ~= 0);
reserve_x = x_min(non_zero_indices);
reserve_x = reserve_x(sign(reserve_x - initial_point) == flag_1);
differences = abs(reserve_x - initial_point);
[~, index] = min(differences);
reserve_x =reserve_x(index);
reserve_y = formula_(reserve_x);

    elseif flag == 'D_1'

    if ~ismember(exclusion, 1)
        formula_ = matlabFunction(equation_combine);
        parameter_21 = @(x) abs(power(2 * x ./ ((1 + 1) *2* r_starb - (1 - 1) *2* r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 + 1) *2* r_fort - (1 - 1) * 2*r_aft), 2) - 1);
        lb = 0;
        ub = x_position_21;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(1), feval(1)] = fminbnd(parameter_21, lb, ub, options);
    end


if ~ismember(exclusion , 2)
        formula_ = matlabFunction(equation_combine);
        parameter_21 = @(x) abs(power(2 * x ./ ((1 - 1) *2* r_starb - (1 + 1) *2* r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 + 1) *2* r_fort - (1 - 1) *2* r_aft), 2) - 1);
        lb = x_position_22;
        ub = 0;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(2), feval(2)] = fminbnd(parameter_21, lb, ub, options);
end
if ~ismember(exclusion , 3)
        formula_ = matlabFunction(equation_combine);
        parameter_21 = @(x) abs(power(2 * x ./ ((1 - 1) *2* r_starb - (1 + 1) *2* r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 - 1) *2* r_fort - (1 + 1) *2* r_aft), 2) - 1);
        lb = x_position_22;
        ub = 0;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(3), feval(3)] = fminbnd(parameter_21, lb, ub, options);
end
if ~ismember(exclusion , 4)
        formula_ = matlabFunction(equation_combine);
        parameter_21 = @(x) abs(power(2 * x ./ ((1 + 1) *2* r_starb - (1 - 1) *2* r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 - 1) *2* r_fort - (1 + 1) *2* r_aft), 2) - 1);
        lb = 0;
        ub = x_position_21;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(4), feval(4)] = fminbnd(parameter_21, lb, ub, options);
end
non_zero_indices = find(feval ~= 0);
feval = feval(non_zero_indices);
non_zero_indices = find(x_min ~= 0);
reserve_x = x_min(non_zero_indices);
reserve_x = reserve_x(sign(reserve_x - initial_point) == flag_1);
differences = abs(reserve_x - initial_point);
[~, index] = min(differences);
reserve_x =reserve_x(index);
reserve_y = formula_(reserve_x);

    elseif flag == 'R_1'
    if ~ismember(exclusion, 1)
        formula_ = matlabFunction(equation_orthographic);
        parameter_11 = @(x) abs(power(2 * x ./ ((1 + 1) * r_starb - (1 - 1) * r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 + 1) * r_fort - (1 - 1) * r_aft), 2) - 1);
        lb = 0;
        ub = x_position_11;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(1), feval(1)] = fminbnd(parameter_11, lb, ub, options);
    end


if ~ismember(exclusion , 2)
        formula_ = matlabFunction(equation_orthographic);
        parameter_11 = @(x) abs(power(2 * x ./ ((1 - 1) * r_starb - (1 + 1) * r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 + 1) * r_fort - (1 - 1) * r_aft), 2) - 1);
        lb = x_position_12;
        ub = 0;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(2), feval(2)] = fminbnd(parameter_11, lb, ub, options);
end
if ~ismember(exclusion , 3)
        formula_ = matlabFunction(equation_orthographic);
        parameter_11 = @(x) abs(power(2 * x ./ ((1 - 1) * r_starb - (1 + 1) * r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 - 1) * r_fort - (1 + 1) * r_aft), 2) - 1);
        lb = x_position_12;
        ub = 0;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(3), feval(3)] = fminbnd(parameter_11, lb, ub, options);
end
if ~ismember(exclusion , 4)
        formula_ = matlabFunction(equation_orthographic);
        parameter_11 = @(x) abs(power(2 * x ./ ((1 + 1) * r_starb - (1 - 1) * r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 - 1) * r_fort - (1 + 1) * r_aft), 2) - 1);
        lb = 0;
        ub = x_position_11;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(4), feval(4)] = fminbnd(parameter_11, lb, ub, options);
end
non_zero_indices = find(x_min ~= 0);
x_min = x_min(non_zero_indices);
differences = abs(x_min - intersection);
[~, index] = min(differences);
reserve_x = x_min(index);
reserve_y = formula_(reserve_x);

    elseif flag == 'R_2'

if ~ismember(exclusion, 1)
        formula_ = matlabFunction(equation_orthographic);
        parameter_21 = @(x) abs(power(2 * x ./ ((1 + 1) *2* r_starb - (1 - 1) *2* r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 + 1) *2* r_fort - (1 - 1) * 2*r_aft), 2) - 1);
        lb = 0;
        ub = x_position_21;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(1), feval(1)] = fminbnd(parameter_21, lb, ub, options);
end


if ~ismember(exclusion , 2)
        formula_ = matlabFunction(equation_orthographic);
        parameter_21 = @(x) abs(power(2 * x ./ ((1 - 1) *2* r_starb - (1 + 1) *2* r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 + 1) *2* r_fort - (1 - 1) *2* r_aft), 2) - 1);
        lb = x_position_22;
        ub = 0;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(2), feval(2)] = fminbnd(parameter_21, lb, ub, options);
end
if ~ismember(exclusion , 3)
        formula_ = matlabFunction(equation_orthographic);
        parameter_21 = @(x) abs(power(2 * x ./ ((1 - 1) *2* r_starb - (1 + 1) *2* r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 - 1) *2* r_fort - (1 + 1) *2* r_aft), 2) - 1);
        lb = x_position_22;
        ub = 0;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(3), feval(3)] = fminbnd(parameter_21, lb, ub, options);
end
if ~ismember(exclusion , 4)
        formula_ = matlabFunction(equation_orthographic);
        parameter_21 = @(x) abs(power(2 * x ./ ((1 + 1) *2* r_starb - (1 - 1) *2* r_port), 2) + ...
            power(2 * formula_(x) ./ ((1 - 1) *2* r_fort - (1 + 1) *2* r_aft), 2) - 1);
        lb = 0;
        ub = x_position_21;
        options = optimset('MaxFunEvals', 100, 'MaxIter', 10);
        [x_min(4), feval(4)] = fminbnd(parameter_21, lb, ub, options);
end
non_zero_indices = find(x_min ~= 0);
x_min = x_min(non_zero_indices);
differences = abs(x_min - intersection);
[~, index] = min(differences);
reserve_x = x_min(index);
reserve_y = formula_(reserve_x);
end
end