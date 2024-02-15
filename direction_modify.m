%Author: Zhang Ziqing
%University: Wuhan University of Technology
%Contact Email: 322980@whut.edu.cn
%Feedback and Support
%If you have any questions, suggestions, or would like to contribute, feel free to contact me 
% or raise an issue on the GitHub Issues page.
function [U, u_DCPA, u_TCPA, u_D, u_C, u_K, distance_D]= direction_modify(ship_opposite_rotated_x, ship_opposite_rotated_y, ship_origin_v, v_opposite, slope_opposite, flag) 
% clc;
% clear;
% clearvars , close;
% num = 20;
% number_total = linspace(1 , 1000 , 20);
% 
% u_K_number = zeros(1,num);
% u_C_number = zeros(1,num);
% u_D_number = zeros(1,num);
% u_DCPA_number = zeros(1,num);
% u_TCPA_number = zeros(1,num);

% for p = 1 : num
% number = number_total(p);
global call_counter
call_counter = call_counter + 1;
    
disp(['Function called for the ', num2str(call_counter), ' time.']);

syms t x

L = 7;

global r_fort r_aft r_starb r_port
k_ad = 10.^(0.3591 * log10(ship_origin_v) + 0.0952);
k_dt = 10.^(0.5441 * log10(v_opposite) - 0.0795);
r_fort = (1 + 1.34 * sqrt(k_ad.^2 + k_dt.^2)) * L;
r_aft = (1 + 0.67 * sqrt(k_ad.^2 + k_dt.^2)) * L;
r_starb = (0.2 + k_ad) * L;
r_port = (0.2 + 0.75 * k_ad) * L;

global equation_orthographic equation_combine

% 分段函数定义
% segmentedFunction = @(m) (m < 0) * -1 + (m >= 0) * 1;

% 参数方程
% parameter_1 = @(x, y)  power(2 * x / ((1 + segmentedFunction(x)) * r_starb - (1 - segmentedFunction(x)) * r_port), 2) + ...
%     power(2 * y / ((1 + segmentedFunction(y)) * r_fort - (1 - segmentedFunction(y)) * r_aft), 2) - 1;
% parameter_2 = @(x, y)  power(2 * x / ((1 + segmentedFunction(x)) * 2 * r_starb - (1 - segmentedFunction(x)) * 2 * r_port), 2) + ...
%     power(2 * y / ((1 + segmentedFunction(y)) * 2 * r_fort - (1 - segmentedFunction(y)) * 2 * r_aft), 2) - 1;

% 创建figure

% % 在figure中绘制参数方程曲线
% h1 = ezplot(parameter_1 ,[-500, 500, -500, 1000]);
% hold on;
% h2 = ezplot(parameter_2 ,3 .* [-500, 500, -500, 1000]);
% hold on;
% % 设置曲线颜色和线宽
% set(h1, 'Color', 'black', 'LineWidth', 1);
% set(h2, 'Color', 'black', 'LineWidth', 1);
% 
% title('Ellipse');
% xlabel('x');
% ylabel('y');
% axis equal;  % 保持坐标轴比例一致
% axis on
% % 去除网格线
% grid off;


% trajectory_opposite = scatter_point (flag) ;
syms t x
% fplot(trajectory_opposite, [-1000, 1500],LineStyle="-",Color='r');

% 计算新坐标

angle_with_x_positive = atan2(slope_opposite, 1);
dx = cos(angle_with_x_positive) * v_opposite*flag;
if slope_opposite<=0
   dy = -ship_origin_v - abs(sin(angle_with_x_positive)) * v_opposite*flag;
end
if slope_opposite>0
    dy = -ship_origin_v + abs(sin(angle_with_x_positive)) * v_opposite*flag;
end
slope_combine = dy/dx;
if slope_opposite>0
    if flag>0
        angle_Cts = pi/2 - atan(slope_opposite);
    else
        angle_Cts = -pi/2 - atan(slope_opposite);
    end
else
    if flag<0
        angle_Cts = -(pi/2 + atan(slope_opposite));
    else
        angle_Cts = pi/2 - atan(slope_opposite);
    end
end
angle_Cts = rad2deg(angle_Cts);

equation_combine = slope_combine*(x - ship_opposite_rotated_x) + ship_opposite_rotated_y;
equation_orthographic = -1./slope_combine*x;
vertical_intercept = ship_opposite_rotated_y - slope_combine*ship_opposite_rotated_x;
horizontal_intercept = -ship_opposite_rotated_y/slope_combine+ship_opposite_rotated_x;

if (-1/slope_combine)>0
    flag_equation_orthographic = [2,4];
else
    flag_equation_orthographic = [1,3];
end

if vertical_intercept>=0&&horizontal_intercept>=0
    flag_equation_combine = 3;
elseif vertical_intercept>=0&&horizontal_intercept<0
    flag_equation_combine = 4;
elseif vertical_intercept<0&&horizontal_intercept<0
    flag_equation_combine = 1;
elseif vertical_intercept<0&&horizontal_intercept>=0
    flag_equation_combine = 2;
end

% fplot (equation_combine)
% hold on;
% fplot (equation_orthographic)
% hold on;

intersection_combine_orthographic_x = solve(equation_combine==equation_orthographic , x);

global x_position_11 x_position_12 x_position_21 x_position_22;

x_position_11 = (r_starb);
x_position_12 = (- r_port);

x_position_21 = (r_starb*2);
x_position_22 = (- r_port*2);

[intersection_D_0_x , intersection_D_0_y] = piecewise_equ_modify('D_0' , intersection_combine_orthographic_x , ship_opposite_rotated_x , flag_equation_combine, flag);
[intersection_D_1_x , intersection_D_1_y] = piecewise_equ_modify('D_1' , intersection_combine_orthographic_x , ship_opposite_rotated_x , flag_equation_combine, flag);

x_position_11 = (r_starb);
x_position_12 = (- r_port);

x_position_21 = (r_starb*2);
x_position_22 = (- r_port*2);

[intersection_R_1_x , intersection_R_1_y] = piecewise_equ_modify('R_1' , intersection_combine_orthographic_x , ship_opposite_rotated_x , flag_equation_orthographic, flag);
[intersection_R_2_x , intersection_R_2_y] = piecewise_equ_modify('R_2' , intersection_combine_orthographic_x , ship_opposite_rotated_x , flag_equation_orthographic, flag);

m = slope_combine;
b = -slope_combine * ship_opposite_rotated_x + ship_opposite_rotated_y;
distance_DCPA = norm([intersection_combine_orthographic_x , -1./slope_combine*intersection_combine_orthographic_x] , 2);
distance_DCPA = double(distance_DCPA);

distance_D = norm([ship_opposite_rotated_x , ship_opposite_rotated_y] , 2);
distance_D_0 = norm([intersection_D_0_x , intersection_D_0_y] , 2);
distance_D_1 = norm([intersection_D_1_x , intersection_D_1_y] , 2);
distance_R_1 = norm([intersection_R_1_x , intersection_R_1_y] , 2);
distance_R_2 = norm([intersection_R_2_x , intersection_R_2_y] , 2);

nq = 1;
ns = 1;
if isempty(intersection_D_1_x)
    ns = 0;
end
if isempty(intersection_D_0_x)
    nq = 0;
end


if nq > 0 && ns > 0
    u_DCPA = 1;
elseif ns == 0
    u_DCPA = 0;
else
    u_DCPA = 1/2 - 1/2 * sind(180 / (distance_R_2 - distance_R_1) * (distance_DCPA - (distance_R_2 + distance_R_1) / 2));
end
disp([nq, ns]);

TCPA = sqrt(distance_D^2 - distance_DCPA^2)/ship_origin_v;
T0 = sqrt(distance_D_1^2 - distance_DCPA^2)/ship_origin_v;
if ns > 0
    u_TCPA = exp(-log(2) * (abs(TCPA) / T0)^2);
else 
    u_TCPA = 0;
end
% u_TCPA_number(p) = subs(u_TCPA , 't' , number);


if ns > 0
    u_D = exp(-log(2) * ((distance_D - distance_D_0) / (distance_D_1 - distance_D_0))^2);
else
    u_D = 0;
end
% u_D_number(p) = subs(u_D , 't' , number);

angle_rad = atan2(m, 1); 
angle_deg = rad2deg(angle_rad);

if ship_opposite_rotated_x>0
    if m>0
        angle_C = 90 - angle_deg;
    else
        angle_C = 180-(angle_deg+90);
    end
elseif ship_opposite_rotated_x<0
    if m>0
        angle_C = -(90+angle_deg);
    else
        angle_C = -(90 + angle_deg);
    end
end


if ns > 0
    u_C = (17/44) * (cosd(angle_C + 161) + sqrt(cosd(angle_C + 161)^2 + 440/289));
else
    u_C = 0;
end
% u_C_number(p) = subs(u_C , 't' , number);


K = ship_origin_v/v_opposite;
if ns > 0
    u_K = 1 / (1 + 2 / (K * sqrt(K^2 + 1 + 2 * K * abs(sind(abs(angle_Cts))))));
else
    u_K = 0;
end
% u_K_number(p) = subs(u_K , 't' , number);

U = 0.36*u_DCPA+0.32*u_TCPA+0.14*u_D+0.1*u_C+0.08*u_K;
% % U = subs(U , 't' , number);
% U_number(p) = U;

% [ optimised_parameters ] = Particle_Swarm_Optimization (Bird_in_swarm, Number_of_quality_in_Bird, MinMaxRange, U, availability_type, 2, 2, 2, 0.4, 0.9, max_iteration);
 
end
% end
% figure;
% plot(number_total, u_K_number, 'LineWidth', 1, 'Color', rand(1,3));
% hold on
% plot(number_total, u_C_number, 'LineWidth', 1, 'Color', rand(1,3));
% hold on
% plot(number_total, u_D_number, 'LineWidth', 1, 'Color', rand(1,3));
% hold on
% plot(number_total, u_DCPA_number, 'LineWidth', 1, 'Color', rand(1,3));
% hold on
% plot(number_total, u_TCPA_number, 'LineWidth', 1, 'Color', rand(1,3));

% Adding legend
% legend('u\_K\_number', 'u\_C\_number', 'u\_D\_number', 'u\_DCPA\_number', 'u\_TCPA\_number');
% 
% % Other settings if needed
% xlabel('X-axis Label');
% ylabel('Y-axis Label');
% title('Your Plot Title');
% grid on;

