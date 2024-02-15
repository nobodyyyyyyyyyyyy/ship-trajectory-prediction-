clc,clear;
global call_counter
call_counter = 0;
for n = 1:50
    basic_velocity = 5 + randn() * 1.5;
    basic_time = 30;
    

    a = basic_velocity;
    tmax = basic_time;
    t = linspace(19.5, tmax + 3, 20);

    ship_origin_y = a * t;
    ship_origin_x = 0;
    theta = 0;
    flag = 1;

    % 船对1的数据
    ship_opposite_1_x = -a * tmax + randn() * (-a * tmax * 0.1);
    ship_opposite_1_y = a * tmax * 2 + randn() * 2 * a * tmax * 0.1;
    ship_opposite_1_angel = -135 + 5 * (2 * rand - 1);
    ship_opposite_1_y_present = ship_opposite_1_y + a * sqrt(2) * t * cosd(ship_opposite_1_angel);
    ship_opposite_1_x_present = ship_opposite_1_x - a * sqrt(2) * t * sind(ship_opposite_1_angel);
    ship_opposite_1_gradient = diff(ship_opposite_1_y_present) ./ diff(ship_opposite_1_x_present);

    for i = 1:20
    ship_opposite_translate_1_x = ship_opposite_1_x_present - ship_origin_x;
    ship_opposite_translate_1_y = ship_opposite_1_y_present - ship_origin_y(i);
    ship_opposite_rotated_1_x = ship_opposite_translate_1_x * cosd(theta) - ship_opposite_translate_1_y * sind(theta);
    ship_opposite_rotated_1_y = ship_opposite_translate_1_x * sind(theta) + ship_opposite_translate_1_y * cosd(theta);

    flag_matrix = diff(ship_opposite_rotated_1_x);
    ship_opposite_1_gradient = diff(ship_opposite_rotated_1_y) ./ flag_matrix;
    flag_matrix = [flag_matrix(1), circshift(flag_matrix, [0, 1])];
    ship_opposite_1_gradient = [ship_opposite_1_gradient(1), circshift(ship_opposite_1_gradient, [0, 1])];
    if flag_matrix(i) <0
        flag = -1;
    end
        [dataStruct(n).U(i), dataStruct(n).u_DCPA(i), dataStruct(n).u_TCPA(i), dataStruct(n).u_D(i), dataStruct(n).u_C(i), dataStruct(n).u_K(i),dataStruct(n).D(i)] = direction_modify(ship_opposite_rotated_1_x(i), ship_opposite_rotated_1_y(i), a, a * sqrt(2), ship_opposite_1_gradient(i),flag);
    end

    % 船对2的数据
    ship_opposite_2_x = -a * tmax + randn() * -a * tmax * 0.1;
    ship_opposite_2_y = a * tmax + randn() * a * tmax * 0.1;
    ship_opposite_2_angel = 90 + 5 * (2 * rand - 1);
    ship_opposite_2_y_present = ship_opposite_2_y + a * t * cosd(ship_opposite_2_angel);
    ship_opposite_2_x_present = ship_opposite_2_x + a * t * sind(ship_opposite_2_angel);

    for i = 1:20
    ship_opposite_translate_2_x = ship_opposite_2_x_present - ship_origin_x;
    ship_opposite_translate_2_y = ship_opposite_2_y_present - ship_origin_y(i);
    ship_opposite_rotated_2_x = ship_opposite_translate_2_x * cosd(theta) - ship_opposite_translate_2_y * sind(theta);
    ship_opposite_rotated_2_y = ship_opposite_translate_2_x * sind(theta) + ship_opposite_translate_2_y * cosd(theta);

    flag_matrix = diff(ship_opposite_rotated_2_x);
    ship_opposite_2_gradient = diff(ship_opposite_rotated_2_y) ./ flag_matrix;
    flag_matrix = [flag_matrix(1), circshift(flag_matrix, [0, 1])];
    ship_opposite_2_gradient = [ship_opposite_2_gradient(1), circshift(ship_opposite_2_gradient, [0, 1])];
    if flag_matrix(i) <0
        flag = -1;
    end
        [dataStruct2(n).U(i), dataStruct2(n).u_DCPA(i), dataStruct2(n).u_TCPA(i), dataStruct2(n).u_D(i), dataStruct2(n).u_C(i), dataStruct2(n).u_K(i), dataStruct2(n).D(i)] = direction_modify(ship_opposite_rotated_2_x(i), ship_opposite_rotated_2_y(i), a, a , ship_opposite_2_gradient(i),flag);
    end

    % 船对3的数据
    ship_opposite_3_x = -a * tmax * 2 + randn() * -a * tmax * 2 * 0.1;
    ship_opposite_3_y = -a * tmax + randn() * -a * tmax * 0.1;
    ship_opposite_3_angel = 45 + 5 * (2 * rand - 1);
    ship_opposite_3_y_present = ship_opposite_3_y + a * sqrt(2) * 2 * t * cosd(ship_opposite_3_angel);
    ship_opposite_3_x_present = ship_opposite_3_x + a * sqrt(2) * 2 * t * sind(ship_opposite_3_angel);

   for i = 1:20
    ship_opposite_translate_3_x = ship_opposite_3_x_present - ship_origin_x;
    ship_opposite_translate_3_y = ship_opposite_3_y_present - ship_origin_y(i);
    ship_opposite_rotated_3_x = ship_opposite_translate_3_x * cosd(theta) - ship_opposite_translate_3_y * sind(theta);
    ship_opposite_rotated_3_y = ship_opposite_translate_3_x * sind(theta) + ship_opposite_translate_3_y * cosd(theta);

    flag_matrix = diff(ship_opposite_rotated_3_x);
    ship_opposite_3_gradient = diff(ship_opposite_rotated_3_y) ./ flag_matrix;
    flag_matrix = [flag_matrix(1), circshift(flag_matrix, [0, 1])];
    ship_opposite_3_gradient = [ship_opposite_3_gradient(1), circshift(ship_opposite_3_gradient, [0, 1])];
    if flag_matrix(i) <0
        flag = -1;
    end
        [dataStruct3(n).U(i), dataStruct3(n).u_DCPA(i), dataStruct3(n).u_TCPA(i), dataStruct3(n).u_D(i), dataStruct3(n).u_C(i), dataStruct3(n).u_K(i),dataStruct3(n).D(i)] = direction_modify(ship_opposite_rotated_3_x(i), ship_opposite_rotated_3_y(i), a, a * sqrt(2), ship_opposite_3_gradient(i),flag);
    end
end
table1 = struct2table(dataStruct);
table2 = struct2table(dataStruct2);
table3 = struct2table(dataStruct3);
file = 'combined_dataStructs.xlsx';
writetable(table1, file, 'Sheet', 'DataStruct1');
writetable(table2, file, 'Sheet', 'DataStruct2');
writetable(table3, file, 'Sheet', 'DataStruct3');
disp('Data has been written to combined_dataStructs.xlsx with different sheets.');

figure;

% % 绘制所有船只的散点
% scatter(ship_opposite_rotated_1_x, ship_opposite_translate_1_y, 'o', 'DisplayName', 'Ship Opposite 1');
% hold on;
% scatter(ship_opposite_2_x_present, ship_opposite_2_y_present, 'o', 'DisplayName', 'Ship Opposite 2');
% scatter(ship_opposite_3_x_present, ship_opposite_3_y_present, 'o', 'DisplayName', 'Ship Opposite 3');
% scatter(zeros(1, 20), ship_origin_y, 'o', 'DisplayName', 'Ship Origin y=30');
% 
% xlabel('X');
% ylabel('Y');
% title('Scatter Plot of Ships');
% legend('Location', 'best'); % 选择最佳位置显示图例
% 
% grid on;
% hold off;

