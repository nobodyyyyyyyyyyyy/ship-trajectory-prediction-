%Author: Zhang Ziqing
%University: Wuhan University of Technology
%Contact Email: 322980@whut.edu.cn
%Feedback and Support
%If you have any questions, suggestions, or would like to contribute, feel free to contact me 
% or raise an issue on the GitHub Issues page.
function [ship_opposite_rotated_1_x , ship_opposite_rotated_1_y , ship_opposite_1_gradient] = transition_matrix(ship_opposite_1_x_present, ship_opposite_1_y_present, ship_origin_x, ship_origin_y, theta)
ship_opposite_translate_1_x = [ship_opposite_1_x_present(end-1), ship_opposite_1_x_present(end)] - ship_origin_x;
ship_opposite_translate_1_y = [ship_opposite_1_y_present(end-1), ship_opposite_1_y_present(end)] - ship_origin_y;
ship_opposite_rotated_1_x = ship_opposite_translate_1_x .* [cosd(theta(end-1)), cosd(theta(end))] - ship_opposite_translate_1_y .* [sind(theta(end-1)) ,sind(theta(end))] ;
ship_opposite_rotated_1_y = ship_opposite_translate_1_x .* [sind(theta(end-1)) ,sind(theta(end))] + ship_opposite_translate_1_y .* [cosd(theta(end-1)), cosd(theta(end))];

ship_opposite_1_gradient = diff(ship_opposite_rotated_1_y) ./ diff(ship_opposite_rotated_1_x);
ship_opposite_rotated_1_x = ship_opposite_rotated_1_x(end);
ship_opposite_rotated_1_y = ship_opposite_rotated_1_y(end);
end