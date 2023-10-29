function [X] = grid_points(division_number)

% compute grid points
x = linspace(0, 1, 2 * division_number + 1);
x = x(2:2:2 * division_number + 1);
[grid_points_x, grid_points_y] = meshgrid(x, x);
X = [reshape(grid_points_x, [], 1), reshape(grid_points_y, [], 1)]';
X = reshape(X, [], 1);

end