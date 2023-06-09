% Define the ranges for the design variables
x1_range = linspace(0.0625, 99.99, 100);  % Range for x1
x2_range = linspace(0.0625, 99.99, 100);  % Range for x2
x3_range = linspace(10, 200, 100);        % Range for x3
x4_range = linspace(10, 200, 100);        % Range for x4

% Create a grid of design variable values
[x3_grid, x4_grid] = meshgrid(x3_range, x4_range);

% Initialize the cost values matrix
cost_values = zeros(size(x3_grid));

% Calculate the cost values for each combination of x3 and x4
for i = 1:size(x3_grid, 1)
    for j = 1:size(x3_grid, 2)
        % Set the current values of x3 and x4
        x3 = x3_grid(i, j);
        x4 = x4_grid(i, j);

        % Create the design vector with x3 and x4 and fill the remaining x1 and x2 with their initial values
        x = [1; 1; x3; x4];

        % Calculate the cost value for the current design vector
        cost_values(i, j) = objective_function(x);
    end
end

% Create the contour plot
contourf(x3_grid, x4_grid, cost_values, 20, 'LineWidth', 0.5);
colorbar;
xlabel('x3');
ylabel('x4');
title('Cost Landscape: Pressure Vessel Optimization');

% Display the plot
figure;

