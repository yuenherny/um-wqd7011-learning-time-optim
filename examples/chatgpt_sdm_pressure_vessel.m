1;
pkg load optim;
function [x_opt, f_opt] = steepest_descent()
    % Define the lower bounds for the design variables
    lb = [0.0625; 0.0625; 10; 10];

    % Define the upper bounds for the design variables
    ub = [99.99; 99.99; 200; 200];

    % Initialize an empty matrix to store the initial design vectors
    x_init = [];

    % Generate 10 random initial design vectors within the bounds
    for i = 1:10
        x_rand = lb + rand(size(lb)) .* (ub - lb);
        x_init = [x_init x_rand];
    end

    % Perform steepest descent for each initial design vector
    for i = 1:size(x_init, 2)
        x = x_init(:, i);

        % Set the convergence criteria
        epsilon = 1e-6;
        max_iterations = 1000;

        % Set the step size (learning rate)
        alpha = 0.01;

        % Perform steepest descent
        for iter = 1:max_iterations
            % Calculate the gradient
            g = calculate_gradient(x);

            % Update the design vector
            x_new = x - alpha * g;

            % Apply constraints
            x_new = apply_constraints(x_new, lb, ub);

            % Calculate the objective function value for the new design
            f_new = objective_function(x_new);

            % Check for convergence
            if norm(x_new - x) < epsilon || abs(f_new - objective_function(x)) < epsilon
                break;
            end

            % Update the design vector and objective function value
            x = x_new;
            f = f_new;
        end

        % Store the optimal design vector and objective function value for each initial vector
        x_opt(:, i) = x;
        f_opt(i) = f;
    end
end

function g = calculate_gradient(x)
    % Calculate the partial derivatives of the objective function
    df_dx1 = -0.0625 / (pi * x(1) * x(2)^3) - 1 / (pi * x(2)^2) * (x(2) - x(1));
    df_dx2 = -0.0625 / (pi * x(1)^3 * x(2)) - 1 / (pi * x(1)^2) * (x(1) - x(2));
    df_dx3 = -0.0625 / (pi * x(3) * x(4)^3) - 1 / (pi * x(4)^2) * (x(4) - x(3));
    df_dx4 = -0.0625 / (pi * x(3)^3 * x(4)) - 1 / (pi * x(3)^2) * (x(3) - x(4));

    % Construct the gradient vector
    g = [df_dx1; df_dx2; df_dx3; df_dx4];
end

function f = objective_function(x)
    % Calculate the objective function value
    f = 0.6224 * x(1) * x(3) * x(4) + 1.7781 * x(2) * x(3)^2 + 3.1661 * x(1)^2 * x(4) + 19.84 * x(1)^2 * x(3);
end

function c = calculate_constraints(x)
    % Calculate the constraint values
    c1 = x(1) - 0.0193 * x(3);
    c2 = x(2) - 0.00954 * x(3);
    c3 = -pi * x(3)^2 * x(4) - (4/3) * pi * x(3)^3 + 1296000;
    c4 = (x(4) / (x(3) * sqrt(2))) - 1;

    % Combine the constraints into a column vector
    c = [c1; c2; c3; c4];
end

function x_constrained = apply_constraints(x, lb, ub)
    % Apply constraints to the design vector
    c = calculate_constraints(x);

    % Check if any constraints are violated
    if any(c > 0)
        % Project the design vector onto the feasible region using a projection method
        options = optimoptions('fmincon', 'Display', 'off');
        x_constrained = fmincon(@(x) norm(x - x_constrained), x, [], [], [], [], lb, ub, @(x) constraint_function(x), options);
    else
        x_constrained = x;
    end
end

function [c, ceq] = constraint_function(x)
    % Calculate the inequality and equality constraints
    c = [x(1) - 0.0193 * x(3);
         x(2) - 0.00954 * x(3);
         -pi * x(3)^2 * x(4) - (4/3) * pi * x(3)^3 + 1296000;
         (x(4) / (x(3) * sqrt(2))) - 1];
    ceq = [];
end

% Call the steepest descent function to obtain the optimal solution
[x_opt, f_opt] = steepest_descent();

% Display the results
disp('Optimal Design Vector:');
disp(x_opt);
disp('Optimal Objective Function Value:');
disp(f_opt);

