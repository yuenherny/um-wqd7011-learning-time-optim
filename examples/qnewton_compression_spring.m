% benchmark
clc; clear;

pkg load symbolic % install first
pkg load statistics
source("../src/solvers.m")
source("../src/hessian_approx.m")
source("../src/constraints.m")

% initialize symbol
syms x y z;
% objection function
f = (z + 2)*y*x^2
f_sol = 0.012665;
c_sol = [0.051749, 0.358179, 11.203763];

% derive gradient
gradient_mat = gradient(f);
df = function_handle(gradient_mat)
f = function_handle(f);

% optim param
iterations = 5000;
alpha = 0.0001;

% Simple bounds of the search domain
% Lower bounds and upper bounds
Lb = [0.05; 0.25; 2];
Ub = [2; 1.3; 15];

% Random initial solutions
for i=1:length(Lb),
  c_init(i)=Lb(i)+(Ub(i)-Lb(i))*rand(1);
end

c_curr = transpose(c_init)
c_best = c_curr;
f_curr = f(c_curr(1), c_curr(2), c_curr(3));
f_best = f_curr;
f_next_array = zeros(iterations,1);
f_next_array(1) = f_curr;

B_curr = eye(3);
disp("Starting iterations...\n")
for iter = 1:iterations

  c_next = quasi_newton_sr1_trivariate(c_curr, df, alpha, B_curr);
  c_next = simplebounds(c_next, Lb, Ub)
  f_next = f(c_next(1), c_next(2), c_next(3));
  f_next_array(iter+1) = f_next;
  f_next = f_next + apply_compression_spring_constraints(c_next)

  if f_next < f_best
    c_best = c_next;
    f_best = f_next;
   endif

  if abs(f_best - f_sol) <= 1e-5
    disp("Minimal error condition met. Iteration halted.")
    break
  endif

  % for next iteration
  B_next = sr1_trivariate(B_curr, c_next, c_curr, df);
  c_curr = c_next;
  f_curr = f_next;

endfor

printf("Ending after of %i iterations...\n", iter)

% Best solution is c_best=[1,...,1] and f_best=0
% Require iter=1e5 to meet minimal error condition
disp("Best x:")
disp(c_best)
printf("Best f: %f\n", f_best)
