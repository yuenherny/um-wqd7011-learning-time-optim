% benchmark
clc; clear;

pkg load symbolic % install first
source("../src/solvers.m")
source("../src/hessian_approx.m")

% initialize symbol
syms x y;
% objection function
f = 100*(y - x^2)^2 + (1-x)^2

% derive gradient
gradient_mat = gradient(f);
df = function_handle(gradient_mat);

% optim param
f = function_handle(f);
iterations = 1000;
alpha = 0.0001;
c_init = [rand(1) * 5; rand(1) * 5]
c_best = c_init;
f_best = f(c_init(1), c_init(2))

% pass gradient into solver function and return next point
c_curr = c_init;
disp("Starting iterations...\n")
for iter = 1:iterations

  c_next = steepest_descent_bivariate(c_curr, df, alpha);
  f_next = f(c_next(1), c_next(2));

  if f_next < f_best
    c_best = c_next;
    f_best = f_next;
  endif

  if abs(f_best) <= 1e-5
    break
  endif

  c_curr = c_next;

endfor
disp("End of iterations...\n")

% Best solution is c_best=[1,...,1] and f_best=0
disp("Best x:")
disp(c_best)
printf("Best f: %f\n", f_best)
