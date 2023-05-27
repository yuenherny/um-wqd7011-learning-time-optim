% benchmark
clc; clear;

pkg load symbolic % install first
source("solvers.m")
source("hessian_approx.m")

% initialize symbol
syms x y;
% objection function
f = x^2;
% constraints
##constraint_1 = 214*x + 616*y; % <= 4
##constraint_2 = 589*x + 209*y; % <= 5.5
##constraint_3 = 55*x + 33*y; % <= 7
##
##constraints = [constraint_1; constraint_2; constraint_3];

% derive gradient
gradient_mat = gradient(f);
df = function_handle(gradient_mat);

% optim param
f = function_handle(f);
iterations = 1000;
alpha = 0.01;
c_init = rand(1) * 10000;
c_best = c_init;
f_best = f(c_init);
c_minima = 0;

% pass gradient into solver function and return next point
c_curr = c_init;
disp("Starting iterations")
for iter = 1:iterations

  c_next = steepest_descent_univariate(c_curr, df, alpha);
  f_next = f(c_next);

  if f_next < f_best
    c_best = c_next;
    f_best = f_next
  endif

  if abs(f_best) <= 1e-5
    break
  endif

  c_curr = c_next;

endfor

disp(["Best x: ", num2str(c_best)])
disp(["Best f: ", num2str(f_best)])
