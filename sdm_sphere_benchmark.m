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
iterations = 100;
alpha = 0.01;
c_init = rand(1) * 10000;
c_best = c_init;
f_best = f(c_init);

% pass gradient into solver function and return next point
c_curr = c_init;
for iter = 1:iterations

  c_next = steepest_descent(c_curr, df, alpha);
  f_next = f(c_next(1), c_next(2));

  if f_next < f_best
    c_best = c_next;
    f_best = f_next;
  endif

  c_curr = c_next;

endfor

disp(c_best)
disp(f_best)
