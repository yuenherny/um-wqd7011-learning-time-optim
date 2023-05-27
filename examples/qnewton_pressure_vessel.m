% benchmark
clc; clear;

pkg load symbolic % install first
source("../src/solvers.m")
source("../src/hessian_approx.m")

% initialize symbol
syms w x y z;
% objection function
f = 0.6224*w*y*z + 1.7781*x*y^2 + 3.1161*w^2*z + 19.84*w^2*y

% derive gradient
gradient_mat = gradient(f);
df = function_handle(gradient_mat);

% derive hessian
hessian_mat = hessian(f);
ddf = function_handle(hessian_mat);

% optim param
f = function_handle(f);
iterations = 10;
alpha = 0.0001;
c_init = rand([1 4]) * 200
c_best = c_init;
f_best = f(c_init(1), c_init(2), c_init(3), c_init(4));
B_init = eye(4);

% pass gradient into solver function and return next point
c_curr = c_init;
B_curr = B_init;
disp("Starting iterations...\n")
for iter = 1:iterations

  c_next = quasi_newton_sr1_quadravariate(c_curr, df, alpha, B_curr);
  f_next = f(c_next(1), c_next(2), c_next(3), c_next(4));

  if f_next < f_best
    c_best = c_next;
    f_best = f_next;
   endif

  if abs(f_best) <= 1e-5
    disp("Minimal error condition met. Iteration halted.")
    break
  endif

  B_next = sr1_quadravariate(B_curr, c_next, c_curr, df);
  c_curr = c_next;

endfor

disp("End of iterations...\n")

% Best solution is c_best=[1,...,1] and f_best=0
% Require iter=1e5 to meet minimal error condition
disp("Best x:")
disp(c_best)
printf("Best f: %f\n", f_best)
