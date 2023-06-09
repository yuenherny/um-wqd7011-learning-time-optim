% benchmark
clc; clear;

pkg load symbolic % install first
pkg load statistics
source("../src/solvers.m")
source("../src/hessian_approx.m")
source("../src/constraints.m")

% initialize symbol
syms w x y z;
% objection function
f = 0.6224*w*y*z + 1.7781*x*y^2 + 3.1161*w^2*z + 19.84*w^2*y
f_sol = 6059.714339;
c_sol = [0.8125; 0.4375; 42.098446; 176.636596];

% derive gradient
gradient_mat = gradient(f);
df = function_handle(gradient_mat)
f = function_handle(f);

% optim param
iterations = 1000;
alpha = 0.01;

% Simple bounds of the search domain
% Lower bounds and upper bounds
Lb = [0.1; 0.1; 10; 100];
Ub = [1; 1; 100; 200];

% Random initial solutions
##for i=1:length(Lb),
##  c_init(i)=Lb(i)+(Ub(i)-Lb(i))*rand(1);
##end
c_init = [1 1 100 200]

c_curr = transpose(c_init)
for i = 1:2
  t = c_curr(i);
  c_curr(i) = compute_metal_thickness(t);
endfor
c_best = c_curr;
f_curr = f(c_curr(1), c_curr(2), c_curr(3), c_curr(4));
f_curr = f_curr + apply_pressure_vessel_constraints(c_curr);

f_best = f_curr;
f_best_array = zeros(iterations,1);
f_best_array(1) = f_best;

B_curr = eye(4);
disp("Starting iterations...\n")
for iter = 1:iterations

  c_next = quasi_newton_quadravariate(c_curr, df, alpha, B_curr);
  c_next = simplebounds(c_next, Lb, Ub);
  for i = 1:2
    t = c_curr(i);
    c_curr(i) = compute_metal_thickness(t);
  endfor
  f_next = f(c_next(1), c_next(2), c_next(3), c_next(4));
##  f_best_array(iter+1) = f_next;
  f_next = f_next + apply_pressure_vessel_constraints(c_next)
  f_best_array(iter+1) = f_next;

  if f_next < f_best
    c_best = c_next;
    f_best = f_next;
   endif

  if abs(f_best - f_sol) <= 1e-5
    disp("Minimal error condition met. Iteration halted.")
    break
  endif

  % for next iteration
  B_next = bfgs_quadravariate(B_curr, c_next, c_curr, df);
  c_curr = c_next;

endfor

printf("Ending after of %i iterations...\n", iter)

% Best solution is c_best=[1,...,1] and f_best=0
% Require iter=1e5 to meet minimal error condition
disp("Best x:")
disp(c_best)
printf("Best f: %f\n", f_best)
