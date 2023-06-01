% benchmark
clc; clear;

pkg load symbolic % install first
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
iterations = 1e4;
alpha = 0.01;

% Simple bounds of the search domain
% Lower bounds and upper bounds
Lb = [0.1; 0.1; 10; 10];
Ub = [99; 99; 200; 200];

% Random initial solutions
for i=1:length(Lb),
  c_init(i)=Lb(i)+(Ub(i)-Lb(i))*rand(1);
end

c_curr = transpose(c_init);
c_best = c_curr;
f_best = f(c_curr(1), c_curr(2), c_curr(3), c_curr(4));

disp("Starting iterations...\n")
for iter = 1:iterations

  c_next = steepest_descent_quadravariate(c_curr, df, alpha);
  f_next = f(c_next(1), c_next(2), c_next(3), c_next(4));
  f_next = f_next + apply_constraints(c_next)

  if f_next < f_best
    c_best = c_next;
    f_best = f_next;
   endif

  if abs(f_best - f_sol) <= 1e-5
    disp("Minimal error condition met. Iteration halted.")
    break
  endif

  c_next = simplebounds(c_next, Lb, Ub);
  c_curr = c_next;

endfor

printf("Ending after of %i iterations...\n", iter)

% Best solution is c_best=[1,...,1] and f_best=0
% Require iter=1e5 to meet minimal error condition
disp("Best x:")
disp(c_best)
printf("Best f: %f\n", f_best)
