% main
clc; clear;

pkg load symbolic % install first
source("solvers.m")

% initialize symbol
syms x y;
% objection function
cost_func = 3*x + 5*y;
% constraints
constraint_1 = 214*x + 616*y; % <= 4
constraint_2 = 589*x + 209*y; % <= 5.5
constraint_3 = 55*x + 33*y; % <= 7

constraints = [constraint_1; constraint_2; constraint_3];

% pass cost function and constraints into solver function and return optimal point

