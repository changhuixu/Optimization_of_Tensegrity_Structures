clc;
clear;
% start stop-watch timer
tic

% load constants
data_sun_river;
global fixed_x;
% split initial coordinates into variables and fixed points
[x,fixed_x] = FixPoints(x);
% join x variables and fixed points. for the purpose of plotting and so on
JoinFixedPoints(x, fixed_x);

% calculate the system energy
Initial_Energy = energy(x)
% calculate the gradient of energy
gF = gradientE(x);
test_gE_val = gEtest(x)

% calculate the hessian of energy
hessF = hessianE(x);
test_hE_val = hEtest(x)

% calculate the equality constraints.
constraintE(x);
% calculate the gradient of equality constraints
gradientC(x);
test_gC_val = gCtest(x)

% calculate the hessian of equality constraints
hessianC(x);
test_hC_val = hCtest(x)

% compute lambda
[lambda, gLagrange]=findLambda(x);
% initial lambda vector is [8.6019;2.9503;2.9503;6.7804]
norm_lambda = norm(lambda)

% compute Lagrange Function
Lagrange(x,lambda);  % result is 9.7703
% compute the gradient of Lagrange Function
gLagrange1 = gradientLagrange(x,lambda);
test_gL_val = gLtest(x,lambda)
gL_diff = gLagrange - gLagrange1; % result is a zero vector

% compute the hessian of Lagrange Function
hessianLagrange(x,lambda);
test_hL_val = hLtest(x,lambda)

% stop timer
timer = toc;
fprintf('Total running time is %d seconds.\n',timer);

