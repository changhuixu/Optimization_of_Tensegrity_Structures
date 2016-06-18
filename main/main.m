%%%%%% function main

clc;
clear;
% start stop-watch timer
tic;

%initialize the system with random coordinates (for test purpose)
x = rand(3,8);

% load constants provided by professor.
data_sun_river;
PlotTensegrity(x);
global fixed_x;
% split initial coordinates into variables and fixed points
[x,fixed_x] = FixPoints(x);

% set parameters for Newton Trust Region method
Delta_max	= 1;
Delta       = 0.05;
eta         = 0.15;
tol         = 1e-5;
trace       = 1;


% do the optimization
 [x,fval,xs,Ds] = ...
    newton_tr(x,Delta,Delta_max,eta,tol,trace);

% join x variables and fixed points. for the purpose of plotting and so on

x = JoinFixedPoints(x, fixed_x);

fprintf('main: optimized system energy = %g\n', fval);

% plot the final configuration of the system (x)
PlotTensegrity(x);

% stop timer
timer = toc;
fprintf('\nTotal running time is %d seconds.\n',timer);

%%%%%% end