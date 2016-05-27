function [hLagrange] = hessianLagrange(x,lambda)
% function [hLagrange] = hessianLagrange(x,lambda)

% Perform hessian calculations for the Lagrange function 
% Input:
%	x:          the coordinates matrix of the system
%               (3*5)
%   lambda:     the lambda vector in Lagrange function
%               (4*1)
% Output:
%	hLagrange:	the hessian (wrt x) of Lagrange function
%               (15*15)
global tube_pts;

number_of_tubes = length(tube_pts);
number_of_points = length(x);

hessFx = hessianE(x);
hessCx = hessianC(x);
hLagrange = zeros(number_of_points*3,number_of_points*3);

for i = 1:number_of_tubes
    hLagrange = hLagrange + lambda(i)*hessCx(:,:,i);
end

hLagrange = hessFx - hLagrange;

end