function [gLagrange] = gradientLagrange(x,lambda)
% function [gLagrange] = gradientLagrange(x,lambda)

% Perform gradient calculations for the Lagrange function 
% Input:
%	x:          the coordinates matrix of the system
%               (3*5)
%   lambda:     the lambda vector in Lagrange function
%               (4*1)
% Output:
%	gLagrange:	the gradient (wrt x) of Lagrange function
%               (15*1)

gradFx = gradientE(x);
gradCx = gradientC(x);

gLagrange = gradFx - gradCx.'*lambda;

end