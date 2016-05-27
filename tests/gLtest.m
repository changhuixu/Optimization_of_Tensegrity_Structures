function [gLtest] = gLtest(x,lambda)
% function [gLtest] = gLtest(x,lambda)

% Perform gradient test for the Lagrange function numerically
% Input:
%	x:          xyz coordinates of ending points of tubes/cables
%               (3*5)
%   lambda:     the lambda in Lagrange function
%               (4*1)
% Output:
%	gLtest:     the gradient of the system total energy


gLvec = gradientLagrange(x,lambda);

% initialize a small matrix d with elements 1e-6
d = randn(size(x))*1e-6;
dvec = d(:);

gLtest = abs(Lagrange(x+d,lambda)-Lagrange(x-d,lambda)-2*(gLvec.')*dvec)/norm(d);

end