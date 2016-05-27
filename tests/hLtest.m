function [hLtest] = hLtest(x,lambda)
% function [hLtest] = hLtest(x,lambda)

% Perform hessian test for the Lagrange function numerically
% Input:
%	x:          xyz coordinates of ending points of tubes/cables
%               (3*5)
%   lambda:     the lambda in Lagrange function
%               (4*1)
% Output:
%	hLtest:     the gradient of the system total energy


hLvec = hessianLagrange(x,lambda);

% initialize a small matrix d with elements 1e-6
d = randn(size(x))*1e-6;
dvec = d(:);

gLvec1 = gradientLagrange(x+d,lambda);
gLvec2 = gradientLagrange(x-d,lambda);

temp = gLvec1 - gLvec2 - 2*hLvec*dvec;
hLtest = norm(temp)/norm(d);

end