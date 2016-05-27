function [hEtest] = hEtest(x)
% function [hEtest] = hEtest(x)

% Perform hessian calculations for the Energy function numerically
% Input:
%	x:       xyz coordinates of ending points of tubes/cables
%           (3*5)
% Output:
%	hEtest:  the hessian of the system total energy


% initialize a small matrix d with elements 1e-6
d = randn(size(x))*1e-6;
dvec = d(:);

gradE1 = gradientE(x+d);
gradE2 = gradientE(x-d);

hessE = hessianE(x);

temp = gradE1-gradE2-2*(hessE.')*dvec;
hEtest = norm(temp)/norm(d);
end