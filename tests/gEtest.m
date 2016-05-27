function [gEtest] = gEtest(x)
% function [gEtest] = gEtest(x)

% Perform gradient calculations for the Energy function numerically
% Input:
%	x:       xyz coordinates of ending points of tubes/cables
%           (3*5)
% Output:
%	gEtest:  the gradient of the system total energy


gEvec = gradientE(x);
%gEvec = gradientFixedPoints(gEvec);

% initialize a small matrix d with elements 1e-6
d = randn(size(x))*1e-6;
dvec = d(:);

gEtest = abs(energy(x+d)-energy(x-d)-2*(gEvec.')*dvec)/norm(d);
end