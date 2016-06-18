function [lambda, gLagrange] = findLambda(x)
% function [lambda, gLagrange] = findLambda(x)

% Perform QR factorization to find lambda 
% Input:
%	x:          the coordinates matrix of the system
%               (3*5)
% Output:
%	lambda:     the lambda in Lagrange function
%               (4*1)
%   gLagrange:  the minimization function value
%               gE - B*lambda
%               (15*1)

gC= gradientC(x);
B = gC.';
gE= gradientE(x);
r = gE;
% Assumes columns of B are linearly independent
% Uses QR factorization
[~,n] = size(B);
[Q,R] = qr(B);
Q1 = Q(:,1:n);
lambda = R(1:n,1:n) \ Q1'*r;
gLagrange = gE - B*lambda;

end