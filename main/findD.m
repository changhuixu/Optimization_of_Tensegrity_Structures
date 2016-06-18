function [d] = findD(x)
% function [d] = findD(x)

% approximately solve min_phat norm(phat)
%                       subject to B' * d = r
%
% Input:
%	x:      the coordinates matrix of the system
%               (3*5)
% Output:
%	d:	the phat value in model function
%               (15*1)


B = gradientC(x).';   % (4*15) matrix
r = -constraintE(x);    % (4*1) vector

% Assumes columns of B are linearly independent
% Uses QR factorization
[~,n] = size(B);
[Q,R] = qr(B);
Q1 = Q(:,1:n);

d = Q1*(R(1:n,1:n)' \ r);

end