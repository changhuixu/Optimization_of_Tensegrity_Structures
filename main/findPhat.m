function [phat,Q2] = findPhat(x)
% function [phat,Q2] = findPhat(x)

% approximately solve min_phat norm(phat)
%                       subject to B^T * phat = r
%
% Input:
%	x:      the coordinates matrix of the system
%               (3*5)
% Output:
%	phat:	the phat value in model function
%               (15*1)
%   Q2:     the orthonormal basis for the kernel of B^T
%               (15*11)

B = gradientC(x).';   % (4*15)' matrix
r = -constraintE(x);    % (4*1) vector

% Assumes columns of B are linearly independent
% Uses QR factorization
[m,n] = size(B);        % m=15, n=4
[Q,R] = qr(B);
Q1 = Q(:,1:n);
Q2 = Q(:,n+1:m);
phat = Q1*(R(1:n,1:n)' \ r);

% phat_opt = -r+B'*phat

end