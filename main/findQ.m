function [q,at_bdry] = findQ(phat,N,hessL,gradientF,Delta,trace)
% function [q,at_bdry] = findQ(x,phat,N,Delta,tol,trace)

% approximately solve min_phat m(phat+N*q)
%                       subject to ||q|| <= sqrt(delta^2 - norm(phat)^2)
% 
% Input:
%	x:      the coordinates matrix of the system
%               (3*5)
%   phat:   size (15*1)
%   N:      size (15*11)
%
% Output:
%	q:      the q in model function m(phat+N*q)
%               (11*1)

B = N.'*hessL*N;
g = N.'*gradientF + N.'*hessL*phat;

% delta for q, calculate when norm(phat)<=Delta
Delta1 = sqrt(Delta^2 - norm(phat)^2);

%[q,at_bdry] = solve_tr(g,B,Delta1,1e-2,trace);
[q,at_bdry] = dogleg_tr(g,B,Delta1,trace);

end