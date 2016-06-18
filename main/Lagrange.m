function [Lval] = Lagrange(x,lambda)
% function [Lval] = Lagrange(x,lambda)

% compute the value of Lagrange function

Lval= energy(x) - constraintE(x).'*lambda;

end