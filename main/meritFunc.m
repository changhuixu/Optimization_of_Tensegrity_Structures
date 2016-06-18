function [fval] = meritFunc(x)
% function [fval] = meritFunc(x)

global mu;

fval = energy(x) + mu * norm(constraintE(x));

end