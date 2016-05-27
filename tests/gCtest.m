function [conGtest] = gCtest(x)
% function [conGtest] = gCtest(x)

% Perform tests for gradient of constraint functions
% Input:
%	x:  xyz coordinates of ending points of tubes/cables
% Output:
%	conGtest:  the test results of gradients of constraints 


conGvec = gradientC(x);

% initialize a small matrix d with elements 1e-6
d = randn(size(x))*1e-6;
dvec = d(:);

conGtest = abs(constraintE(x+d)-constraintE(x-d)-2*(conGvec)*dvec)/norm(d);


end