function [conHtest] = hCtest(x)
% function [conHtest] = hCtest(x)

% Perform tests for Hessians of constraint functions
% Input:
%	x:  xyz coordinates of ending points of tubes/cables
% Output:
%	conHtest:  the test results of hessians of constraints 

global tube_pts;
number_of_tubes = length(tube_pts);

conHtest = zeros(number_of_tubes,1);
conHvec  = hessianC(x);

% initialize a small matrix d with elements 1e-6
d = randn(size(x))*1e-6;
dvec = d(:);

conGvec1 = gradientC(x+d);
conGvec2 = gradientC(x-d);

for i = 1:number_of_tubes
    temp = conGvec1(i,:)'-conGvec2(i,:)'-2*conHvec(:,:,i)*dvec;
	conHtest(i) = norm(temp)/norm(d);
end

end