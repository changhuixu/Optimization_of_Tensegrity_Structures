function [C] = constraintE(x)
% function [C] = constraintE(x)

% Perform constraint calculations
% Input:
%	x:  xyz coordinates of ending points of tubes/cables
%       (3*5)
% Output:
%	C:  the array of constraints the system
%       (4*1)

global tube_pts tube_lengths fixed_x;

x = JoinFixedPoints(x, fixed_x);
number_of_tubes = length(tube_pts);

C = zeros(number_of_tubes, 1);

for i = 1:number_of_tubes
    xdiff = x(:,tube_pts(i,1)) - x(:,tube_pts(i,2));
    C(i,1) = norm(xdiff)^2 - tube_lengths(i)^2;
end

end