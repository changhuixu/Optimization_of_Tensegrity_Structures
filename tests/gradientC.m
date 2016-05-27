function [gradientC] = gradientC(x)
% function [gradientC] = gradientC(x)

% Perform gradient of constraint functions calculation
% Input:
%	x:              xyz coordinates of ending points of tubes/cables
%                   (3*5)
% Output:
%	gradientC:      the matrix of gradients of constraints the system
%                   (4*15)

global tube_pts fixed_pts fixed_x;

x = JoinFixedPoints(x, fixed_x);
number_of_tubes = length(tube_pts);
number_of_fixed_pts = length(fixed_pts);
number_of_pts = length(x);
number_of_pts_out = number_of_pts - number_of_fixed_pts;
conGradient = zeros(number_of_tubes, 3*number_of_pts);
gradientC = zeros(number_of_tubes, 3*number_of_pts_out);

for i = 1:number_of_tubes
    xs = tube_pts(i,1);
    xe = tube_pts(i,2);
    xdiff = x(:,xs) - x(:,xe);
    conGradient(i, 3*xs-2:3*xs) = 2*xdiff;
    conGradient(i, 3*xe-2:3*xe) = -2*xdiff;
end

j = 1;
for i = 1:number_of_pts
    if any(i==fixed_pts)
        continue;   % ignore the gradients of fixed points
    else
        gradientC(:,3*j-2:3*j) = conGradient(:,3*i-2:3*i);
        j = j + 1;
    end
end

end