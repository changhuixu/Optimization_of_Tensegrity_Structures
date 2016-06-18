function [hFixedPoints] = hessianE(x)
% function [hFixedPoints] = hessianE(x)

% Perform hessian calculations for the Energy function 
% Input:
%	x:   xyz coordinates of ending points of tubes/cables
%           (3*5)
% Intermediate Variable:
%	hE:  the hessian of the system total energy
%           (24*24)
% Output:
%	hFixedPoints:	the corrected Hessian taking care of fixed points
%                   (15*15)

global cable_pts cable_len fixed_pts;
global Es fixed_x;

x = JoinFixedPoints(x, fixed_x);
number_of_cables = length(cable_pts);
number_of_pts = length(x);

hE = zeros(number_of_pts*3,number_of_pts*3);

for c = 1:number_of_cables
    % d is the cable length after stretched
    xs = cable_pts(c,1);
    xe = cable_pts(c,2);
    xdiff = x(:,xs) - x(:,xe);
    d = norm(xdiff);

    if d - cable_len(c) > 0
        hess = (Es/cable_len(c)-Es/d)*eye(3)+Es*(xdiff)*(xdiff).'/d^3;
        hE(3*xs-2:3*xs, 3*xs-2:3*xs) = hE(3*xs-2:3*xs, 3*xs-2:3*xs) + hess;
        hE(3*xe-2:3*xe, 3*xe-2:3*xe) = hE(3*xe-2:3*xe, 3*xe-2:3*xe) + hess;
        hE(3*xs-2:3*xs, 3*xe-2:3*xe) = hE(3*xs-2:3*xs, 3*xe-2:3*xe) - hess;
        hE(3*xe-2:3*xe, 3*xs-2:3*xs) = hE(3*xe-2:3*xe, 3*xs-2:3*xs) - hess;
    end
end


number_of_fixed_pts = length(fixed_pts);
number_of_pts_out = number_of_pts - number_of_fixed_pts;
hFixedPoints = zeros(3*number_of_pts_out,3*number_of_pts_out);
tempHess = zeros(3*number_of_pts_out,3*number_of_pts);

% ignore the hessian of fixed points
j = 1;
for i = 1:number_of_pts
    if any(i==fixed_pts)
        continue;   % ignore the row of hessian of fixed points
    else
        tempHess(3*j-2:3*j,:) = hE(3*i-2:3*i,:);
        j = j + 1;
    end
end

j = 1;
for i = 1:number_of_pts
    if any(i==fixed_pts)
        continue;   % ignore the row of hessian of fixed points
    else
        hFixedPoints(:,3*j-2:3*j) = tempHess(:,3*i-2:3*i);
        j = j + 1;
    end
end


end