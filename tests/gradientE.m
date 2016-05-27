function [gFixedPoints] = gradientE(x)
% function [gFixedPoints] = gradientE(x)

% Perform gradient calculations for the Energy function 
% Input:
%	x:   xyz coordinates of ending points of tubes/cables
%           (3*5)
% Intermediate variable:
%	gE:  the gradient of the system total energy
%           (24*1)
% Output:
%	gFixedPoints:	the corrected gradient taking care of fixed points
%                   (15*1)

global tube_pts cable_pts tube_masses cable_len fixed_pts;
global g Es fixed_x;

x = JoinFixedPoints(x, fixed_x);
number_of_cables = length(cable_pts);
number_of_pts = length(x);

gE = zeros(3*number_of_pts,1);

for i = 1:number_of_pts
    [t,~] = find(tube_pts==i);
    gE(3*i,1) = tube_masses(t)*g*0.5;
end

for c = 1:number_of_cables
    % d is the cable length after stretched
    xs = cable_pts(c,1);
    xe = cable_pts(c,2);
    xdiff = x(:,xs) - x(:,xe);
    d = norm(xdiff);

    if d - cable_len(c) > 0
        gradient = (Es/cable_len(c) - Es/d)*xdiff;
        gE(3*xs-2:3*xs,1) = gE(3*xs-2:3*xs,1) + gradient;
        gE(3*xe-2:3*xe,1) = gE(3*xe-2:3*xe,1) - gradient;
    end
end


% taking care of fixed points
% eliminate those entries in gradient vector
number_of_fixed_pts = length(fixed_pts);
number_of_pts_out = number_of_pts - number_of_fixed_pts;
gFixedPoints = zeros(3*number_of_pts_out,1);

j = 1;

for i = 1:number_of_pts
    if any(i==fixed_pts)
        continue;   % ignore the gradients of fixed points
    else
        gFixedPoints(3*j-2:3*j,1) = gE(3*i-2:3*i,1);
        j = j + 1;
    end
end

end