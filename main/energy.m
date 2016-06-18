function [E] = energy(x)
% function [E] = energy(x)

% Perform Energy calculation for the whole system
% Input:
%	x:  xyz coordinates of ending points of tubes/cables
%       (3*5)
% Output:
%	E:  the total energy of the system
%       (1*1) scalar

global tube_pts cable_pts tube_masses cable_len;
global e3 g Es fixed_x;

x = JoinFixedPoints(x, fixed_x);

number_of_tubes = length(tube_pts);
number_of_cables = length(cable_pts);

% calculate gravitational energy of tubes
E1 = 0;
for t = 1:number_of_tubes
    E1 = E1 + tube_masses(t)*(e3.')* ...
        (x(:,tube_pts(t,1)) + x(:,tube_pts(t,2)));
end
E1 = E1*g*0.5;

% calculate elastic energy of cables
E2 = 0;
for c = 1:number_of_cables
    % d is the cable length after stretched
    xs = cable_pts(c,1);
    xe = cable_pts(c,2);
    xdiff = x(:,xs) - x(:,xe);
    d = norm(xdiff);
    l = d - cable_len(c);
    if l > 0
        E2 = E2 + (l^2)/cable_len(c);
    end
end
E2 = E2*Es*0.5;


% total energy E 
E = E1 + E2;

end
