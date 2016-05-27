function [hessianC] = hessianC(x)
% function [hessianC] = hessianC(x)

% Perform Hessian calculation of constraint functions 
% Input:
%	x:              xyz coordinates of ending points of tubes/cables
%                   (3*5)
% Output:
%	hessianC:       the matrix of Hessiants of constraints the system
%                   (15*15*4)

global tube_pts fixed_pts fixed_x;
x = JoinFixedPoints(x, fixed_x);
number_of_tubes = length(tube_pts);
number_of_pts = length(x);
number_of_fixed_pts = length(fixed_pts);
number_of_pts_out = number_of_pts - number_of_fixed_pts;

conHessian = zeros(3*number_of_pts, 3*number_of_pts, number_of_tubes);
tempHess = zeros(3*number_of_pts_out, 3*number_of_pts, number_of_tubes);
hessianC = zeros(3*number_of_pts_out, 3*number_of_pts_out, number_of_tubes);

for i = 1:number_of_tubes
    xs = tube_pts(i,1);
    xe = tube_pts(i,2);
    conHessian(3*xs-2:3*xs, 3*xs-2:3*xs, i) = ...
        conHessian(3*xs-2:3*xs, 3*xs-2:3*xs, i) + 2* eye(3);
    conHessian(3*xe-2:3*xe, 3*xe-2:3*xe, i) = ...
        conHessian(3*xe-2:3*xe, 3*xe-2:3*xe, i) + 2* eye(3);
    conHessian(3*xs-2:3*xs, 3*xe-2:3*xe, i) = ...
        conHessian(3*xs-2:3*xs, 3*xe-2:3*xe, i) - 2* eye(3);
    conHessian(3*xe-2:3*xe, 3*xs-2:3*xs, i) = ...
        conHessian(3*xe-2:3*xe, 3*xs-2:3*xs, i) - 2* eye(3);
end


% ignore the hessian of fixed points
for t = 1:number_of_tubes
    j = 1;
    for i = 1:number_of_pts
        if any(i==fixed_pts)
            continue;   % ignore the row of hessian of fixed points
        else
            tempHess(3*j-2:3*j,:,t) = conHessian(3*i-2:3*i,:,t);
            j = j + 1;
        end
    end

    j = 1;
    for i = 1:number_of_pts
        if any(i==fixed_pts)
            continue;   % ignore the col of hessian of fixed points
        else
            hessianC(:,3*j-2:3*j,t) = tempHess(:,3*i-2:3*i,t);
            j = j + 1;
        end
    end
end

end