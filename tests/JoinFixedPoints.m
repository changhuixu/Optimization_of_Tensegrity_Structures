function [x_out] = JoinFixedPoints(x,fixed_x)
% function [x_out] = JoinFixedPoints(x,fixed_x)

% Ignore fixed coordinates, don't consider them in future optimization
% Input:
%	x_fixed:	the coordinates not taking account of fixed points
%                   (3*3)
%	x:          the coordinates not taking account of fixed points
%                   (3*5)
% Output:
%	x_out:      the coordinates of all points
%                   (3*8)

global fixed_pts;

number_of_pts = length(x);
number_of_fixed_pts = length(fixed_pts);
number_of_pts_out = number_of_pts + number_of_fixed_pts;

x_out = zeros(3, number_of_pts_out);
j = 1;
k = 1;

for i = 1:number_of_pts_out
    if any(i==fixed_pts)
        x_out(:,i) = fixed_x(:,j);
        j = j + 1;
    else
        x_out(:,i) = x(:,k);
        k = k + 1;
    end
end

end