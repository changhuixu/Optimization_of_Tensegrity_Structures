function [x_out,x_fixed] = FixedPoints(x)
% function [x_out,x_fixed] = FixedPoints(x)

% Ignore fixed coordinates, don't consider them in future optimization
% Input:
%	x:          the coordinates of all points
%                   (3*8)
% Output:
%	x_out:      the coordinates not taking account of fixed points
%                   (3*5)
%	x_fixed:	the coordinates not taking account of fixed points
%                   (3*3)

global fixed_pts;

number_of_pts = length(x);
number_of_fixed_pts = length(fixed_pts);
number_of_pts_out = number_of_pts - number_of_fixed_pts;

x_out   = zeros(3, number_of_pts_out);
x_fixed = zeros(3, number_of_fixed_pts);

j = 1;
k = 1;
for i = 1:number_of_pts
    if any(i==fixed_pts)
        x_fixed(:,k) = x(:,i);
        k = k + 1;
    else
        x_out(:,j) = x(:,i);
        j = j + 1;
    end
end

if j ~= number_of_pts_out + 1 || k ~= number_of_fixed_pts + 1
    fprintf('\nError: Check the number of total/fixed points.\n\n');
    return
end

end