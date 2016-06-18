function PlotTensegrity(x)

global tube_pts cable_pts;

number_of_tubes = length(tube_pts);
number_of_cables = length(cable_pts);

figure

% Plot tubes
for i = 1:number_of_tubes
    xs = tube_pts(i,1);
    xe = tube_pts(i,2);
    pts = [x(:,xs) x(:,xe)];
    plot3(pts(1,:), pts(2,:), pts(3,:),'-b','LineWidth',5);
    axis equal;
    hold on;
end

% Plot cables
for i = 1:number_of_cables
    xs = cable_pts(i,1);
    xe = cable_pts(i,2);
    pts = [x(:,xs) x(:,xe)];
    plot3(pts(1,:), pts(2,:), pts(3,:),'--r','LineWidth',1);
    axis equal;
    hold on;
end

hold off;
% choose a view angle to show the Tensegrity structure
view([133,2]);
title('Tensegrity Structure Optimization Result', 'FontSize', 18);

end
