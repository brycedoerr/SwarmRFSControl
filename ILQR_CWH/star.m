function [xdes,num_density_des] = star(pt_density)
%
% Function to create evenly spaced points for a star shape
%
% Synopsis:
%     [xdes,num_density_des] = star(pt_density)
%
% Input:
%     pt_density          =   density between each point in the star shape
%
% Output:
%     xdes       =   target states (# target states x # target densities)
%     num_density_des = # target densities
%
% By: Bryce Doerr -- Aug. 2018
% Adapted from Sattari:
% https://www.mathworks.com/matlabcentral/fileexchange/41454-grid-of-points-within-a-polygon

%star
xv = [0.5;0.2;1.0;0;0.8;0.5]*2-1;
yv = [1.0;0.1;0.7;0.7;0.1;1]*2-1;

N = sqrt(pt_density);

%Find the bounding rectangle
lower_x = min(xv);
higher_x = max(xv);

lower_y = min(yv);
higher_y = max(yv);

%Create a grid of points within the bounding rectangle
inc_x = 1/N;
inc_y = 1/N;

interval_x = lower_x:inc_x:higher_x;
interval_y = lower_y:inc_y:higher_y;
[bigGridX, bigGridY] = meshgrid(interval_x, interval_y);
	
%Filter grid to get only points in polygon
in = inpolygon(bigGridX(:), bigGridY(:), xv, yv);

%Return the co-ordinates of the points that are in the polygon
inPoints = [bigGridX(in), bigGridY(in)];
    

figure (1)
plot(inPoints(:,1),inPoints(:,2),'.r','LineWidth',2)

num_density_des=length(inPoints);

xdes=[inPoints,zeros(num_density_des,4)]';