function test_cost_shape(~,xdes,fh_cost)
%
% Function to plot surface plot of different divergences
%
% Synopsis:
%     test_cost_shape(~,xdes,fh_cost)
%
% Input:
%     xdes       =   target states (# target states x # target densities)
%     dist_name          =   string for labelling distances
%
%
% By: Bryce Doerr -- Aug. 2018

%Initialization
R=eye(16,16);

%Create a mesh grid [-12,12] by [-12,12]
[X,Y]=meshgrid(linspace(-12,12,50),linspace(-12,12,50));
cost_out=zeros(50,50);
for i=1:50
    for j=1:50
        %Move a point throughout the mesh grid with all other densities as
        %stationary
        x0=[X(i,j) Y(i,j) 0 0; -3 3 0 0; -3 -3 0 0; 3 -3 0 0]';
        cost_out(i,j) = fh_cost(zeros(16,1),reshape(x0,4,4),xdes,R);
    end
end

%Plot Surface Plot
figure
surf(cost_out);hold on;
xlabel('x')
ylabel('y')
title('Position Shape');
