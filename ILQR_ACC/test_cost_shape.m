function test_cost_shape(x,xdes,size_ctrl,fh_cost)
%
% Function to plot surface plot of different divergences
%
% Synopsis:
%     test_cost_shape(x,xdes,size_ctrl,fh_cost)
%
% Input:
%     x          =   current intensities (# states x # current densities)
%     xdes       =   target states (# target states x # target densities)
%     size_ctrl  =   length of control sequence
%     fh_cost          =   function handle for cost
%
%
% By: Bryce Doerr -- Aug. 2018

%Initialization
x0=x;
size_u=size_ctrl*size(x,2);
R=eye(size_u);

%Create a mesh grid [-12,12] by [-12, 12]
[X,Y]=meshgrid(linspace(-12,12,50),linspace(-12,12,50));
cost_out=zeros(50,50);
for i=1:50
    for j=1:50
        %Move a point throughout the mesh grid with all other densities as
        %stationary
        x0(1,1)=X(i,j);
        x0(2,1)=Y(i,j);
        cost_out(i,j) = fh_cost(zeros(size_u,1),x0,xdes,R);
    end
end

%Plot Surface Plot
figure; 
surf(cost_out);hold on;
xlabel('x')
ylabel('y')
title('Position Shape');

