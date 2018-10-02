function cost=objective(u,x0,A,B,n_inputs,nt,xdi,R,fh_cost)
%
% Function to get cost for entire trajectory in time
%
% Synopsis:
%     cost=objective(u,x0,A,B,n_inputs,nt,xdi,R,fh_cost)
%
% Input:
%     u           =   control input for all densities for prediction horizon(column vector)
%     x0          =   states for each density (# states x # densities)
%     A           =   Block A matrix for densities
%     B           =   Block B matrix for densities
%     n_inputs    =   # of control inputs
%     nt          =   length of time
%     xdi        =   target states (# target states x # target densities)
%     R           =   Control Weight Matrix (matrix size u)
%     fh_cost     = function handle for cost
%
% Output:
%     cost          =   scalar cost through time
%
% By: Bryce Doerr -- Aug. 2018

%Initialization
x_sim2=zeros(4*4,nt);
u_sim=reshape(u,n_inputs,nt);
x_sim2(:,1)=x0;
cost=0;

%Determine cost from k=1
for k=1:nt-1
    cost=cost+fh_cost(u_sim(:,k),reshape(x_sim2(:,k),4,4),reshape(xdi(:,k),4,size(xdi,1)/4),R);
    x_sim2(:,k+1)=A*x_sim2(:,k)+B*u_sim(:,k);
end	




