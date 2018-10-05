%Main ILQR (Run this one)
%
% Main script to run ILQR code
% By: Bryce Doerr -- Aug. 2018
% Copyright <2018> <Bryce Doerr>
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

clc;
close all;
clear;

%Config
dist_name = 'l2q'; % can be l2, l2q, cs
dt=1;
t=0:dt:40;

%Satellite Parameters
n=.00113;%ang freq
m0=1;%kg

%Make a star shape
[xdes,num_density_des]=star(9);
xdes=10*xdes;

%Rotate Star Shape
for i=1:length(t)-1
    for j=1:num_density_des
    xdes(1:2,j,i+1)=[cos(n*t(i+1)) -sin(n*t(i+1)); sin(n*t(i+1)) cos(n*t(i+1))]*xdes(1:2,j,i);
    end
end

%Initial Conditions for densitiees
x0=-1 + (1+1)*rand(num_density_des,1);
y0=-1 + (1+1)*rand(num_density_des,1);
z0=zeros(num_density_des,1);
xdot0=rand(num_density_des,1)*0;
ydot0=-2*n*x0;
zdot=zeros(num_density_des,1);
x0=[x0 y0 z0 xdot0 ydot0 zdot]';

%Size of Densities
num_density=size(x0,2);
size_state=size(x0,1);
size_ctrl=3;

%Determine Appropriate Grad_Hess for different objective functions
if strcmp(dist_name,'l2q')
        fh_cost = str2func(strcat('cost_nl_',dist_name));
        fh_cost_grad_hess=str2func(strcat('dj1'));


end
%Switch between using an initial guess or initialize to zero for control
%input
u_ctrl=CWH_dyn_example(x0,xdes,t,dt,size_ctrl);
% u_ctrl=zeros(size_ctrl*size(x0,2),length(t));

%Plot surface plot of different distances
test_cost_shape(x0,xdes,size_ctrl,fh_cost)

%Determine ILQR Response
ilqr_script(x0,xdes,t,dt,u_ctrl,fh_cost,fh_cost_grad_hess,dist_name)

%Plot response
plot_map(x0,xdes,dist_name)

%Plot GM of ILQR response
contour_map(x0,xdes,dist_name)

%Plot control input and mass loss
mass_script(m0,t,x0,xdes,dist_name)




