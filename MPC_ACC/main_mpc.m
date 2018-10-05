%Main MPC (Run this one)
%
% Main script to run MPC code
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
clear all;

%Config
dist_name = 'l2q'; % can be l2, l2q, cs
dt=.008;
t=[0:dt:.4];

x0=[3 3 2 2; -3 3 -2 2; -3 -3 -2 -2; 3 -3 2 -2]';

%Change # of target densities
xdes=[1 1 0 0; -1 1 0 0; -1 -1 0 0;1 -1 0 0;]'; %case1
%  xdes=[1 1 0 0; -1 1 0 0; -1 -1 0 0;]';%case2
% xdes=[1 1 0 0; -1 1 0 0; -1 -1 0 0;1 -1 0 0; 0 0 0 0]';%case3

num_density_des=length(xdes);
num_density=length(x0);

%Determine Appropriate Grad_Hess for different objective functions
if strcmp(dist_name,'l2q')
        fh_cost = str2func(strcat('cost_nl_',dist_name));

        
elseif strcmp(dist_name,'l2')
        fh_cost = str2func(strcat('cost_nl_',dist_name));

        
elseif strcmp(dist_name,'cs')
        fh_cost = str2func(strcat('cost_nl_',dist_name));


end

%Plot surface plot of different distances
test_cost_shape(x0,xdes,fh_cost)

%Determine MPC response
mpc_script(x0,xdes,t,dt,dist_name,fh_cost)

%Plot GM of ILQR response
contour_map(x0,xdes,dist_name)




