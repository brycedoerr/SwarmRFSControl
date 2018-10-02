function cost = cost_nl_cs(u,xi,xdi,R)
%
% Function to determine cost using Cauchy-Schwarz Divergence
%
% Synopsis:
%     cost = cost_nl_cs(u,xi,xdi,R)
%
% Input:
%     u           =   control input for all densities(column vector)
%     xi          =   states for each density (# states x # densities)
%     xdi         =   target states (# target states x # target densities)
%     R           =   Control Weight Matrix (matrix size u)
%
% Output:
%     cost          =   scalar cost using Cauchy-Schwarz Divergence
%
%
% By: Bryce Doerr -- Aug. 2018

%Initialization
xj=xi;
xdj=xdi;
J1=0;
J2=0;
S=0.5^2*diag([1 1 1 1]);

for i=1:length(xi)
    for j=1:length(xi)
       J1=mvnpdf(xi(:,i), xj(:,j), 2*S)+J1;
    end
end

for i=1:length(xi)
    for dj=1:length(xdj)
        J2=mvnpdf(xi(:,i), xdj(:,dj), 2*S)+J2;
    end
end

for di=1:length(xdi)
    for dj=1:length(xdj)
    end
end

%Output cost function
cost=1/2*log(J1)-log(J2)+u'*R*u;
