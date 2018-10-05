function cost = cost_nl_l2(u,xi,xdi,R)
%
% Function to determine cost using L2^2 Divergence
%
% Synopsis:
%     cost = cost_nl_l2(u,xi,xdi,R)
%
% Input:
%     u           =   control input for all densities(column vector)
%     xi          =   states for each density (# states x # densities)
%     xdi         =   target states (# target states x # target densities)
%     R           =   Control Weight Matrix (matrix size u)
%
% Output:
%     cost          =   scalar cost using L2^2 Divergence
%
%
% By: Bryce Doerr -- Aug. 2018

%Initialization
xj=xi;
xdj=xdi;
J1=0;
J1s=0;
J2=0;
S=0.5^2*diag([1 1]);

for i=1:length(xi)
    for j=1:length(xi)       
       J1=mvnpdf(xi(1:2,i), xj(1:2,j), 2*S)+J1;       
    end
end

for i=1:length(xi)
    for dj=1:size(xdj,2)
        J2=mvnpdf(xi(1:2,i), xdj(1:2,dj), 2*S)+J2;
        J1s=-1*(-1/2*(xi(1:2,i)-xdj(1:2,dj))'*inv(2*S)*(xi(1:2,i)-xdj(1:2,dj))-log(sqrt((2*pi)^(size(xi,1)/2)*det(2*S))))+J1s;
    end
end

%Output cost function
cost=600*(J1-2*J2)+u'*R*u;
