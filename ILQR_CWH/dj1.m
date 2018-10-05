function [cx, cxx, cu, cuu, cux]= dj1(x,xdes,u,R)
%
% Function to determine grad & hess of cost using L2^2 + Quad Divergence
%
% Synopsis:
%     [cx, cxx, cu, cuu, cux]= dj1(x,xdes,u,R)
%
% Input:
%     x          =   states for each density (# states x # densities)
%     xdes         =   target states (# target states x # target densities)
%     u           =   control input for all densities(column vector)
%     R           =   Control Weight Matrix (matrix size u)
%
% Output:
%     cx          =   gradient of cost in terms of x
%     cxx          =   hessian of cost in terms of xx
%     cu          =   gradient of cost in terms of u
%     cuu          =   hessian of cost in terms of uu
%     cux          =   hessian of cost in terms of ux
%
% By: Bryce Doerr -- Aug. 2018

%Obtain Gradient and Hessian for L2Q in cost_nl_l2q
x = transpose(x);
y = transpose(xdes);
n = size(x, 1); %size distribution
t = size(x, 2); %size state
h = size(y, 1); %size target distribution
S=0.5^2*diag(ones(1, t));
Q = inv(2*S);

%Initialize Gradient and Hessian
DJ1 = zeros(size(x));
DDJ1 = zeros(size(x,1), size(x, 2), size(x,1), size(x,2));

DJ2 = zeros(size(x));
DDJ2 = zeros(size(x,1), size(x, 2), size(x,1), size(x,2));

DJ3 = zeros(size(x));
DDJ3 = zeros(size(x,1), size(x, 2), size(x,1), size(x,2));

%Loop through Gradient and Hessians for number of densities
for i = 1:n
    Qxi = Q*transpose(x(i, :));
    Qxi_xi = x(i, :)*Qxi;
    for j = 1:n
        Qxj = Q*transpose(x(j, :));
        Rij = Qxi_xi - x(i, :)*Qxj - x(j, :)*Qxi + x(j, :)*Qxj;
        Rij = exp(-1/2*Rij);

        chi_ijyz = zeros(n, t);
        
        for l = 1:n
            for m = 1:t
            term21 = 0;
                for p = 1:t
                    term21 = term21 - Q(m, p)*x(j, p) - Q(p, m)*x(j, p) + Q(m, p)*x(i, p) + Q(p, m)*x(i, p);
                end
                chi_ijyz(l, m) = (term21*eq(i, l) - term21*eq(j, l));
            end
        end
        for r = 1:n
            for s = 1:t
                DJ_contrib = -1/2*Rij*chi_ijyz(r, s);
                DJ1(r, s) = DJ1(r, s) + DJ_contrib;
                for l = 1:n
                    for m = 1:t
                        eta = 2*Q(s,m)*eq(i, l)*eq(i, r) + ... 
                              2*Q(s,m)*eq(j, l)*eq(j, r) - ... 
                              2*Q(s, m)*eq(j, l)*eq(i, r)  - ...
                              2*Q(s,m)*eq(i, l)*eq(j, r); 

                        DDJ1(r, s, l, m) = DDJ1(r, s, l, m) + 1/4*Rij*chi_ijyz(l, m)*chi_ijyz(r,s) + -1/2*Rij*eta;
                    end
                end
            end
        end        
    end
end

Qyk_yk = zeros(1, h);
for k = 1:h
    Qyk_yk(k) = y(k, :)*(Q*transpose(y(k, :)));
end
        
for i = 1:n
    Qxi = Q*transpose(x(i, :));
    Qxi_xi = x(i, :)*Qxi;
    for k = 1:h
        Qyk = Q*transpose(y(k, :));
        
        Rik = Qxi_xi - y(k, :)*Qxi - x(i, :)*Qyk + Qyk_yk(k);
        Rik = exp(-1/2*Rik);
        chi_ikyz = zeros(n, t);
        
        for l = 1:n
            for m = 1:t
                term23 = 0;
                for p = 1:t
                    term23 = term23 + Q(p, m)*x(i, p)+Q(m, p)*x(i, p) - Q(p, m)*y(k, p) - Q(m, p)*y(k, p);
                end
                chi_ikyz(l, m) = term23*eq(i, l);
            end
        end
               
   
        for r = 1:n
            for s = 1:t
                DJ2_contrib = -1/2*Rik*chi_ikyz(r, s);
                DJ2(r, s) = DJ2(r, s) + DJ2_contrib;
                DJ3_contrib = -1/2*chi_ikyz(r, s);
                DJ3(r, s) = DJ3(r, s) + DJ3_contrib;
                for l = 1:n
                    for m = 1:t
                        eta_ikrslm = 2*Q(m, s)*eq(i, l)*eq(i, r);
                        DDJ2(r, s, l, m) = DDJ2(r, s, l, m) + 1/4*Rik*chi_ikyz(l, m)*chi_ikyz(r, s) + -1/2*Rik*eta_ikrslm;
                        DDJ3(r, s, l, m) = DDJ3(r, s, l, m) + -1/2*eta_ikrslm;
                    end
                end
            end
        end
    end
end

multVal = (1/sqrt((2*pi)^t*det(2*S)));
DDJ1 = DDJ1*multVal;
DDJ2 = DDJ2*multVal;
DJ1 = DJ1*multVal;
DJ2 = DJ2*multVal;

cx1 = reshape(transpose(DJ1), [], 1);
cx2 = reshape(transpose(DJ2), [], 1);
cx3 = reshape(transpose(DJ3), [], 1);
cx=300*(cx1)-1200*(2*cx2)-(1/30)*cx3;%Gradient cx

cxx1 = zeros(length(cx1), length(cx1));
cxx2 = zeros(length(cx1), length(cx1));
cxx3 = zeros(length(cx1), length(cx1));

count = 0;
for i = 1:n
    for j = 1:t
        count = count + 1;
        cxx1(:, count) = reshape(transpose(DDJ1(:, :, i, j)), [], 1);
        cxx2(:, count) = reshape(transpose(DDJ2(:, :, i, j)), [], 1);
        cxx3(:, count) = reshape(transpose(DDJ3(:, :, i, j)), [], 1);
    end
end
cxx=300*(cxx1)-1200*(2*cxx2)-(1/30)*cxx3;%Hessian cxx

%Other Gradients and Hessians needed
cu=2*R*u;
cuu=2*R;
cux=zeros(size(R,1),size(x,1)*size(x,2));

end