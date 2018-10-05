function u_ctrl= CWH_dyn_example(x0,xdes,t,dt,size_ctrl)
%
% Function to find a good initial guess for ILQR using LQ integrater
%
% Synopsis:
%     u_ctrl= CWH_dyn_example(x0,xdes,t,dt,size_ctrl)
%
% Input:
%     x0          =   states for each density (# states x # densities)
%     xdes       =   target states (# target states x # target densities)
%     t           =  Simulation time for trajectory (vector)
%     dt          =  time-step for simulation
%     size_ctrl   =  size of the control input
%
% By: Bryce Doerr -- Aug. 2018

%Clohessy Wiltshire

%Size of densities
num_density=size(x0,2);

%Orbital Parameters
n=.00113;%ang freq
m=1;%kg

%Build Dynamical System for number of densities
Z= [0     0 0 1 0   0;...
    0     0 0 0 1   0;...
    0     0 0 0 0   1;...
    3*n^2 0 0 0 2*n  0;...
    0    0 0 -2*n 0 0;...
    0 0 -n^2 0 0 0;]; %build the A matrix
A=kron(eye(size(x0,2)),Z);
W=[0 0 0; 0 0 0; 0 0 0; 1/m 0 0; 0 1/m 0; 0 0 1/m;];%build the B matrix
Br=repmat(W,1,size(x0,2));
Bc=mat2cell(Br, size(W,1),repmat(size(W,2),1,size(x0,2)));
B=blkdiag(Bc{:});
V=[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0]; %Build C matrix
Cr=repmat(V,1,size(x0,2));
Cc=mat2cell(Cr, size(V,1),repmat(size(V,2),1,size(x0,2)));
C=blkdiag(Cc{:});
D=zeros(length(A)/2,size(B,2));
sys=ss(A,B,C,D);

%% Discrete
sys_d=c2d(sys, dt);
A=sys_d.A;
B=sys_d.B;
C=sys_d.C;
D=sys_d.D;

%% Initialize Control
ushape=zeros(1,size_ctrl*size(x0,2));
u=repmat(ushape,length(t),1)';
u_ctrl=u;


%% LQR Response
x_ctrl=zeros(length(A),length(t));
x_ctrl(:,1)=reshape(x0,[],1);
xi=zeros(length(A)/2,length(t));
y=zeros(length(A)/2,length(t));

%LQR weights
Q=1*eye(size(A,1)+size(C,1));
R=.0000008*eye(3*size(x0,2));

%Find K Gains using LQI
K=lqi(sys_d,Q,R);

%Initialize Reference
r=xdes;
r(4:6,:,:)=[];

%Simulate
for i=1:length(t)-1
    u_ctrl(:,i)=-K*[x_ctrl(:,i);xi(:,i)];
    x_ctrl(:,i+1)=A*x_ctrl(:,i)+B*u_ctrl(:,i);
    y(:,i)=C*x_ctrl(:,i)+D*u_ctrl(:,i);
    
    xi(:,i+1)=xi(:,i)+dt*(reshape(r(:,:,i),[],1)-y(:,i));
end

%Plot Figures
figure
hold on;
xlim([-5 5]);
ylim([-5 5]);
for i=1:num_density
    plot(x_ctrl(6*i-5,:),x_ctrl(6*i-4,:))
    plot(x_ctrl(6*i-5,1),x_ctrl(6*i-4,1),'or');
    plot(x_ctrl(6*i-5,end),x_ctrl(6*i-4,end),'*g');
end
plot(xdes(1,:),xdes(2,:),'xb')
xlabel('X')
ylabel('Y')
xlim([-15 15]);
ylim([-15 15]);

figure
subplot(2,2,1);
hold on;
for i=1:num_density
    plot(x_ctrl(6*i-5,1),x_ctrl(6*i-4,1),'or');
end
xlabel('x')
ylabel('y')
title('time = 0 s')
grid on;
for i=1:num_density
    plot(xdes(1,i,1),xdes(2,i,1),'xk','LineWidth',2)
end
xlim([-15 15]);
ylim([-15 15]);

subplot(2,2,2);
hold on;
for i=1:num_density
    plot(x_ctrl(6*i-5,16),x_ctrl(6*i-4,16),'or');
end
xlabel('x')
ylabel('y')
title('time = 15 s')
grid on;
for i=1:num_density
    plot(xdes(1,i,16),xdes(2,i,16),'xk','LineWidth',2)
end
xlim([-15 15]);
ylim([-15 15]);

subplot(2,2,3);
hold on;
for i=1:num_density
    plot(x_ctrl(6*i-5,31),x_ctrl(6*i-4,31),'or');
end
xlabel('x')
ylabel('y')
title('time = 30 s')
grid on;
for i=1:num_density
    plot(xdes(1,i,31),xdes(2,i,31),'xk','LineWidth',2)
end
xlim([-15 15]);
ylim([-15 15]);

subplot(2,2,4);
hold on;
for i=1:num_density
    plot(x_ctrl(6*i-5,41),x_ctrl(6*i-4,41),'or');
end
xlabel('x')
ylabel('y')
title('time = 40 s')
grid on;
for i=1:num_density
    plot(xdes(1,i,41),xdes(2,i,41),'xk','LineWidth',2)
end
xlim([-15 15]);
ylim([-15 15]);

