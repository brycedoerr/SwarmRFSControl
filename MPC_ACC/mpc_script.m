function mpc_script(x0,xdes,t,dt,dist_name,fh_cost)
%
% Function to simulate MPC via RFS formulation
%
% Synopsis:
%     mpc_script(x0,xdes,t,dt,dist_name,fh_cost)
%
% Input:
%     x0          =   states for each density (# states x # densities)
%     xdes        =   target states (# target states x # target densities)
%     t           = time sequence (column vector)
%     dist_name   =   string for labelling distances
%     fh_cost     = function handle for cost
%
%
% By: Bryce Doerr -- Aug. 2018

%% Initialize Dynamics with number of densities
tic
Z=[0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0;]; %build the A matrix
z=zeros(4,4); %zero matrix to build the A matrix
A=[Z z z z; z Z z z; z z Z z; z z z Z];
W=[0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0;];%build the B matrix
z1=zeros(4,4);%zero matrix to build the B matrix
B=[W z1 z1 z1; z1 W z1 z1; z1 z1 W z1; z1 z1 z1 W;];
V=[1 0 0 0; 0 1 0 0]; %Build C matrix0
z2=zeros(2,4);%zero matrix to build the C matrix
C=[V z2 z2 z2; z2 V z2 z2; z2 z2 V z2; z2 z2 z2 V];
D=zeros(8,16);
sys=ss(A,B,C,D);
%% Discrete
sys_d=c2d(sys, dt);
A=sys_d.A;
B=sys_d.B;

%% Simulation
R=.0000008*eye(16,16);

ushape=zeros(1,16);
u=repmat(ushape,length(t),1)';
u_sim=zeros(size(u,1),size(u,2));

nt=length(t);
xdi=reshape(xdes,4*size(xdes,2),1);
xdi=repmat(xdi,1,nt);
n_inputs=size(u,1);
x_sim2=zeros(16,nt);
x_sim2(:,1)=reshape(x0,4*4,1);
options = optimoptions('fminunc');

%MPC Simulation
for k=1:length(t)-1
    nt=3; %Prediction Horizon
    uout = fminunc(@objective,zeros(16*nt,1),options,reshape(x_sim2(:,k),4*4,1),A,B,n_inputs,nt,xdi,R,fh_cost);
    u_sim(:,k)=uout(1:16,1);
    x_sim2(:,k+1)=A*x_sim2(:,k)+B*u_sim(:,k);
end

x_sim2=zeros(length(A),length(t));
x_sim2(:,1)=reshape(x0,4*4,1);
for k=1:length(t)-1
    x_sim2(:,k+1)=A*x_sim2(:,k)+B*u_sim(:,k);    
end	
toc

%% Plot Figures
figure
subplot(4,1,1);plot(t,x_sim2(1,:),'-b','LineWidth',2);hold on;plot(t,1*x_sim2(1,:).^0,'--r')
plot(t,x_sim2(2,:),'.g','LineWidth',2);plot(t,1*x_sim2(1,:).^0,'--r')
xlabel('time (s)');ylabel('Pos. 1');grid on;

subplot(4,1,2);plot(t,x_sim2(5,:),'-b','LineWidth',2);hold on;plot(t,-1*x_sim2(1,:).^0,'--r')
plot(t,x_sim2(6,:),'.g','LineWidth',2);plot(t,1*x_sim2(1,:).^0,'--r')
xlabel('time (s)');ylabel('Pos. 2');grid on;

subplot(4,1,3);plot(t,x_sim2(9,:),'-b','LineWidth',2);hold on;plot(t,-1*x_sim2(1,:).^0,'--r')
plot(t,x_sim2(10,:),'.g','LineWidth',2);plot(t,-1*x_sim2(1,:).^0,'--r')
xlabel('time (s)');ylabel('Pos. 3');grid on;

subplot(4,1,4);
p1=plot(t,x_sim2(13,:),'-b','LineWidth',2);
hold on;
p2=plot(t,1*x_sim2(1,:).^0,'--r');
p3=plot(t,x_sim2(14,:),'.g','LineWidth',2);
plot(t,-1*x_sim2(1,:).^0,'--r')
xlabel('time (s)');ylabel('Pos. 4');grid on;

hl=legend([p1,p3,p2],{'X Pos.','Y Pos.','Target Destinations'},'location','southeast');
newPosition = [.75 0.063 0.0 0.0];
newUnits = 'normalized';
set(hl,'Position', newPosition,'Units', newUnits);

%% Save Results
num_density_des=size(xdes,2);
num_density=length(x0);
save(strcat('test_case_',dist_name,'_',int2str(num_density),int2str(num_density_des),'.mat'));