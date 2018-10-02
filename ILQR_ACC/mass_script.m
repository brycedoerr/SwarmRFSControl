function mass_script(m0,t,x0,xdes,dist_name)

%Mass output%
% Function to plot control and change in mass of ILQR via RFS formulation
%
% Synopsis:
%     mass_script(m0,t,x0,xdes,dist_name)
%
% Input:
%     m0         = Initial mass of each agent
%     t          = time sequnce (vector)
%     x0          =   states for each density (# states x # densities)
%     xdes       =   target states (# target states x # target densities)
%     dist_name          =   string for labelling distances
%
%
% By: Bryce Doerr -- Aug. 2018


%Initialize Densities
num_density_des=size(xdes,2);
num_density=size(x0,2);

%Load Data
load(strcat('test_case_',dist_name,'_',int2str(num_density),int2str(num_density_des),'.mat'));

%Set number of densities you want to plot
num_ctrl_plot=5;

%ODE 45 options
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);

%Initialize mass and acceleration
m=zeros(length(t),num_ctrl_plot);
m(1,:)=m0;%kg
u=zeros(num_ctrl_plot,length(t));

%Simulate the mass and acceleration response
for j=1:num_ctrl_plot
    u(j,:)=sqrt(u_bar(3*j-2,:).^2+u_bar(3*j-1,:).^2+u_bar(3*j,:).^2);
    for i=2:length(t)
        [~,x45]=ode45(@fn_mass,[t(i-1),t(i)],m(i-1,j)',opts,u(j,i-1));
        m(i,j)=x45(end,:);
    end
end


%Plot Figures
figure
subplot(2,1,1)
hold on;
grid on;
for j=1:num_ctrl_plot
    plot(t,u(j,:),'LineWidth',2)
end

xlabel('Time (s)')
ylabel('Acceleration (m/s)')

subplot(2,1,2)
hold on;
grid on;
for j=1:num_ctrl_plot
    plot(t,m(:,j),'LineWidth',2)
end

xlabel('Time (s)')
ylabel('Mass (kg)')
end