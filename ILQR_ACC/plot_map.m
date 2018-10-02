function plot_map(x0,xdes,dist_name)
%
% Function to plot the trajectory of ILQR via RFS formulation
%
% Synopsis:
%     plot_map(x0,xdes,dist_name)
%
% Input:
%     x0          =   states for each density (# states x # densities)
%     xdes       =   target states (# target states x # target densities)
%     dist_name          =   string for labelling distances
%
%
% By: Bryce Doerr -- Aug. 2018

%Initialize Size
num_density_des=size(xdes,2);
num_density=size(x0,2);

%Load Data
load(strcat('test_case_',dist_name,'_',int2str(num_density),int2str(num_density_des),'.mat'));

%Plot Figures
figure
hold on;
xlim([-15 15]);
ylim([-15 15]);
for i=1:num_density
    plot(x_bar(6*i-5,:),x_bar(6*i-4,:))
    plot(x_bar(6*i-5,1),x_bar(6*i-4,1),'or');
    plot(x_bar(6*i-5,end),x_bar(6*i-4,end),'*g');
end
plot(xdes(1,:),xdes(2,:),'xb')
xlabel('X')
ylabel('Y')
sum(u_new,2)

figure
subplot(2,2,1);
hold on;
for i=1:num_density
    plot(x_bar(6*i-5,1),x_bar(6*i-4,1),'or');
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
    plot(x_bar(6*i-5,16),x_bar(6*i-4,16),'or');
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
    plot(x_bar(6*i-5,31),x_bar(6*i-4,31),'or');
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
    plot(x_bar(6*i-5,41),x_bar(6*i-4,41),'or');
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

end