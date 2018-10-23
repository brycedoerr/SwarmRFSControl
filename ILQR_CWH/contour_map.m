function contour_map(x0,xdes,dist_name)
%
% Function to plot Gaussian Mixtures of ILQR via RFS formulation
%
% Synopsis:
%     contour_map(x0,xdes,dist_name)
%
% Input:
%     x0          =   states for each density (# states x # densities)
%     xdes       =   target states (# target states x # target densities)
%     dist_name          =   string for labelling distances
%
%
% By: Bryce Doerr -- Aug. 2018

%Initialize Density
num_density_des=size(xdes,2);
num_density=size(x0,2);

%Load Data
load(strcat('test_case_',dist_name,'_',int2str(num_density),int2str(num_density_des),'.mat'));

x_sim2=x_bar;
u_sim=u_bar;
N=400; %number of samples
ts=length(t); %number of time steps

D=zeros(num_density,1);

%Create a mesh grid [-15,15] by [-15,15]
[X,Y]=meshgrid(linspace(-15, 15,100),linspace(-15, 15,100));
k_count=1;

%Initialize Statistics for GM
mu=zeros(num_density,2);
sigma=zeros(2,2,num_density);
for i=1:num_density
    mu(i,:) = x_sim2(6*i-5:6*i-4,1)';
    sigma(:,:,i) = 0.03*diag([1 1]);
end

%Create GM distribution
p = ones(1,num_density)/2;
obj = gmdistribution(mu,sigma,p);
points=random(obj,N);
points_sim=zeros(6,ts,N);

%Initial Conditions for each sample
for k=1:N 
   points_sim(:,1,k)=[points(k,:)';0;0;0;0];
end

for j=1:ts-1
    for i=1:num_density
        mu(i,:) = x_sim2(6*i-5:6*i-4,j);
        sigma(:,:,i) =  0.03*diag([1 1]);
    end
    
    p = ones(1,num_density)/2;
    obj = gmdistribution(mu,sigma,p);
    Fun=@(x,y)pdf(obj,[x y]);
    PDFj=reshape(Fun(reshape(X,100*100,1),reshape(Y,100*100,1)),100,100);
    
    %Loop through samples
    for k=1:N 
        for i=1:num_density
           D(i)= sqrt((points_sim(1:2,j,k)-mu(i,:)')'*inv(sigma(:,:,1))*(points_sim(1:2,j,k)-mu(i,:)'));
        end
        [~,I]=min(D);
        points_sim(:,j+1,k)=A(6*I-5:6*I,6*I-5:6*I)*points_sim(:,j,k)+B(6*I-5:6*I,3*I-2:3*I)*u_sim(3*I-2:3*I,j);
    end
    
    if j==1
        PDFj1=PDFj; %save for time=0
    elseif j==16
        PDFj16=PDFj;
    elseif j==31
        PDFj31=PDFj;
    end    

    %Create GIF
    figure (25)
    contour(X,Y,PDFj,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x (m)')
    ylabel('y (m)')
    title('Clohessy-Wiltshire Relative Motion Control')
    grid on;
    for i=1:num_density_des
        plot(xdes(1,i,j),xdes(2,i,j),'xk','LineWidth',1)
    end
    xlim([-15 15]);
    ylim([-15 15]);
    hold off;

    mm(:,k_count) = getframe(gcf);
    k_count=k_count+1;
    pause(1);
end

%Save contours for Final Time
for i=1:num_density
        mu(i,:) = x_sim2(6*i-5:6*i-4,end);
        sigma(:,:,i) =  0.03*diag([1 1]);
end
p = ones(1,num_density)/2;
obj = gmdistribution(mu,sigma,p);
Fun=@(x,y)pdf(obj,[x y]);
PDFj=reshape(Fun(reshape(X,100*100,1),reshape(Y,100*100,1)),100,100);
PDFj41=PDFj;

%Plot Figures
figure
subplot(2,2,1);
contour(X,Y,PDFj1,linspace(0.006,0.218,100),'linewidth',2)
hold on;
set(gcf,'color',[1 1 1]);
xlabel('x')
ylabel('y')
title('time = 0 s')
grid on;
for i=1:num_density_des
    plot(xdes(1,i,1),xdes(2,i,1),'xk','LineWidth',1)
end
xlim([-15 15]);
ylim([-15 15]);
hold off;

subplot(2,2,2);
contour(X,Y,PDFj16,linspace(0.006,0.218,100),'linewidth',2)
hold on;
set(gcf,'color',[1 1 1]);
xlabel('x')
ylabel('y')
title('time = 15 s')
grid on;
for i=1:num_density_des
    plot(xdes(1,i,16),xdes(2,i,16),'xk','LineWidth',1)
end
xlim([-15 15]);
ylim([-15 15]);
hold off;

subplot(2,2,3);
contour(X,Y,PDFj31,linspace(0.006,0.218,100),'linewidth',2)
hold on;
set(gcf,'color',[1 1 1]);
xlabel('x')
ylabel('y')
grid on;
title('time = 30 s')
for i=1:num_density_des
    plot(xdes(1,i,31),xdes(2,i,31),'xk','LineWidth',1)
end
xlim([-15 15]);
ylim([-15 15]);
hold off;

subplot(2,2,4);
contour(X,Y,PDFj41,linspace(0.006,0.218,100),'linewidth',2)
hold on;
set(gcf,'color',[1 1 1]);
xlabel('x')
ylabel('y')
grid on;
title('time = 40 s')
for i=1:num_density_des
    plot(xdes(1,i,41),xdes(2,i,41),'xk','LineWidth',1)
end
xlim([-15 15]);
ylim([-15 15]);
hold off;

hp4 = get(subplot(2,2,4),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.00  hp4(2)+.015  0.03  hp4(2)+hp4(3)*2.1-.02])

 
    
%Save GIF    
movie2gif(mm,'ggg3.gif')