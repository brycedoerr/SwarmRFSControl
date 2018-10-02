function contour_map(x0,xdes,dist_name)
%
% Function to plot Gaussian Mixtures of MPC via RFS formulation
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
num_density=length(x0);

%Load Data
load(strcat('test_case_',dist_name,'_',int2str(num_density),int2str(num_density_des),'.mat'));

N=400; %number of samples
ts=51; %number of time steps

%Create a mesh grid [-5, 5] by [-5, 5]
[X,Y]=meshgrid(linspace(-5, 5,100),linspace(-5, 5,100));
k_count=1;

%Initialize Statistics for GM
mu=zeros(4,2);
sigma=zeros(2,2,4);
for i=1:4
    mu(i,:) = x_sim2(4*i-3:4*i-2,1);
    sigma(:,:,i) = [0.5 0;0 0.5];
end

%Create GM Distribution
p = ones(1,4)/2;
obj = gmdistribution(mu,sigma,p);
points=random(obj,N);
points_sim=zeros(4,ts,N);

%Initial Conditions for each sample
for k=1:N %Initial Conditions for each sample
   points_sim(:,1,k)=[points(k,:)';0;0];
end

%Simulate Samples
for j=1:ts-1
    for i=1:4
        mu(i,:) = x_sim2(4*i-3:4*i-2,j);
        sigma(:,:,i) = [0.5 0;0 0.5];
    end
    
    p = ones(1,4)/2;
    obj = gmdistribution(mu,sigma,p);
    Fun=@(x,y)pdf(obj,[x y]);
    PDFj=reshape(Fun(reshape(X,100*100,1),reshape(Y,100*100,1)),100,100);
    %Loop through samples
    for k=1:N
        %Determine which samples go to which density
        if points_sim(1,1,k)>0 && points_sim(2,1,k)>0
            points_sim(:,1,k)=[points_sim(1:2,1,k);x_sim2(3:4,1)];
            points_sim(:,j+1,k)=A(1:4,1:4)*points_sim(:,j,k)+B(1:4,1:4)*u_sim(1:4,j);
        elseif points_sim(1,1,k)<0 && points_sim(2,1,k)>0
            points_sim(:,1,k)=[points_sim(1:2,1,k);x_sim2(7:8,1)];
            points_sim(:,j+1,k)=A(5:8,5:8)*points_sim(:,j,k)+B(5:8,5:8)*u_sim(5:8,j);
        elseif points_sim(1,1,k)<0 && points_sim(2,1,k)<0
            points_sim(:,1,k)=[points_sim(1:2,1,k);x_sim2(11:12,1)];
            points_sim(:,j+1,k)=A(9:12,9:12)*points_sim(:,j,k)+B(9:12,9:12)*u_sim(9:12,j);
        elseif points_sim(1,1,k)>0 && points_sim(2,1,k)<0
            points_sim(:,1,k)=[points_sim(1:2,1,k);x_sim2(15:16,1)];
            points_sim(:,j+1,k)=A(13:16,13:16)*points_sim(:,j,k)+B(13:16,13:16)*u_sim(13:16,j);
        end
    end
    
    %Save snapshots at certain time-steps
    if j==1
        PDFj1=PDFj; %save for time=0
    elseif j==7
        PDFj7=PDFj;
    elseif j==14
        PDFj14=PDFj;
    elseif j==50
        PDFj50=PDFj;
    end    
     k_count=k_count+1;
end

%Plot Figures
if num_density_des==4
    figure
    subplot(2,2,1);
    contour(X,Y,PDFj1,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    title('time = 0.00 s')
    grid on;
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(xdes(1,4),xdes(2,4),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,1,:)),squeeze(points_sim(2,1,:)),'.r');
    hold off;
    
    subplot(2,2,2);
    contour(X,Y,PDFj7,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    title('time = 0.05 s')
    grid on;
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(xdes(1,4),xdes(2,4),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,7,:)),squeeze(points_sim(2,7,:)),'.r')
    hold off;
    
    subplot(2,2,3);
    contour(X,Y,PDFj14,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    grid on;
    title('time = 0.10 s')
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(xdes(1,4),xdes(2,4),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,14,:)),squeeze(points_sim(2,14,:)),'.r');
    hold off;
    
    subplot(2,2,4);
    contour(X,Y,PDFj50,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    grid on;
    title('time = 0.40 s')
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(xdes(1,4),xdes(2,4),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,50,:)),squeeze(points_sim(2,50,:)),'.r');
    hold off;
    
    hp4 = get(subplot(2,2,4),'Position');
    colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.03  hp4(2)+hp4(3)*2.1-.03])
end
if num_density_des==3
    figure
    subplot(2,2,1);
    contour(X,Y,PDFj1,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    title('time = 0.00 s')
    grid on;
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,1,:)),squeeze(points_sim(2,1,:)),'.r');
    hold off;
    
    subplot(2,2,2);
    contour(X,Y,PDFj7,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    title('time = 0.05 s')
    grid on;
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,7,:)),squeeze(points_sim(2,7,:)),'.r');
    hold off;
    
    subplot(2,2,3);
    contour(X,Y,PDFj14,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    grid on;
    title('time = 0.10 s')
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,14,:)),squeeze(points_sim(2,14,:)),'.r');
    hold off;
    
    subplot(2,2,4);
    contour(X,Y,PDFj50,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    grid on;
    title('time = 0.40 s')
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,50,:)),squeeze(points_sim(2,50,:)),'.r');
    hold off;
    
    hp4 = get(subplot(2,2,4),'Position');
    colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.03  hp4(2)+hp4(3)*2.1-.032])
end
if num_density_des==5
    figure
    subplot(2,2,1);
    contour(X,Y,PDFj1,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    title('time = 0.00 s')
    grid on;
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(xdes(1,4),xdes(2,4),'xk','LineWidth',2)
    plot(xdes(1,5),xdes(2,5),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,1,:)),squeeze(points_sim(2,1,:)),'.r');
    hold off;
    
    subplot(2,2,2);
    contour(X,Y,PDFj7,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    title('time = 0.05 s')
    grid on;
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(xdes(1,4),xdes(2,4),'xk','LineWidth',2)
    plot(xdes(1,5),xdes(2,5),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,7,:)),squeeze(points_sim(2,7,:)),'.r');
    hold off;
    
    subplot(2,2,3);
    contour(X,Y,PDFj14,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    grid on;
    title('time = 0.10 s')
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(xdes(1,4),xdes(2,4),'xk','LineWidth',2)
    plot(xdes(1,5),xdes(2,5),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,14,:)),squeeze(points_sim(2,14,:)),'.r');
    hold off;
    
    subplot(2,2,4);
    contour(X,Y,PDFj50,linspace(0.006,0.218,100),'linewidth',2)
    hold on;
    set(gcf,'color',[1 1 1]);
    xlabel('x')
    ylabel('y')
    grid on;
    title('time = 0.40 s')
    plot(xdes(1,1),xdes(2,1),'xk','LineWidth',2)
    plot(xdes(1,2),xdes(2,2),'xk','LineWidth',2)
    plot(xdes(1,3),xdes(2,3),'xk','LineWidth',2)
    plot(xdes(1,4),xdes(2,4),'xk','LineWidth',2)
    plot(xdes(1,5),xdes(2,5),'xk','LineWidth',2)
    plot(squeeze(points_sim(1,50,:)),squeeze(points_sim(2,50,:)),'.r');
    hold off;
    hp4 = get(subplot(2,2,4),'Position');
    colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.03  hp4(2)+hp4(3)*2.1-.03]);
    
end