function ilqr_script(x0,xdes,t,dt,fh_cost,grad_hess_symbolic,dist_name)
%
% Function to simulate ILQR via RFS formulation
%
% Synopsis:
%     ilqr_script(x0,xdes,t,dt,fh_cost,grad_hess_symbolic,dist_name)
%
% Input:
%     x0          =   states for each density (# states x # densities)
%     xdes        =   target states (# target states x # target densities)
%     t           = time sequence (column vector)
%     dt          = time-step
%     fh_cost     = function handle for cost
%     grad_hess_symbolic     = function handle for gradient & hessian cost
%     dist_name   =   string for labelling distances

%
%
% By: Bryce Doerr -- Aug. 2018
% ILQR Based on Tassa's MATLAB Algorithm:
% https://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization
%% Initialize Dynamics with number of densities
tic
Z=[0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0;]; %build the A matrix
z=zeros(4,4); %zero matrix to build the A matrix
A=[Z z z z; z Z z z; z z Z z; z z z Z];
W=[0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0;];%build the B matrix
z1=zeros(4,4);%zero matrix to build the B matrix
B=[W z1 z1 z1; z1 W z1 z1; z1 z1 W z1; z1 z1 z1 W;];
V=[1 0 0 0; 0 1 0 0]; %Build C matrix
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
x_sim=zeros(length(A),length(t));

ushape=zeros(1,16);
u=repmat(ushape,length(t),1)';

x_sim(:,1)=reshape(x0,[],1);

%Initialize Quadratized Cost Function
l_x=zeros(length(A),length(t));
l_xx=zeros(length(A),length(A),length(t));
l_u=zeros(length(B),length(t));
l_uu=zeros(length(B),length(B),length(t));
l_ux=zeros(length(B),length(A),length(t));
V_x=zeros(length(A),length(t));
V_xx=zeros(length(A),length(A),length(t));

Q_x=zeros(length(A),length(t));
Q_u=zeros(length(B),length(t));
Q_xx=zeros(length(A),length(A),length(t));
Q_ux=zeros(length(B),length(A),length(t));
Q_uu=zeros(length(B),length(B),length(t));
Q_ux_reg=zeros(length(B),length(A),length(t));
Q_uu_reg=zeros(length(B),length(B),length(t));

%Initialize Control
k_control=zeros(length(B),length(t));
K=zeros(length(B),length(A),length(t));

%Converging Tolerance
tol=1e-3;
tolGrad=1e-3;

%Regularization
delta_0=2;%factor
delta=1;
lamb_max=1e10;
lamb_min=1e-6;
lamb=1;

%Current Response + Current cost
u_bar=zeros(size(u,1),size(u,2));
x_bar=zeros(length(A),length(t));
l_bar=zeros(1,length(t));
x_bar(:,1)=x_sim(:,1);

for k=1:length(t)-1
    x_bar(:,k+1)=A*x_bar(:,k)+B*u_bar(:,k);
    l_bar(:,k)=fh_cost(u_bar(:,k),reshape(x_bar(:,k),4,4),xdes,R);
end
l_bar(:,end)=fh_cost(u_bar(:,end),reshape(x_bar(:,end),4,4),xdes,R);
cost_bar=sum(l_bar(1:end-1));

%Line Search parameters
alpha=10.^linspace(0,-3,11);
count=0;
zMin=0;
n_iter=100;

%Start ILQR Iteration
for iter=1:n_iter
    count=count+1;
    
    %Differentiate cost along new trajectory
    for k=1:length(t)-1
        [l_x(:,k),l_xx(:,:,k),l_u(:,k),l_uu(:,:,k),l_ux(:,:,k)]=grad_hess_symbolic(x_bar(:,k),reshape(xdes,size(xdes,1)*size(xdes,2),1),u_bar(:,k),R);
    end
    [l_x(:,end),l_xx(:,:,end),l_u(:,end),l_uu(:,:,end),l_ux(:,:,end)]=grad_hess_symbolic(x_bar(:,end),reshape(xdes,size(xdes,1)*size(xdes,2),1),u_bar(:,end),R);

    %Backward pass
    backPassDone=0;
    while~backPassDone
        dV=[0 0];
        diverge =0;
        V_x(:,end)=l_x(:,end);
        V_xx(:,:,end)=l_xx(:,:,end);
            
        for k=length(t)-1:-1:1
            Q_x(:,k)=l_x(:,k)+A'*V_x(:,k+1);
            Q_u(:,k)=l_u(:,k)+B'*V_x(:,k+1);
            Q_xx(:,:,k)=l_xx(:,:,k)+A'*V_xx(:,:,k+1)*A;
            Q_ux(:,:,k)=l_ux(:,:,k)+B'*V_xx(:,:,k+1)*A;
            Q_uu(:,:,k)=l_uu(:,:,k)+B'*V_xx(:,:,k+1)*B;

            Q_ux_reg(:,:,k)=l_ux(:,:,k)+B'*(V_xx(:,:,k+1)+lamb*eye(length(A)))*A;
            Q_uu_reg(:,:,k)=l_uu(:,:,k)+B'*(V_xx(:,:,k+1)+lamb*eye(length(A)))*B; 
            
            [R_chol,d]=chol(Q_uu_reg(:,:,k));
            if d~=0
                diverge=k;
                break;
            end
            kK=-R_chol\(R_chol'\[Q_u(:,k) Q_ux_reg(:,:,k)]);
            k_control(:,k)=kK(:,1);
            K(:,:,k)=kK(:,2:length(A)+1);

            dV= dV + [k_control(:,k)'*Q_u(:,k)  .5*k_control(:,k)'*Q_uu(:,:,k)*k_control(:,k)];
            V_x(:,k)=Q_x(:,k)+K(:,:,k)'*Q_uu(:,:,k)*k_control(:,k)+K(:,:,k)'*Q_u(:,k)+Q_ux(:,:,k)'*k_control(:,k);
            V_xx(:,:,k)=Q_xx(:,:,k)+K(:,:,k)'*Q_uu(:,:,k)*K(:,:,k)+K(:,:,k)'*Q_ux(:,:,k)+Q_ux(:,:,k)'*K(:,:,k);
            V_xx(:,:,k)=.5*(V_xx(:,:,k)+V_xx(:,:,k)');
        end
        
        if diverge
            fprintf('Cholesky failed at timestep %d.\n',diverge);
            %increase regularization term
            delta=max(1*delta_0,delta*delta_0);
            lamb=max(lamb_min,lamb*delta);
            if lamb>lamb_max
                fprintf( 'Exceeds regularization limit at iteration %d\n',count)
                break;
            end
            continue
        end
    backPassDone=1;
    end
    
    g_norm=mean(max(abs(k_control) ./ (abs(u_bar)+1),[],1));
    if g_norm < tolGrad && lamb < 1e-5
        fprintf('\nSUCCESS: gradient norm < tolGrad\n');
        break;
    end
    
    %Forward Pass
    fwdPassDone = 0;
    if backPassDone
        u_new=zeros(size(u,1),size(u,2));
        x_new=zeros(length(A),length(t));
        x_new(:,1)=x_sim(:,1);
        l_new=zeros(1,length(t));
        
        for alp=alpha
            for k=1:length(t)-1
                u_new(:,k)=u_bar(:,k)+alp*(k_control(:,k))+K(:,:,k)*(x_new(:,k)-x_bar(:,k));
                x_new(:,k+1)=A*x_new(:,k)+B*u_new(:,k);
                l_new(:,k)=fh_cost(u_new(:,k),reshape(x_new(:,k),4,4),xdes,R);
            end
            l_new(:,end)=fh_cost(u_new(:,end),reshape(x_new(:,end),4,4),xdes,R);
            
            cost_new=sum(l_new(1:end-1));
            dcost=cost_bar-cost_new;
            expected = -alp*(dV(1) + alp*dV(2));
            if expected > 0
                z = dcost/expected;
            else
                z = sign(dcost);
                warning('non-positive expected reduction: should not occur');
            end
            if (z > zMin)
                fwdPassDone = 1;
                break;
            end
        end
    end
    if fwdPassDone
        %decrease regularization term
        delta=min(1/delta_0,delta/delta_0);
        lamb=lamb*delta;
        if lamb<=lamb_min
           lamb=0; 
        end
        u_bar=u_new;
        x_bar=x_new;
        cost_bar=cost_new;
        
        if dcost<tol
            formatSpec = 'Converged at iteration = %d\n Cost = %f\n logLambda = %f\n alpha = %f\n';
            fprintf(formatSpec,count,cost_new,log10(lamb),alp);
            break;
        end
        formatSpec = 'Iteration = %d\n Cost = %f\n logLambda = %f\n alpha = %f\n';
        fprintf(formatSpec,count,cost_new,log10(lamb),alp);    


    else
        %increase regularization term
        delta=max(1*delta_0,delta*delta_0);
        lamb=max(lamb_min,lamb*delta);
        if lamb>lamb_max
            fprintf( 'Exceeds regularization limit at iteration %d\n',count)
            break;
        end
        fprintf('Reject the control perturbation. Increase the regularization term.\n');
    end
end
toc

%% Plot Results
figure
subplot(4,1,1);plot(t,x_bar(1,:),'-b','LineWidth',2);hold on;plot(t,1*x_bar(1,:).^0,'--r')
plot(t,x_bar(2,:),'.g','LineWidth',2);plot(t,1*x_bar(1,:).^0,'--r')
xlabel('time (s)');ylabel('Pos. 1');grid on;

subplot(4,1,2);plot(t,x_bar(5,:),'-b','LineWidth',2);hold on;plot(t,-1*x_bar(1,:).^0,'--r')
plot(t,x_bar(6,:),'.g','LineWidth',2);plot(t,1*x_bar(1,:).^0,'--r')
xlabel('time (s)');ylabel('Pos. 2');grid on;

subplot(4,1,3);plot(t,x_bar(9,:),'-b','LineWidth',2);hold on;plot(t,-1*x_bar(1,:).^0,'--r')
plot(t,x_bar(10,:),'.g','LineWidth',2);plot(t,-1*x_bar(1,:).^0,'--r')
xlabel('time (s)');ylabel('Pos. 3');grid on;

subplot(4,1,4);
plot(t,x_bar(13,:),'-b','LineWidth',2);
hold on;
plot(t,1*x_bar(1,:).^0,'--r');
plot(t,x_bar(14,:),'.g','LineWidth',2);
plot(t,-1*x_bar(1,:).^0,'--r');
xlabel('time (s)');ylabel('Pos. 4');grid on;

figure
subplot(4,1,1);plot(t,x_bar(3,:),'-b','LineWidth',2);hold on;plot(t,1*x_bar(3,:).^0,'--r')
plot(t,x_bar(4,:),'.g','LineWidth',2);plot(t,1*x_bar(4,:).^0,'--r')
xlabel('time (s)');ylabel('Vel. 1');grid on;

subplot(4,1,2);plot(t,x_bar(7,:),'-b','LineWidth',2);hold on;plot(t,-1*x_bar(7,:).^0,'--r')
plot(t,x_bar(8,:),'.g','LineWidth',2);plot(t,1*x_bar(8,:).^0,'--r')
xlabel('time (s)');ylabel('Vel. 2');grid on;

subplot(4,1,3);plot(t,x_bar(11,:),'-b','LineWidth',2);hold on;plot(t,-1*x_bar(11,:).^0,'--r')
plot(t,x_bar(12,:),'.g','LineWidth',2);plot(t,-1*x_bar(12,:).^0,'--r')
xlabel('time (s)');ylabel('Vel. 3');grid on;

subplot(4,1,4);
plot(t,x_bar(15,:),'-b','LineWidth',2);
hold on;
plot(t,1*x_bar(15,:).^0,'--r');
plot(t,x_bar(16,:),'.g','LineWidth',2);
plot(t,-1*x_bar(16,:).^0,'--r');
xlabel('time (s)');ylabel('Vel. 4');grid on;


%% Save Results
num_density_des=size(xdes,2);
num_density=size(x0,2);
save(strcat('test_case_',dist_name,'_',int2str(num_density),int2str(num_density_des),'.mat'));