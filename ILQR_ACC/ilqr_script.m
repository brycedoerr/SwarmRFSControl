function ilqr_script(x0,xdes,t,dt,u_ctrl,fh_cost,grad_hess_symbolic,dist_name)
%
% Function to simulate ILQR via RFS formulation
%
% Synopsis:
%     ilqr_script(x0,xdes,t,dt,u_ctrl,fh_cost,grad_hess_symbolic,dist_name)
%
% Input:
%     x0          =   states for each density (# states x # densities)
%     xdes        =   target states (# target states x # target densities)
%     t           = time sequence (column vector)
%     dt          = time-step
%     u_ctrl      = Initial control input
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
n=.00113;%ang freq
Z= [0     0 0 1 0   0;...
    0     0 0 0 1   0;...
    0     0 0 0 0   1;...
    3*n^2 0 0 0 2*n  0;...
    0    0 0 -2*n 0 0;...
    0 0 -n^2 0 0 0;]; %build the A matrix
A=kron(eye(size(x0,2)),Z);
W=[0 0 0; 0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1;];%build the B matrix. Control input u is m/s^2
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

%% Simulation
R=8*eye(3*size(x0,2));
x_sim=zeros(length(A),length(t));

ushape=zeros(1,3*size(x0,2));
u=repmat(ushape,length(t),1)';

x_sim(:,1)=reshape(x0,[],1);

%Initialize Quadratized Cost Function
l_x=zeros(length(A),length(t));
l_xx=zeros(length(A),length(A),length(t));
l_u=zeros(size(B,2),length(t));
l_uu=zeros(size(B,2),size(B,2),length(t));
l_ux=zeros(size(B,2),length(A),length(t));

%Initialize Control
k_control=zeros(size(B,2),length(t));
K=zeros(size(B,2),length(A),length(t));

V_x=zeros(length(A),length(t));
V_xx=zeros(length(A),length(A),length(t));

%Converging Tolerance
tol=1e-3;
tolGrad=1e-3;

%Regularization
delta_0=1.6;%factor
delta=1;
lamb_max=1e10;
lamb_min=1e-6;
lamb=1;

%Current Response + Current cost
u_bar=u_ctrl;
x_bar=zeros(length(A),length(t));
l_bar=zeros(1,length(t));
x_bar(:,1)=x_sim(:,1);

for k=1:length(t)-1
    x_bar(:,k+1)=A*x_bar(:,k)+B*u_bar(:,k);
    l_bar(:,k)=fh_cost(u_bar(:,k),reshape(x_bar(:,k),size(x0,1),size(x0,2)),xdes(:,:,k),R);
end
l_bar(:,end)=fh_cost(u_bar(:,end),reshape(x_bar(:,end),size(x0,1),size(x0,2)),xdes(:,:,end),R);
cost_bar=sum(l_bar(1:end));

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
        [l_x(:,k),l_xx(:,:,k),l_u(:,k),l_uu(:,:,k),l_ux(:,:,k)]=grad_hess_symbolic(reshape(x_bar(:,k),size(x0,1),size(x0,2)),xdes(:,:,k),u_bar(:,k),R);
    end
    [l_x(:,end),l_xx(:,:,end),l_u(:,end),l_uu(:,:,end),l_ux(:,:,end)]=grad_hess_symbolic(reshape(x_bar(:,end),size(x0,1),size(x0,2)),xdes(:,:,end),u_bar(:,end),R);
    
    %Backward pass 
    backPassDone=0;
    while~backPassDone
        dV=[0 0];
        V_x(:,end)=l_x(:,end);
        V_xx(:,:,end)=l_xx(:,:,end);
        diverge =0;
        
        for k=length(t)-1:-1:1
            Q_x=l_x(:,k)+A'*V_x(:,k+1);
            Q_u=l_u(:,k)+B'*V_x(:,k+1);
            Q_xx=l_xx(:,:,k)+A'*V_xx(:,:,k+1)*A;
            Q_ux=l_ux(:,:,k)+B'*V_xx(:,:,k+1)*A;
            Q_uu=l_uu(:,:,k)+B'*V_xx(:,:,k+1)*B;
            
            Q_ux_reg=l_ux(:,:,k)+B'*(V_xx(:,:,k+1)+lamb*eye(length(A)))*A;
            Q_uu_reg=l_uu(:,:,k)+B'*(V_xx(:,:,k+1)+lamb*eye(length(A)))*B;%reg n
            
            [R_chol,d]=chol(Q_uu_reg);
            if d~=0
                diverge=k;
                break;
            end
            kK=-R_chol\(R_chol'\[Q_u Q_ux_reg]);
            k_control(:,k)=kK(:,1);
            K(:,:,k)=kK(:,2:length(A)+1);

             dV= dV + [k_control(:,k)'*Q_u  .5*k_control(:,k)'*Q_uu*k_control(:,k)];
             V_x(:,k)=Q_x+K(:,:,k)'*Q_uu*k_control(:,k)+K(:,:,k)'*Q_u+Q_ux'*k_control(:,k);
             V_xx(:,:,k)=Q_xx+K(:,:,k)'*Q_uu*K(:,:,k)+K(:,:,k)'*Q_ux+Q_ux'*K(:,:,k);
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
                l_new(:,k)=fh_cost(u_new(:,k),reshape(x_new(:,k),size(x0,1),size(x0,2)),xdes(:,:,k),R);
            end
            l_new(:,end)=fh_cost(u_new(:,end),reshape(x_new(:,end),size(x0,1),size(x0,2)),xdes(:,:,k),R);
            
            cost_new=sum(l_new(1:end));
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
        lamb=lamb*delta*(lamb>lamb_min);

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

subplot(4,1,2);plot(t,x_bar(7,:),'-b','LineWidth',2);hold on;plot(t,-1*x_bar(1,:).^0,'--r')
plot(t,x_bar(8,:),'.g','LineWidth',2);plot(t,1*x_bar(1,:).^0,'--r')
xlabel('time (s)');ylabel('Pos. 2');grid on;

subplot(4,1,3);plot(t,x_bar(13,:),'-b','LineWidth',2);hold on;plot(t,-1*x_bar(1,:).^0,'--r')
plot(t,x_bar(14,:),'.g','LineWidth',2);plot(t,-1*x_bar(1,:).^0,'--r')
xlabel('time (s)');ylabel('Pos. 3');grid on;

subplot(4,1,4);
plot(t,x_bar(19,:),'-b','LineWidth',2);
hold on;
plot(t,1*x_bar(1,:).^0,'--r');
plot(t,x_bar(20,:),'.g','LineWidth',2);
plot(t,-1*x_bar(1,:).^0,'--r');
xlabel('time (s)');ylabel('Pos. 4');grid on;

figure
subplot(4,1,1);plot(t,x_bar(4,:),'-b','LineWidth',2);hold on;plot(t,1*x_bar(3,:).^0,'--r')
plot(t,x_bar(5,:),'.g','LineWidth',2);plot(t,1*x_bar(4,:).^0,'--r')
xlabel('time (s)');ylabel('Vel. 1');grid on;

subplot(4,1,2);plot(t,x_bar(10,:),'-b','LineWidth',2);hold on;plot(t,-1*x_bar(7,:).^0,'--r')
plot(t,x_bar(11,:),'.g','LineWidth',2);plot(t,1*x_bar(8,:).^0,'--r')
xlabel('time (s)');ylabel('Vel. 2');grid on;

subplot(4,1,3);plot(t,x_bar(16,:),'-b','LineWidth',2);hold on;plot(t,-1*x_bar(11,:).^0,'--r')
plot(t,x_bar(17,:),'.g','LineWidth',2);plot(t,-1*x_bar(12,:).^0,'--r')
xlabel('time (s)');ylabel('Vel. 3');grid on;

subplot(4,1,4);
plot(t,x_bar(22,:),'-b','LineWidth',2);
hold on;
plot(t,1*x_bar(15,:).^0,'--r');
plot(t,x_bar(23,:),'.g','LineWidth',2);
plot(t,-1*x_bar(16,:).^0,'--r');
xlabel('time (s)');ylabel('Vel. 4');grid on;


%% Save Results
num_density_des=size(xdes,2);
num_density=size(x0,2);
save(strcat('test_case_',dist_name,'_',int2str(num_density),int2str(num_density_des),'.mat'));