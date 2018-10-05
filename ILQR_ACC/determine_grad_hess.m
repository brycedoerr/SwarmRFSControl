function determine_grad_hess(num_density,num_density_des,size_state,size_ctrl,fh,dist_name)
% Function to determine gradient and hessian of cost function. Makes Fcn.
%
% Synopsis:
%     determine_grad_hess(num_density,num_density_des,size_state,size_ctrl,fh,dist_name)
%
% Input:
%     num_density        =   length current intensities 
%     num_density_des    =   length of target intensities
%     size_state         =   length of State 
%     size_ctrl          =   length of control sequence
%     fh                 =   function handle for cost
%     dist_name          =   label Distance Name is function generated
%
% By: Bryce Doerr -- Aug. 2018

x = cell(size_state,num_density);
for i = 1:size_state
    for j = 1:num_density
        x{i,j} = sprintf('x%d%d',i,j);
    end
end
x = x(:); % now x is a 16-by-1 vector
x = sym(x, 'real');

size_state_des=size_state;
xdes = cell(size_state_des,num_density_des);
for i = 1:size_state_des
    for j = 1:num_density_des
        xdes{i,j} = sprintf('xdes%d%d',i,j);
    end
end
xdes = xdes(:); % now x is a 16-by-1 vector
xdes = sym(xdes, 'real');

R = cell(size_ctrl,num_density);
for i = 1:size_ctrl*num_density
    for j = 1:num_density*size_ctrl
        R{i,j} = sprintf('R%d%d',i,j);
    end
end
%R = R(:); % now x is a 16-by-1 vector
R = sym(R, 'real');


u = cell(size_ctrl,num_density);
for i = 1:size_ctrl
    for j = 1:num_density
        u{i,j} = sprintf('u%d%d',i,j);
    end
end
u = u(:); % now x is a 16-by-1 vector
u = sym(u, 'real');

cost_x=gradient(fh(u,reshape(x,size_state,num_density),reshape(xdes,size_state_des,num_density_des),R),x);
cost_xx=jacobian(cost_x,x);
cost_u=gradient(fh(u,reshape(x,size_state,num_density),reshape(xdes,size_state_des,num_density_des),R),u);
cost_uu=jacobian(cost_u,u);
cost_ux=jacobian(cost_u,x);

matlabFunction(cost_x,cost_xx,cost_u,cost_uu,cost_ux,'file',strcat('grad_hess_symbolic_',dist_name,'_',int2str(num_density),int2str(num_density_des)),'vars',{x,xdes,u,R});


% % %old
% cost_x=gradient(cost_nl_l2_wo_control(reshape(x,size_state,num_density),reshape(xdes,size_state_des,num_density_des)),x);
% cost_xx=jacobian(cost_x,x);
% matlabFunction(cost_x,cost_xx,'file','grad_hess_symbolic','vars',{x,xdes});

fprintf('************Gradients and Hessians Calculated!**************\n');
