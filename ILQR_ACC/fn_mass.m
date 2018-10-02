function dm=fn_mass(t,m,u)
%
% Equations of motion for a mass system
%
% Synopsis:
%     fn_mass(t,m,u)
%
% Input:
%     t       =   time
%     m       =   mass of agent 
%     u       =   control input (acceleration)
%
% Output:
%     dm      = mass at next time-step
%
% By: Bryce Doerr -- Aug. 2018

    %Mass Equations of Motion
    % u is m/s^2
    g=9.8; % m/s^2
    Isp=300; % sec
    dm=-m*u/(g*Isp);% units of dm is kg/s


end