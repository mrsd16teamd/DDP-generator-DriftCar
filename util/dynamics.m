function [dx] = dynamics(x, u, params)

% Drifting dynamics developed based on 
% Dynamics And Control Of Drifting In Automobiles, Hindiyeh 2013

% The parameters used in this example refer to a 1/10 scale RC car
% The model is expected to work for full scale vehicle as well

% ----------------------------------------
% --------------Model Params--------------
% ----------------------------------------
m = 2.35;           % mass (kg)
g = 9.81;
L = 0.257;          % wheelbase (m)
b = 0.14328;        % CoG to rear axle
a = L-b;            % CoG to front axle
G_front = m*g*b/L;   % calculated load or specify front rear load directly
G_rear = m*g*a/L;

C_x = 116;   %330       % longitudinal stiffness
C_alpha = 197; %300      % lateral stiffness
Iz = 0.025; %0.02065948883 % rotational inertia
mu = 1.31; %5.2/G_rear   
mu_spin = 0.55; %4.3/G_rear

if exist('params','var')
    % params = [C_alpha, C_x, Iz, mu, mu_spin]
    C_alpha = params.c_a;
    C_x = params.c_x;
    Iz = params.Iz;
    mu = params.mu;
    mu_spin = params.mu_s;
end

% ----------------------------------------
% ------States/Inputs Interpretation------
% ----------------------------------------
pos_x = x(1,:);
pos_y = x(2,:);
pos_phi = x(3,:);

Ux = x(4,:);
Uy = x(5,:);
r = x(6,:);

Ux_cmd = u(1,:);
delta = u(2,:);

% ----------------------------------------
% --------------Tire Dynamics-------------
% ----------------------------------------
% lateral slip angle alpha
% if Ux == 0 && Uy == 0   % vehicle is still no slip
%     alpha_F = 0;
%     alpha_R = 0;
% elseif Ux == 0      % perfect side slip
%     alpha_F = pi/2*sign(Uy)-delta;
%     alpha_R = pi/2*sign(Uy);
% elseif Ux < 0    % rare ken block situations
%     alpha_F = (sign(Uy)*pi)-atan((Uy+a*r)/abs(Ux))-delta;
%     alpha_R = (sign(Uy)*pi)-atan((Uy-b*r)/abs(Ux));
% else                % normal situation
%     alpha_F = atan((Uy+a*r)/abs(Ux))-delta;
%     alpha_R = atan((Uy-b*r)/abs(Ux));
% end
alpha_F = atan((Uy+a*r)/(abs(Ux)+1e-3))-(Ux/(abs(Ux)+1e-3))*delta;

alpha_R = atan((Uy-b*r)/(abs(Ux)+1e-3));

% safety that keep alpha in valid range
% alpha_F = wrapToPi(alpha_F);
% alpha_R = wrapToPi(alpha_R);
% disp([alpha_F, alpha_R]);
% disp([Ux, Ux_cmd, mu, mu_spin, G_rear, C_x, C_alpha]);

[Fxf,Fyf] = tire_dyn(Ux, Ux, mu, mu_spin, G_front, C_x, C_alpha, alpha_F);
[Fxr,Fyr] = tire_dyn(Ux, Ux_cmd, mu, mu_spin, G_rear, C_x, C_alpha, alpha_R);

% disp([Fxf, Fyf, Fxr, Fyr])
% ----------------------------------------
% ------------Vehicle Dynamics------------
% ----------------------------------------
% ddx
r_dot = (a*Fyf*cos(delta)-b*Fyr)/Iz;
Ux_dot = (Fxr-Fyf*sin(delta))/m+r*Uy;
Uy_dot = (Fyf*cos(delta)+Fyr)/m-r*Ux;

% translate dx to terrain frame
U = sqrt(Ux^2+Uy^2+1e-6);

if Ux < 0
    beta = sign(Uy)*pi-atan(Uy/(abs(Ux)+1e-3));
else
    beta = atan(Uy/abs(Ux));
end
% beta = wrapToPi(beta);

Ux_terrain = U*cos(beta+pos_phi);
Uy_terrain = U*sin(beta+pos_phi);
dx = [Ux_terrain;Uy_terrain;r;Ux_dot;Uy_dot;r_dot];
end
