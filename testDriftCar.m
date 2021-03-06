%% parameters
p.dt= 0.02;

p.m = 2.35;           % mass (kg)
g = 9.81;
L = 0.257;          % wheelbase (m)
p.b = 0.14328;        % CoG to rear axle
p.a = L-p.b;            % CoG to front axle
p.G_f = p.m*g*p.b/L;   % calculated load or specify front rear load directly
p.G_r = p.m*g*p.a/L;

p.c_x = 116;          % longitude stiffness
p.c_a = 197;      % laternal stiffness
p.Iz = 0.025; % roatation inertia
p.mu = 1.31;   
p.mu_s = 0.55;

p.cu= 1e-1*[1 1];
p.cdu= 1e-2*[4 6];

p.cf= [ 15 15 5 10 .1 .1];
p.pf= [ .01 .01 .1 .1 .1 .1];

p.cx  = 1e-2*[15 5 3];          % running cost coefficients 
p.cdx = 1e-3*[1 5 2];
p.px  = [.01 .01 .1];   % smoothness scales for running cost

p.cdrift = -0.01;

p.limThr= [0 4];
p.limSteer= [-0.77  0.77];

p.xDes = [5 0 0 0 0 0]; % Moose Test

p.lane_center = 0.34;
p.lane_thres = 0.30;
p.croad = 1;
% p.xDes = [1.2 0.5 pi 0 0 0]; % Parallel Park


p.k_pos = 0.5;
p.k_vel = 0;
p.d_thres = 0.3;
p.Obs = [1 0];

%% initial conditions
T= 101;              % horizon

x0= [0;0;0;3;0;0;3;0];   % initial state - Moose Test
% x0= [0;0;0;3.5;0;0;3.5;0;0;0];   % initial state - Parallel Park


x00 = x0;

u0(1,:) = 0.25*randn(1,T) + 3; % commanded speed
u0(2,:) = 0.1*randn(1,T) + 0.2; % steering

Op.max_iter= 5000;

figure(1)
init_plot(x0,p.xDes,p.Obs);

X = [];
U = [];
% while pdist([x0(1:2)';p.xDes(1:2)]) > 0.15
tic
[success, x, u, cost]= iLQGDriftCar(x0, u0, p, Op);
X = [X x(:,:)];
U = [U u(:,:)];

% p.cf = 5*p.cf;
% p.cdu = 2*p.cdu;
% x0 = x(:,T/2);
% u0 = [u(:,T/2+1:end),.1*randn(2,T/2)];
toc
% end

figure(1)
car_plot(X,U);
rerun(x00,U,p);

figure(2)
subplot(2,1,1)
title('Throttle')
plot(U(1,:))
subplot(2,1,2)
title('Steering')
plot(U(2,:))


