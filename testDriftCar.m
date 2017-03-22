%% parameters
p.h= 0.02;

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
p.cdu= 1e-1*[1 15];

p.cf= [ 10 10 1 .1 .1 .1];
p.pf= [ .01 .01 .1 .1 .1 .1];

p.cx  = 1e-2*[5  5 4];          % running cost coefficients 
p.cdx = 1e-1*[5 0.5 0.2];
p.px  = [.01 .01 .1];   % smoothness scales for running cost

p.cdrift = -0.001;

p.limThr= [-1 4];
p.limSteer= [-0.68  0.76];

p.xDes = [5 0 0 3 0 0]; % Moose Test
% p.xDes = [1.2 0.5 pi 0 0 0]; % Parallel Park


p.k_pos = 0.1;
p.k_vel = 0;
p.d_thres = 0.3;
p.Obs = [1 0];

%% initial conditions
T= 50;              % horizon
t= (1:501)*p.h;
x0= [0;0;0;3;0;0;3;0;0;0];   % initial state - Moose Test
% x0= [0;0;0;3.5;0;0;3.5;0;0;0];   % initial state - Parallel Park


x00 = x0;
u0(1,:) = 0.25*randn(1,T) + 3; % commanded speed
u0(2,:) = 0.1*randn(1,T) + 0.1; % steering

Op.max_iter= 100;

figure(1)
init_plot(x0,p.xDes,p.Obs);

X = [];
U = [];
while pdist([x0(1:2)';p.xDes(1:2)]) > 0.15
tic
[success, x, u, cost]= iLQGDriftCarFlat(x0, u0, p, Op);
X = [X x(:,1:T/2)];
U = [U u(:,1:T/2)];

p.cf = 5*p.cf;
p.cdu = 2*p.cdu;
x0 = x(:,T/2);
u0 = [u(:,T/2+1:end),.1*randn(2,T/2)];
toc
end

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


