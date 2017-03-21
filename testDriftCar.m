%% parameters
p.h= 0.05;

p.m = 2.35;           % mass (kg)
g = 9.81;
L = 0.257;          % wheelbase (m)
p.b = 0.14328;        % CoG to rear axle
p.a = L-p.b;            % CoG to front axle
p.G_f = p.m*g*p.b/L;   % calculated load or specify front rear load directly
p.G_r = p.m*g*p.a/L;

p.c_x = 45;          % longitude stiffness
p.c_a = 50;      % laternal stiffness
p.Iz = 0.045; % roatation inertia
p.mu = 0.75;   
p.mu_s = 0.2;

p.cu= 1e-2*[1 .1];
p.cdu= 1e-1*[.01 1];

p.cf= [ 10 10 1 .1 .1 .1];
p.pf= [ .01 .01 .1 .1 .1 .1];

p.cx  = 1e-1*[5  5 4];          % running cost coefficients 
p.cdx = 1e-2*[0.5 0.5 0.2];
p.px  = [.01 .01 .1];   % smoothness scales for running cost

p.cdrift = -0.001;

p.limThr= [-1 4];
p.limSteer= [-0.68  0.76];

p.xDes = [3 0 0 0 0 0];

p.k_pos = 0.5;
p.k_vel = 0.1;
p.d_thres = 0.3;
p.Obs = [1 0];

%% initial conditions
T= 50;              % horizon
t= (1:501)*p.h;
x0= [0;0;0;3;0;0;0;0;0;0];   % initial state
u0(1,:) = 0.25*randn(1,T) +2; % commanded speed
u0(2,:) = 0.1*randn(1,T); % steering

Op.max_iter= 100;

tic
[success, x, u, cost]= iLQGDriftCar(x0, u0, p, Op);
toc

figure
init_plot(x0,p.xDes,p.Obs);
car_plot(x,u);
