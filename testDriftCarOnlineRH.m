function testDriftCarOnlineRH

%% Initialize Global MATLAB-ROS node if not already active

tx1 = rosdevice('192.168.1.118','ubuntu','ubuntu');

try
    rosinit('http://192.168.1.118:11311')
catch EXP
    if string(EXP.identifier)=='robotics:ros:node:GlobalNodeNotRunning'
        disp('ROS Master not running!')
    elseif string(EXP.identifier)=='robotics:ros:node:NoMasterConnection'
        disp('Cannot connect to ROS master! Check network')
        return
    end
end

%% Initialize Publishers, Subscribers, and Messages

obs_sub = rossubscriber('/ccs', 'std_msgs/Float32MultiArray', 'BufferSize', 10);
obs_msg = obs_sub.LatestMessage;

amcl_sub = rossubscriber('/amcl_pose_echo', 'geometry_msgs/PoseWithCovarianceStamped', 'BufferSize', 10);
amcl_msg = amcl_sub.LatestMessage;

twist_sub = rossubscriber('/cmd_vel', 'geometry_msgs/Twist', 'BufferSize', 10);
twist_msg_sub = twist_sub.LatestMessage;

odom_sub = rossubscriber('/odometry/filtered', 'nav_msgs/Odometry', 'BufferSize',10);
odom_msg = odom_sub.LatestMessage;

twist_chatpub = rospublisher('/cmd_vel','geometry_msgs/Twist');
twist_msg = rosmessage(twist_chatpub);

traj_chatpub = rospublisher('/traj_sim','nav_msgs/Odometry');
traj_msg = rosmessage(traj_chatpub);

disp('Created following subscribers:');
disp(obs_sub);

disp('Created following pubishers:');
disp(twist_chatpub);
disp(traj_chatpub);


%% Parameters

disp('Initializing iLQR parameters')

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
p.cdu= 1e-1*[10 15];

p.cf= [ 10 10 1 .1 .1 .1];
p.pf= [ .01 .01 .1 .1 .1 .1];

p.cx  = 1e-2*[5  5 4];          % running cost coefficients 
p.cdx = 1e-1*[5 0.5 0.2];
p.px  = [.01 .01 .1];   % smoothness scales for running cost

p.cdrift = -0.001;

p.limThr= [-1 4];
p.limSteer= [-0.68  0.76];

p.xDes = [3 0 0 3 0 0]; % Moose Test
% p.xDes = [1.2 0.5 pi 0 0 0]; % Parallel Park


p.k_pos = 0;
p.k_vel = 0;
p.d_thres = 0.3;
p.Obs = [1 0];

%% Initial conditions
T = 50;              % horizon
t = (1:501)*p.h;
x0= [0;0;0;3;0;0;3;0;0;0];   % initial state - Moose Test
% x0= [0;0;0;3;0;0;3;0;0;0];   % initial state - Parallel Park

x00 = x0;
% load('moose_test_sim.mat','u')
% u0 = u; % Load pre-generated trajectory as initialization
u0(1,:) = 0.25*randn(1,T) + 3; % commanded speed
u0(2,:) = 0.1*randn(1,T) + 0.1; % steering

Op.max_iter= 100;

%% Initialise path and control

X = [];%zeros(10,100);
U = [];%zeros(2,100);
count = 1;

%% Start moving and build up speed

rampTraj(3, 2); % (target_speed, ramp_time);

%% Get obstacle location and generate trajectory

disp('Waiting for obstacle data...')

p.Obs = getObstacleLocation;

% TODO: Incorporate live states

while pdist([x0(1:2)';p.xDes(1:2)]) > 0.15
    tic
    [success, x, u, cost]= iLQGDriftCar(x0, u0, p, Op);
    
    X = [X x(:,:)];
    U = [U u(:,:)];

    p.cf = 5*p.cf;
    p.cdu = 2*p.cdu;
    x0 = x(:,T/2);
    u0 = [u(:,T/2+1:end),.1*randn(2,T/2)];
    toc
end

% Execute & plot trajectory

publishTrajectory;

disp('Plotting executed trajectory')

figure(1)
init_plot(x0,p.xDes,p.Obs);

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

%% Close MATLAB-ROS global node connection

rosshutdown

%% Nested Helper functions

function obs = getObstacleLocation
    while(1)
        obs_msg = receive(obs_sub,100);
        data = obs_msg.Data(1:2);
        obs = double(data(1:2)');
        if ((obs(1)>0 || obs(2)>0) && obs(1)<2 && abs(obs(2))<1)
            obs
            break;
        end
    end
    twist_msg_sub = twist_sub.LatestMessage;
    odom_msg = odom_sub.LatestMessage;
%     amcl_msg = amcl_sub.LatestMessage;
    
    x0(4) = odom_msg.Twist.Twist.Linear.X;
    x0(5) = odom_msg.Twist.Twist.Linear.Y;
    x0(6) = odom_msg.Twist.Twist.Angular.Z;
    x0(7) = twist_msg_sub.Linear.X;
    x0(8) = twist_msg_sub.Angular.Z;

end

function rampTraj(target_speed, ramp_time)

    disp('Building speed')
    rate = rosrate(T);
    reset(rate);
    
    init_odom_msg = odom_sub.LatestMessage;
    init_x = init_odom_msg.Pose.Pose.Position.X;
    init_y = init_odom_msg.Pose.Pose.Position.Y;
    init_ori = init_odom_msg.Pose.Pose.Orientation;
    
    for i=1:T*ramp_time    
        twist_msg.Linear.X = target_speed;
        twist_msg.Angular.Z = 0;

        traj_msg.Header.FrameId = '/map';
        traj_msg.Twist.Twist = twist_msg;
        traj_msg.Pose.Pose.Position.X = target_speed/T + init_x;
        traj_msg.Pose.Pose.Position.Y = init_y;
        traj_msg.Pose.Pose.Orientation = init_ori;
       
        send(twist_chatpub,twist_msg);
        send(traj_chatpub,traj_msg);
        waitfor(rate);
    end
end

function publishTrajectory
    disp('Publishing trajectory')
    rate = rosrate(T);
    reset(rate);
    
    init_odom_msg = odom_sub.LatestMessage;
    init_x = init_odom_msg.Pose.Pose.Position.X;
    init_y = init_odom_msg.Pose.Pose.Position.Y;
    init_ori = init_odom_msg.Pose.Pose.Orientation;
    
    for i=1:size(U,2)    
        twist_msg.Linear.X = U(1,i);
        twist_msg.Angular.Z = U(2,i);

        traj_msg.Header.FrameId = '/map';
        traj_msg.Twist.Twist = twist_msg;
        traj_msg.Pose.Pose.Position.X = X(1,i) + init_x;
        traj_msg.Pose.Pose.Position.Y = X(2,i) + init_y;
        traj_msg.Pose.Pose.Orientation = init_ori;

        send(twist_chatpub,twist_msg);
        send(traj_chatpub,traj_msg);
        waitfor(rate);
    end

    twist_msg.Linear.X = 0;
    twist_msg.Angular.Z = 0;
    send(twist_chatpub,twist_msg);
end

end



