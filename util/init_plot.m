function init_plot(x0,x_des,obs)
% set(gca,'xlim',[-4 4],'ylim',[-4 4],'DataAspectRatio',[1 1 1])
grid on
box on
hold all

plot([0 5],[-0.3,-0.3]);
plot([0 5],[0.9,0.9])

% Make boxes to represent start and end car poses
P = [-0.15  -0.15  0.15  0.15  -0.15; -0.08  0.08  0.08  -0.08  -0.08; 1 1 1 1 1];
start_x = x0(1);
start_y = x0(2);
start_phi = wrapToPi(x0(3));
tar_x = x_des(1);
tar_y = x_des(2);
tar_phi = wrapToPi(x_des(3));
start_A = [cos(start_phi) -sin(start_phi) start_x; sin(start_phi) cos(start_phi) start_y; 0 0 1];
tar_A = [cos(tar_phi) -sin(tar_phi) tar_x; sin(tar_phi) cos(tar_phi) tar_y; 0 0 1];
start = start_A*P;
tar = tar_A*P;

% Plot start and end
axis auto equal
plot(start(1,:),start(2,:),'color','b','linewidth',2);
plot(tar(1,:),tar(2,:),'color','r','linewidth',2);

% Plot obstacle and range at which obstacle cost is imposed
if ~isempty(obs)
    plot(obs(1),obs(2),'.','Color','r','MarkerSize',50)
    plot(obs(1),obs(2),'o','MarkerSize',75)
end
end