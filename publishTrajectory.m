function publishTrajectory(twist_chatpub,twist_msg,traj_chatpub,traj_msg,X,U)
    disp('Publishing trajectory')
    T=50;
    rate = rosrate(T);
    reset(rate);
    for i=1:size(U,2)
        i
        twist_msg.Linear.X = U(1,i);
        twist_msg.Angular.Z = U(2,i);

        traj_msg.Header.FrameId = '/map';
        traj_msg.Twist.Twist = twist_msg;
        traj_msg.Pose.Pose.Position.X = X(1,i);
        traj_msg.Pose.Pose.Position.Y = X(2,i);

        send(twist_chatpub,twist_msg)
        send(traj_chatpub,traj_msg)
        waitfor(rate);
    end

    twist_msg.Linear.X = 0;
    twist_msg.Angular.Z = 0;
    send(twist_chatpub,twist_msg);
end
