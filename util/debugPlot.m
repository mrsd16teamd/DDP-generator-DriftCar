function [ output_args ] = debugPlot(X,U,Y,dt)
%DEBUGPLOT Summary of this function goes here
%   Detailed explanation goes here
figure(3)
subplot(3,3,1)
hold on
plot(X(1,:))
plot(Y(1,:))

subplot(3,3,2)
hold on
plot(X(2,:))
plot(Y(2,:))

subplot(3,3,3)
hold on
plot(X(3,:))
plot(Y(3,:))

subplot(3,3,4)
hold on
plot(X(4,:))
plot(Y(4,:))

subplot(3,3,5)
hold on
plot(X(5,:))
plot(Y(5,:))

subplot(3,3,6)
hold on
plot(X(6,:))
plot(Y(6,:))

subplot(3,3,7)
hold on
plot(X(4,2:end)-X(4,1:end-1))
plot(Y(4,2:end)-Y(4,1:end-1))

subplot(3,3,8)
hold on
plot(X(5,2:end)-X(5,1:end-1))
plot(Y(5,2:end)-Y(5,1:end-1))

subplot(3,3,9)
hold on
plot(X(6,2:end)-X(6,1:end-1))
plot(Y(6,2:end)-Y(6,1:end-1))

end

