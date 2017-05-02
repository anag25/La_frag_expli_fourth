function [x,y,theta,DS] = multimode_disturbed_circle(R0,N)

[x1,y1,~,~] = rp_input_disturbed_circle(R0,N,1,0.005);
[x2,y2,~,~] = rp_input_disturbed_circle(R0,N,5,0.005);
[x3,y3,~,~] = rp_input_disturbed_circle(R0,N,11,0.005);
[x4,y4,~,~] = rp_input_disturbed_circle(R0,N,15,0.005);
[x5,y5,~,~] = rp_input_disturbed_circle(R0,N,23,0.005);
[x6,y6,~,~] = rp_input_disturbed_circle(R0,N,31,0.005);

x_mean=mean([x1;x2;x3;x4;x5;x6]);
y_mean=mean([y1;y2;y3;y4;y5;y6]);
%theta_mean=mean([theta1;theta2;theta3;theta4;theta5;theta6]);
[x,y,DS]=amar_redistribute(x_mean,y_mean);

xf=[x(2:end),x(2)];yf=[y(2:end),y(2)];
xb=[x(end-1),x(1:end-1)];yb=[y(end-1),y(1:end-1)];

dx=0.5*(xf-xb);dy=0.5*(yf-yb);
theta=atan2(dy,dx)+2*pi;
dtheta_f = theta(2:N+1) - theta(1:N);
dtheta_f(N+1) = dtheta_f(1);
meandiff = mean(abs(dtheta_f));

temp = 0;
for i = 1:N+1
    theta(i) = theta(i)+temp;
    if dtheta_f(i) < -10*meandiff
        temp = 2*pi;
    elseif dtheta_f(i) > 10*meandiff
        temp = 0;
    end
end

theta=reshape(theta,[1,length(theta)]);