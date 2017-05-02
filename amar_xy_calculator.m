function [X,Y] = amar_xy_calculator(x,y,theta,theta_old,v_n,dt,DS)
%#codegen
N=length(theta);
n=1:1:N;

x0=x(1)+dt*v_n(1)*sin(theta_old(1));
y0=y(1)-dt*v_n(1)*cos(theta_old(1));

% xend=x(end)+dt*v_n(end)*sin(theta(end));
% yend=y(end)-dt*v_n(end)*cos(theta(end));
xend=x0;
yend=y0;

c= cos(theta);
s=sin(theta);


Ex = xend -x0 -(DS)*trapz(c);
Ey = yend -y0 -(DS)*trapz(s);



X=x0+DS*cumtrapz(c)+(n-1).*Ex./(N-1);
Y=y0+DS*cumtrapz(s)+(n-1).*Ey./(N-1);