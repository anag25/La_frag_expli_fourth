function [xaug,yaug,theta,ds] = rp_input_disturbed_circle(R0,N,m,h)


%#codegen
%% function to generate input x,y, theta and ds 
%R - radius of the base circle
%N - number of nodes on the circle
%m - number of lobes on the disturbed circle
%a - amplitude of disturbance and it is percentage of R  and should be less %than 1






quadfun1 = @(pt) R0 * sqrt(1 + h^2*(sin(m*pt)).^2 + 2*h*sin(m*pt) + h^2*m^2*(cos(m*pt)).^2);
L0 = integral(quadfun1,0,2*pi,'RelTol',1e-10);
quadfun2 = @(pt) 0.5 * R0^2 * (1 +  2*h*sin(m*pt) + h^2*(sin(m*pt)).^2);
A0 = integral(quadfun2,0,2*pi,'RelTol',1e-10);

A = A0;
L = L0;
Rcircle = (A0/pi)^0.5; % Radius of circle of initial area
Lcircle = 2*pi*Rcircle;% Perimeter of that circle 

% Redistribute by solving a differential equation to high precision
ds = L/N;
s = 0:ds:L;
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[~, ptheta] = ode45(@(ss,pt) 1./quadfun1(pt), s, 0, options);
ptheta = ptheta';
%ptheta(N+1)/(2*pi) - 1
r = R0*(1 + h*sin(m*ptheta));
%polar(ptheta1, r1)

xaug = r.*cos(ptheta); yaug = r.*sin(ptheta); %x & y coordinates; augmented vectors of size N+1 including N+1-th repeat node
%figure(1);plot(xaug,yaug,'ro')

dxi = 2*pi/N;   % Delta xi
sxi = L/(2*pi); % ds/ d xi 

% analytical result for tangent vector components
drdptheta = R0*h*m*cos(m*ptheta);
dpthetads = 1./quadfun1(ptheta);
tx = (drdptheta.*cos(ptheta) - yaug).*dpthetads; 
ty = (drdptheta.*sin(ptheta) + xaug).*dpthetads;
% define theta variable based on tangent vector components
theta = atan2(ty,tx) + 2*pi;

% making theta a continuous function 
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
% S=1:1:N; 
% theta=(S-1)*2*pi/N; % angle made by the each node
% A=a*R; % percentage of R: amplitude of the distubance on the circle
% r=R*(1+sin(theta*m)*A); % radius vector magnitude of each node.
% x=r.*cos(theta); % abscissa of nodes in cartisian coordinate system.
% y=r.*sin(theta);
% 
% x=[x,x(1)]; % x coordinates of nodes
% y=[y,y(1)];
% [x,y,DS]=amar_redistribute(x,y);
% 
% xf=[x(2:end),x(2)];
% yf=[y(2:end),y(2)];
%  
% xb=[x(end-1),x(1:end-1)];
% yb=[y(end-1),y(1:end-1)];
% 
% dx=0.5*(xf-xb);
% dy=0.5*(yf-yb);
% theta=atan2(dy,dx);%
% %theta=atan2(y,x);%  - this gives in the range of -pi to pi beacuse of the nature of inverse tangent.
% % but we want in 0 - 2pi range. We do the following adjustment
% %theta(theta < 0 )= theta(theta <0)+ 2*pi;
% % theta which we got here is normal angle i.e, angle made by normal at each node with repect to x axis. But we need tangnet angle as Hou et al algo is based on tangent angle.
% % In order to covert it to tangent angle add \pi/2 to it.
% %theta=theta+pi/2;
% %theta(end)=theta(1)+2*pi;
% thetaf=[theta(2:N+1),theta(2)];
% thetab=[theta(N),theta(1:N)];
% dt=0.5*(thetaf-thetab);
% 
% p = find(dt < -pi/2);
% theta(p(1)+1:N+1)=theta(p(1)+1:N+1)+2*pi;
% x=reshape(x,[1,length(x)]);
% y=reshape(y,[1,length(y)]);
% theta=reshape(theta,[1,length(theta)]);
