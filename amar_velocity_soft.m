function [v_n,v_t,H,w_n,B,dl_dt] = amar_velocity_soft(theta,DS,c1,c2,c3,Acalc,A0,k_a)
%#codegen
% function to calculate normal velocity , tangential velocityand curvature

% v_t: tangentail velocity
% v_n: normal velocity
% H:   curvature
% dL_t: rate of change of contour length
% dA_t: rate of change of contour enclosed area


N=length(theta);


 %% calculation of H where H=\d{theta}{s}= \frac{delta theta}{delta s}
 
 
 theta_f=[theta(2:N),theta(2)]; % forward nodes with periodic boundary condition
 theta_b=[theta(N-1),theta(1:N-1)]; % backward nodes with periodic boundary condition
 dtheta= (theta_f-theta_b); % difference of forward nodes and backward nodes for the purpose of central difference scheme
 dtheta(1)=-3*theta(1)+4*theta(2)-theta(3);
 dtheta(N)=dtheta(1);
 H= -dtheta/(2*DS);    % central difference
 
 H(N)=H(1);
%  H(1)=H(N);
%% calculation of v_n
H_f=[H(2:N),H(2)]; % forward nodes
H_b=[H(N-1),H(1:N-1)]; % backward nodes

d2H_ds2 = (H_f-2*H + H_b)/(DS^2); % second order differential - central difference


d2H_ds2(N)=d2H_ds2(1); % periodic boundary

 w_n=c1*H+c2*d2H_ds2+c3*H.^3;
 %B=-trapz(w_n)/(N-1);
 
 
 B=- k_a*(Acalc/A0-1)*DS;

v_n = w_n + B; % normal velocity

%% calculation of v_t

 Hvn=H.*v_n;
 
 n=1:1:N;
 v_t1= cumtrapz(Hvn);
 v_t2=-(n-1).*trapz(Hvn)./(N-1);
 v_t= DS*(v_t1+v_t2);
 v_t(N)=v_t(1);
%% calculate dl_dt
dl_dt=trapz(dtheta.*v_n/2);