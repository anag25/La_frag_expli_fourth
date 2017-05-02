   function [theta,d2theta,d4theta] = amar_evolution_4thorder(theta,H,v_t,c1,c2,c3,DS,dt)
   %#codegen
  % function for estimating theta at next time step based on 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %   dtheta       d{c0}        d^2{theta}         d^4{theta}            d{theta}     d^2{theta}          d{theta}   %
  %  -------- = - ------- + c1*------------- + c2*------------  + 3*c3*{---------}^2*------------- +v_t*----------   %                        
  %    dt           ds             ds^2              ds^4                  ds            ds^2               ds       %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 N=length(theta);
theta=reshape(theta,[1,N]);
H=reshape(H,[1,N]);
H2=H.*H;
%hv=heaviside(v_t);% heaveside function for upwinding purpose. Gives 0 or 1 based on sign on v_t
 hv= (1+sign(v_t))*0.5;


theta_f=[theta(2:N),theta(2)];
theta_b=[theta(N-1),theta(1:N-1)];
offset=theta(N)-theta(1);
% dtheta=zeros(1,N);
% dtheta(1)=0.5*(theta(2)-(theta(N-1)-offset));
% dtheta(2:N-1)= 0.5*((theta(3:N)-theta(1:N-2)));
% dtheta(N)=dtheta(1);


d2theta= zeros(1,N);
d2theta(1)=theta(2)-2*theta(1)+(theta(N-1)-offset);
d2theta(2:N-1)= theta(3:N)-2*theta(2:N-1)+theta(1:N-2);
d2theta(N)=d2theta(1);


d4theta= zeros(1,N);

d4theta(3:N-2)= theta(5:N)-4*theta(4:N-1)+6*theta(3:N-2)-4*theta(2:N-3)+theta(1:N-4);
d4theta(2)=theta(4)-4*theta(3)+6*theta(2)-4*theta(1)+(theta(N-1)-offset);
d4theta(N-1)=(theta(2)+offset)-4*theta(N)+6*theta(N-1)-4*theta(N-2)+theta(N-3);
d4theta(1)=theta(3)-4*theta(2)+6*theta(1)-4*(theta(N-1)-offset)+(theta(N-2)-offset);
d4theta(N)=d4theta(1);



dtheta_f=theta_f -theta;
dtheta_f(N)=dtheta_f(1); % periodic boundary condition

dtheta_b= theta-theta_b;
dtheta_b(1)=dtheta_b(N);


dtheta_upwind=(1-hv).*dtheta_b +(hv).*dtheta_f;


dtheta_dt =   c1*d2theta/(DS^2)+ 3*c3*H2.*d2theta/(DS^2) + c2*d4theta/(DS)^4   + (v_t.*dtheta_upwind)/DS; % since co is 0, it is not taken into account

theta=theta+dt*dtheta_dt;






















%  theta2f=[theta(3:N),theta(2:3)];
% 
%  
%  theta2b=[theta(end-2:end-1),theta(1:end-2)];


% v_nf=[v_n(2:end),v_n(2)];
% v_nb=[v_n(end-1),v_n(1:end-1)];
% 
% dv_n = (v_nf-v_nb)*0.5;

% 
% d1t1=(3*theta-4*theta_b+theta2b)/2; % \d{theta}{zeta} upwind scheme: this is  used when v_t < 0
% d1t1(1)=d1t1(N-1);
% d1t1(N)=d1t1(2);


% d1t2=(-theta2f+4*theta_f-3*theta)/2; % \d{theta}{zeta} upwind scheme: this is used when v_t > 0
% d1t2(1)=d1t2(N-1);
% d1t2(N)=d1t2(2);
