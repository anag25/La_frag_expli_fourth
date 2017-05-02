function [c0,c1,c2,c3]= twod_coeffecients(mu,q,va,fa,R,eta1,eta2)
%#codegen
% function to calculate 1d equation parameters
% eta1 is the fluid1's friction
% eta2 is fluid2'2 friction

delta_eta= eta2-eta1;
eta_bar= 0.5*(eta2+eta1);

beta= (eta1*q)/(eta_bar*2*pi*R^2);
alpha=(delta_eta*mu/(4*eta_bar))+(q/(2*pi*R^2))+(eta2)*(va)/(2*eta_bar);


mm2= (0.75-0.3*R^2)+ sqrt((0.75-0.3*R^2)^2 + 0.4*(alpha-beta)*R^5 +0.2*R^2-0.3); % (maimum disturbace growth causing mode)^2

if mm2 < 0
    
    mm=1;
    
    
else
    mm =sqrt(mm2);
     
end

omega_m= -alpha +(alpha-beta+ (0.5/(R^3))-(0.75/(R^5)))*mm-((0.5/(R^3))-(1.25/(R^5)))*(mm^3)-(0.5/(R^5))*(mm^5);

 if mm == 1
     A=0; % this has to be changed to RP's version based on the notes.
 else
     A= (omega_m+beta)*R^4/((1-mm2)^2);
 end

c0= (eta1*q)/(eta_bar*2*pi*R) + A*fa;

C= -A*(1+((2*(mm2)-2.5)/R^2));

c1=A+C;
c2=-A;
c3=-0.5*A;

