
close all

%% Euler time stepping based code
% over all algorithm is based on Hou et al algorithm for pdes of boundary
% this main.m script is for solving 1D theta evolution equation in zi space
%ds=s_zi*dzi; code uses discretization of zi in such a way that dzi=1;zi_min = 1 (1st node)...
%zi_max = N (last node- which matches the 1st node due to periodic boundary conditions)
% this way of zi usage is different from notes where zi_min=0, zi_max =2*pi
% therefore s_zi = ds; 
% This code is for cell fragments where area of the cell remains constant as no extra mass is being added into it externally. 
%% STEP 1 --input initial coordinates %%
% 
tic 
R=10; % radius of the cicle
R0=R; % storing radius in a different variable
N=2000; % number of nodes on the contour
m=0 ; % no. of sinusolidal lobes on the disturbed circle
a=0.01; % amplitude of the disturbance
% generate contour nodes' x,y, theta and DS for the disturbed contour
%[x,y,theta,DS] = rp_input_disturbed_circle(R0,N,m,a);
 [x,y,theta,DS] = multimode_disturbed_circle(R0,N);
% generate contour nodes' x,y, theta and DS for the base state contour 
[x_unp,y_unp,theta_unp,DS_unp]=rp_input_disturbed_circle(R0,N,0,0);
DS0=DS; %initial DS
L=DS*N; % total length of the contour.
A=amar_area(x,y); 
R=2*A/L;  

Area_base=pi*R0^2; % unperturbed area or base state area.
L_base=2*pi*R0;



%% Step 2--Driving parameters, constants for 2d system input %%
 
v_a=5; % wet active velocity or actin deploymerization velocity at the boundary
 mu=-L_base*v_a/Area_base;   % depolymerization rate = -\mu=-unifo  rm growth rate in our notes
 f_a=0; % dry active pressure or force density 
 q=0;   % point injection rate
 k_a=5000; % stiffness of area potential

 v_mu=mu*R0/2;
 v_q=q/L;  % veclocity due to injection
 [c0,c1,c2,c3]= twod_coeffecients(mu,q,v_a,f_a,R,2,0); % function to calculate 1 dimensional coefficients from 2D linear stability analysis.
 
 
 %% step 3--calculation Zone %%
 % this v_n calculation can be added to for loop instead of clalculating it
 % out side for loop
 %[v_n,v_t,H,w,b,dL_dt]= amar_velocity_repulsion2(theta,x,y,DS,c0,c1,c2,c3,A,Area_base,k_a); % %this function has soft area potential and repulsion as well
 [v_n,v_t,H,w_n,B,dL_dt] = amar_velocity(theta,DS,c0,c1,c2,c3); % this has hard area constraint. 
T0=0; % variable to calcluate time elapsed   
dt_approx= 0.1*min(DS/k_a,min(0.5*DS/max(abs(v_t)),0.25*DS^4/(max(abs(3*c3*H.^2*DS^2  + c1*DS^2 -2*c2)))));
 tmax = 0.15; nmax = round(tmax/dt_approx); nplt = floor((tmax/100)/dt_approx);

 Vn_store = zeros(201,N+1);
 H_store=zeros(201,N+1);
 X_store = zeros(201,N+1);% store x coordinate at each time step.
 Y_store = zeros(201,N+1);% store y coordinate at each time step.
 theta_store = zeros(201,N+1);% store theta coordinate at each time step.
 Area_store=zeros(201,1);
Time_store=zeros(201,1);
L_store=zeros(201,1);
 Vn_store(1,:) = v_n;
 H_store(1,:) = H;
 X_store(1,:) = x;
 Y_store(1,:) = y;
 theta_store(1,:) = theta;
 Area_store(1,1)=A;
 Time_store(1,1)=T0;
 L_store(1,1)=L;
%  Rvobj = VideoWriter('C:\Users\SJ\Google Drive\iitb\Phd\Coding\Hou et al\Cell fragments\explicit code fourth order\videos');
v = VideoWriter(sprintf('Fragment_fourth_order_softarea_filter_off_va_5_multi_%d.avi',N));

open(v);
 
 plot(x,y,'r','LineWidth',1);
 xlim([-35*R0/10 35*R0/10])
 ylim([-35*R0/10 35*R0/10])
 title({'Cell Fragment evolution'})
 xlabel('x \rightarrow')
 ylabel('y \rightarrow')

 Fr = getframe(gcf);
 writeVideo(v,Fr);
%sg_filt=sgolay(1,3);
trfun = exp(-10*((0:N)/N).^25);
p=1;
j=1;
tic
 for i = 1:nmax
%     dt=0.001;
     %estimate DS based on current length of the contour
     DS=L/N;
     if DS > 5* DS0
         fprintf('Current DS is too large when compared to initial DS.')
         break
     end
      %decide time step size based on stability criteria
     dt=0.1*min(DS/k_a,min(0.5*DS/max(abs(v_t)),0.25*DS^4/(max(abs(3*c3*H.^2*DS^2  + c1*DS^2 -2*c2))))); %
    % caluclate new theta
    [theta_new,t2,t4] = amar_evolution_4thorder_mex(theta,H,v_t,c1,c2,c3,DS,dt); % function to calculate theta using 4th order equation
    %theta_new = ifft(fft(theta_new).*trfun,'symmetric');
     %theta_new=sgolayfilt(theta_new,1,3);
    [x,y]=amar_xy_calculator_mex(x,y,theta_new,theta,v_n,dt,DS); % calculate coordinate values at new time step.
    
    theta=theta_new;
    Ar=amar_area_mex(x,y);
    R=2*Ar/L;
    %[c0,c1,c2,c3]= twod_coeffecients(mu,q,v_a,f_a,R,2,0);
%         if rem(i,10)==0
%           Area_base=Area_base +L*(v_mu+v_a+v_q)*dt;
%             [v_n,v_t,H,w,b,dL_dt,near_counter]= amar_velocity_repulsion2(theta,x,y,DS,c0,c1,c2,c3,Ar,Area_base,k_a); % calculate normal velocity, tangential velocity and H at new time step
% %         
%        else
                
%        end
   
    T0=T0+dt;% update time.
     L=L+dL_dt*dt;
    v_q=q/L;
     v_mu=mu*Ar/L;
       % [v_n,v_t,H,w,b,dL_dt]= amar_velocity_repulsion2(theta,x,y,DS,c0,c1,c2,c3,Ar,Area_base,k_a);
     [v_n,v_t,H,w,b,dL_dt]= amar_velocity_soft_mex(theta,DS,c1,c2,c3,Ar,Area_base,k_a);  % it has soft area potential but no  repulsion
  %  [v_n,v_t,H,w_n,B,dL_dt] = amar_velocity(theta,DS,c0,c1,c2,c3);
    if T0 > tmax
        break
    end
    
       if   mod(i,nplt)==0  
            L_store(j+1,1)=L;
            Vn_store(j+1,:)=v_n;
            H_store(j+1,:)=H;
            X_store(j+1,:) = x;
            Y_store(j+1,:) = y;
            theta_store(j+1,:) = theta;
            Area_store(j+1,1)=Ar;
            Time_store(j+1,1)=T0;
             plot(x,y,'g',x_unp,y_unp,'r','LineWidth',1);
             xlim([-35*R0/10 35*R0/10])
             ylim([-35*R0/10 35*R0/10])
             title({['CEll fragment evolution, no. of nodes = ',num2str(N), ', dt = ',num2str(dt)] ;['active force v_a = ',num2str(v_a),', lobes =multi , q=0  time = ',num2str(T0)]});
             xlabel('x \rightarrow')
             ylabel('y \rightarrow')
%              drawnow;
%              refreshdata(pp)
%             %pause(0.001)
            
            Fr = getframe(gcf);
            writeVideo(v,Fr);
            j=j+1;
            
       end
    
 end
 toc
 
 close(v);
 
 
 
%  Vn_storefinal=Vn_store(1:j-1,:);
%  Area_storefinal=Area_store(1:j-1,:);
%  X_storefinal=X_store(1:j-1,:);
%  Y_storefinal=Y_store(1:j-1,:);
%  Time_storefinal=Time_store(1:j-1,:);
 Area_error = (Area_store(1:j-1,:)-Area_base)./Area_base;
 figure
 f1=plot(Time_store(1:j-1,1),Area_error);
 filename = sprintf('Area_error_plot_%dnodesand va_5_.png', N);
 set(f1,'LineWidth',2)
 title({['CEll fragment evolution, no. of nodes = ',num2str(N), ', dt = ',num2str(dt)] ;['active force v_a = ',num2str(v_a),', lobes = multi , q=0  time = ',num2str(T0)]});
 xlabel('Time \rightarrow')
 ylabel('Area_error \rightarrow')
 saveas(f1,filename);
 
 figure
 f2=plot(Time_store(1:j-1,1),Vn_store(1:j-1,1));
 set(f2,'LineWidth',2)
 title({['CEll fragment evolution 1st node v_n plot, no. of nodes = ',num2str(N), ', dt = ',num2str(dt)] ;['active force v_a = ',num2str(v_a),', lobes = multi, q=0  time = ',num2str(T0)]});
 xlabel('Time \rightarrow')
 ylabel('V_n of 1st node \rightarrow')
 filename2 = sprintf('1st_node v_n_plot_%dnodes and va_5_.png', N);
 saveas(f2,filename2) 
 
  figure
 f3=plot(Time_store(1:j-1,1),L_store(1:j-1,1));
 set(f3,'LineWidth',2)
 title({['CEll fragment evolution L plot, no. of nodes = ',num2str(N), ', dt = ',num2str(dt)] ;['active force v_a = ',num2str(v_a),', lobes = multi, q=0  time = ',num2str(T0)]});
 xlabel('Time \rightarrow')
 ylabel('L \rightarrow')
 filename3 = sprintf('L_plot_%dnodes and va_5_.png', N);
 saveas(f3,filename3)