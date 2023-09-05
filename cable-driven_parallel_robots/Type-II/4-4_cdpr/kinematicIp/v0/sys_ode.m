% RSSLM-CDPR-Type-II-KinematicIP sys_ode module. This module finds the joint accelerations of the system under study from the equation of motion given their positions and velocities.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls are made to invdyn_tree_eff; dmdtqdot; torque; jacobian; GIM_tree; 
% acceleration; ddq_tree_eff 

% System: 4-4 CDPR with movements of cables' exit points 

function [ds] =sys_ode(t,y)

% Global variables -- required
global startTime fid1 fid2 fid3 pos jtor cflam posIndex n nc type alp a b th bt r dx dy dz m Icxx Icyy Iczz Icxy Icyz Iczx tf;

% Global variables -- defined/modified
global tu ddq tau_d q dq J lam bitrj;

% Displaying the time instance
q(1:n,1)=y(1:n);
dq(1:n,1)=y(n+1:2*n);
time=t;
fprintf('%.16e\n',time);

% Finding C dq+tug using inverse dynamics algorithm
[tu] = invdyn_tree_eff(q, dq, zeros(n,1), b, th);

% Finding input joint torque
[tau_d] = torque(q, dq, t);

% Caculation of the acceleration
if type == 1
    jacobian();
    % Computing the input accelerations of the cables' exit points
    bitrj = base_trj(time);
    [ddq]=acceleration(q, dq, tu, tau_d);
elseif type == 0
    [ddq]=ddq_tree_eff(q,n,alp,a,b,th,bt,r,dx,dy,dz,m,Icxx,Icyy,Iczz,Icxy,Icyz,Iczx,tau_d-tu);
else
    error('Enter correct type of system type, i.e,. 0 for oper-loop and 1 for closed-loop')
end

% Solving differential equation
ds=zeros(2*n,1);
ds(1:n)=dq;
ds(n+1:2*n)=ddq;

% Derivation of joint energy
%ds(2*n+1)=-tau_d'*ds(1:n);

tempTime = round(time,2);
prevTime = round(pos(1,posIndex),2);
finTime  = round(tf,2);

%disp([time, tempTime, prevTime]);

if prevTime < tempTime
    posIndex=posIndex+1;
end 

% Adding/Updating elements to the arrays pos and jtor.  
pos(1,posIndex) = time;
pos(2:n+1,posIndex) = q;
pos(n+2:3*n+1,posIndex) = ds(1:2*n);
jtor(1:n,posIndex) = tau_d;
jtor(n+1:2*n,posIndex) = J'*lam;
cflam(:,posIndex) = lam;

% Writing data for every 1 hour or if the joint accelerations are computed
% for more than 1000 time instances with a resolution of 0.01 second.
if (toc(startTime)>3600 || posIndex>999 ||tempTime==finTime)
    for ii=1:posIndex
        fprintf(fid1,"%.16e ",pos(:,ii));
        fprintf(fid1,"\n");
        fprintf(fid2,"%.16e ",jtor(:,ii));
        fprintf(fid2,"\n");
        fprintf(fid3,"%.16e ",cflam(:,ii));
        fprintf(fid3,"\n");
    end
    pos = zeros(3*n+1,1000);
    jtor = zeros(2*n,1000);
    cflam = zeros(nc,1000);
    posIndex = 1;
    pos(1,posIndex)=time;
    startTime = tic;
end

end
