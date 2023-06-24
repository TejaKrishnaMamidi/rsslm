% RSSLM-CDPR-Type-II-DynamicIp torque module. The generalised forces due to spring and damper elements of the cables and damping forces of the moving platform are given here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

% System: 4-4 CDPR with cables attached to quadcopters

function [tau_d] = torque(q, dq, t)

% Global variables -- required
global n ali ka ktr kgr ca ctr cgr cdmp fltr nls1 nls2 nls3 nls4 nls5 nls6 nls7 nls8 nls9 nus1 nus2 nus3 nus4 nus5 nus6 nus7 nus8 nus9;

%% Free Simulation
tau_d=zeros(n,1);

%% Spring-Damper actuator element

% External forces/torques due to spring-damper elements

% sub-system 1
for ii = nls1+3:3:nus1
	tau_d(ii)=ktr(1)*(fltr(ii) - q(ii)) - ctr(1)*dq(ii);
	tau_d(ii+1)=ktr(1)*(fltr(ii+1) - q(ii+1)) - ctr(1)*dq(ii+1);
	tau_d(ii+2)=ka(1)*(2*ali(1) - q(ii+2)) - ca(1)*dq(ii+2);
end

% sub-system 2
for jj = nls2+3:3:nus2
	tau_d(jj)=ktr(2)*(fltr(jj) - q(jj)) - ctr(2)*dq(jj);
	tau_d(jj+1)=ktr(2)*(fltr(jj+1) - q(jj+1)) - ctr(2)*dq(jj+1);
	tau_d(jj+2)=ka(2)*(2*ali(2) - q(jj+2)) - ca(2)*dq(jj+2);
end

% sub-system 3
for kk = nls3+3:3:nus3
	tau_d(kk)=ktr(3)*(fltr(kk) - q(kk)) - ctr(3)*dq(kk);
	tau_d(kk+1)=ktr(3)*(fltr(kk+1) - q(kk+1)) - ctr(3)*dq(kk+1);
	tau_d(kk+2)=ka(3)*(2*ali(3) - q(kk+2)) - ca(3)*dq(kk+2);
end

% sub-system 4
for ll = nls4+3:3:nus4
	tau_d(ll)=ktr(4)*(fltr(ll) - q(ll)) - ctr(4)*dq(ll);
	tau_d(ll+1)=ktr(4)*(fltr(ll+1) - q(ll+1)) - ctr(4)*dq(ll+1);
	tau_d(ll+2)=ka(4)*(2*ali(4) - q(ll+2)) - ca(4)*dq(ll+2);
end

% sub-system 5 (moving platform)
tau_d(nls5:nus5)=-cdmp*(dq(nls5:nus5).*dq(nls5:nus5).*sign(dq(nls5:nus5)));
% Ground reaction forces
if q(nls5) <= fltr(nls5)
    tau_d(nls5)=kgr*(fltr(nls5)-q(nls5))-cgr*dq(nls5);
end

% External forces due to quadcopters

% sub-system 6 (quadcopter-1)
%c63 = cos(q(nus6)); s63=sin(q(nus6));
c62 = cos(q(nus6-1)-pi/2); s62=sin(q(nus6-1)-pi/2);
c61 = cos(q(nus6-2)); s61=sin(q(nus6-2)); 
%rotq1 = [c62*c63,               -c62*s63,               s62
%         c63*s61*s62 + c61*s63,  c61*c63 - s61*s62*s63, -s61*c62
%         -(c61*c63*s62)+s61*s63, c63*s61+c61*s62*s63,   c61*c62];

tau_d(nls6:nus6)=force_trj(1, [c61*c62; s62; -s61*c62], t);
%if q(nls6) <= fltr(nls6)
%    tau_d(nls6)=tau_d(nls6)+kgr*(fltr(nls6)-q(nls6))-cgr*dq(nls6);
%end

% sub-system 7 (quadcopter-2)
%c73 = cos(q(nus7)); s73=sin(q(nus7));
c72 = cos(q(nus7-1)-pi/2); s72=sin(q(nus7-1)-pi/2);
c71 = cos(q(nus7-2)); s71=sin(q(nus7-2)); 
%rotq2 = [c72*c73,               -c72*s73,               s72
%         c73*s71*s72 + c71*s73,  c71*c73 - s71*s72*s73, -s71*c72
%         -(c71*c73*s72)+s71*s73, c73*s71+c71*s72*s73,   c71*c72];
tau_d(nls7:nus7)=force_trj(2, [c71*c72; s72; -s71*c72], t);
%if q(nls7) <= fltr(nls7)
%    tau_d(nls7)=tau_d(nls7)+kgr*(fltr(nls7)-q(nls7))-cgr*dq(nls7);
%end

% sub-system 8 (quadcopter-3)
%c83 = cos(q(nus8)); s83=sin(q(nus8));
c82 = cos(q(nus8-1)-pi/2); s82=sin(q(nus8-1)-pi/2);
c81 = cos(q(nus8-2)); s81=sin(q(nus8-2)); 
%rotq3 = [c82*c83,               -c82*s83,               s82
%         c83*s81*s82 + c81*s83,  c81*c83 - s81*s82*s83, -s81*c82
%         -(c81*c83*s82)+s81*s83, c83*s81+c81*s82*s83,   c81*c82];
tau_d(nls8:nus8)=force_trj(3, [c81*c82; s82; -s81*c82], t);
%if q(nls8) <= fltr(nls8)
%    tau_d(nls8)=tau_d(nls8)+kgr*(fltr(nls8)-q(nls8))-cgr*dq(nls8);
%end

% sub-system 9 (quadcopter-4)
%c93 = cos(q(nus9)); s93=sin(q(nus9));
c92 = cos(q(nus9-1)-pi/2); s92=sin(q(nus9-1)-pi/2);
c91 = cos(q(nus9-2)); s91=sin(q(nus9-2)); 
%rotq4 = [c92*c93,               -c92*s93,               s92
%         c93*s91*s92 + c91*s93,  c91*c93 - s91*s92*s93, -s91*c92
%         -(c91*c93*s92)+s91*s93, c93*s91+c91*s92*s93,   c91*c92];
tau_d(nls9:nus9)=force_trj(4, [c91*c92; s92; -s91*c92], t);
%if q(nls9) <= fltr(nls9)
%    tau_d(nls9)=tau_d(nls9)+kgr*(fltr(nls9)-q(nls9))-cgr*dq(nls9);
%end

end

