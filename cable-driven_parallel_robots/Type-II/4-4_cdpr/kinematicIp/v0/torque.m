% RSSLM-CDPR-Type-II-KinematicIp torque module. The generalised forces due to spring and damper elements of the cables and damping forces of the moving platform are given here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

% System: 4-4 CDPR with movements of cables' exit points

function [tau_d] = torque(q, dq, t)

% Global variables -- required
global n ali ka ktr kgr ca ctr cgr cdmp fltr nls1 nls2 nls3 nls4 nls5 nus1 nus2 nus3 nus4 nus5;

%% Free Simulation
tau_d=zeros(n,1);

%% Spring-Damper elements

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

end

