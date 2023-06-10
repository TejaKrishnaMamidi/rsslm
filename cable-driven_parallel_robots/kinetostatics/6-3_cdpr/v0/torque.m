% RSSLM-CDPR-Kinetostatics torque module. The generalised forces due to spring and damper elements of the cables and damping forces of the moving platform are given here.

% Contibutors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

% System: 6-3 CDPR

function [tau_d] = torque(q, dq)

% Global variables -- required
global n ali ka ktr ca ctr fltr cdmp;
global nls1 nls2 nls3 nls4 nls5 nls6 nls7 nus1 nus2 nus3 nus4 nus5 nus6 nus7;

%% Free Simulation
tau_d=zeros(n,1);

%% Spring-damper elements

% External torques due to spring-damper elements

% sub-system 1
for ii = nls1:3:nus1
	tau_d(ii)=ktr*(fltr(ii) - q(ii)) - ctr*dq(ii);
	tau_d(ii+1)=ktr*(fltr(ii+1) - q(ii+1)) - ctr*dq(ii+1);
	%tau_d(ii)=-ctr*dq(ii); 
	%tau_d(ii+1)=-ctr*dq(ii+1); 
	tau_d(ii+2)=ka*(2*ali(1) - q(ii+2)) - ca*dq(ii+2);
end

% sub-system 2
for jj = nls2:3:nus2
	tau_d(jj)=ktr*(fltr(jj) - q(jj)) - ctr*dq(jj);
	tau_d(jj+1)=ktr*(fltr(jj+1) - q(jj+1)) - ctr*dq(jj+1);
	%tau_d(jj)=-ctr*dq(jj); 
	%tau_d(jj+1)=-ctr*dq(jj+1);
	tau_d(jj+2)=ka*(2*ali(2) - q(jj+2)) - ca*dq(jj+2);
end

% sub-system 3
for kk = nls3:3:nus3
	tau_d(kk)=ktr*(fltr(kk) - q(kk)) - ctr*dq(kk);
	tau_d(kk+1)=ktr*(fltr(kk+1) - q(kk+1)) - ctr*dq(kk+1);
	%tau_d(kk)=-ctr*dq(kk); 
	%tau_d(kk+1)=-ctr*dq(kk+1);
	tau_d(kk+2)=ka*(2*ali(3) - q(kk+2)) - ca*dq(kk+2);
end

% sub-system 4
for ll = nls4:3:nus4
	tau_d(ll)=ktr*(fltr(ll) - q(ll)) - ctr*dq(ll);
	tau_d(ll+1)=ktr*(fltr(ll+1) - q(ll+1)) - ctr*dq(ll+1);
	%tau_d(ll)=-ctr*dq(ll); 
	%tau_d(ll+1)=-ctr*dq(ll+1);
	tau_d(ll+2)=ka*(2*ali(4) - q(ll+2)) - ca*dq(ll+2);
end

% sub-system 5
for mm = nls5:3:nus5
	tau_d(mm)=ktr*(fltr(mm) - q(mm)) - ctr*dq(mm);
	tau_d(mm+1)=ktr*(fltr(mm+1) - q(mm+1)) - ctr*dq(mm+1);
	%tau_d(mm)=-ctr*dq(mm); 
	%tau_d(mm+1)=-ctr*dq(mm+1);
	tau_d(mm+2)=ka*(2*ali(5) - q(mm+2)) - ca*dq(mm+2);
end

% sub-system 6
for nn = nls6:3:nus6
	tau_d(nn)=ktr*(fltr(nn) - q(nn)) - ctr*dq(nn);
	tau_d(nn+1)=ktr*(fltr(nn+1) - q(nn+1)) - ctr*dq(nn+1);
	%tau_d(nn)=-ctr*dq(nn); 
	%tau_d(nn+1)=-ctr*dq(nn+1);
	tau_d(nn+2)=ka*(2*ali(6) - q(nn+2)) - ca*dq(nn+2);
end

% sub-system 7 (end-effector)
tau_d(nls7:nus7)=-cdmp*(dq(nls7:nus7).*dq(nls7:nus7).*sign(dq(nls7:nus7)));

end

