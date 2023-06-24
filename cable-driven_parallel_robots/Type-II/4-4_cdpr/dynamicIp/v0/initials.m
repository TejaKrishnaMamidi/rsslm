% RSSLM-CDPR-Type-II-DynamicIp initials module. The initial conditions are defined here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

function[] = initials()

%System: 4-4 CDPR with cables attached to quadcopters

% Global variables -- required
global n ali nls1 nls2 nls3 nls4 nls5 nls6 nls7 nls8 nls9 nus1 nus2 nus3 nus4 nus5 nus6 nus7 nus8 nus9;

% Global variables -- defined
global y0 ti tf incr rtol atol int_type q dq;

% Initial conditions
q = -(pi/2)*ones(n,1);

% Joint values

% sub-system 1
q(nls1:nls1+2) = [19/20; 2/5; sqrt(29/5)/4];
q(nls1+3) = 0;
q(nls1+4) = pi;    

% sub-system 2
q(nls2:nls2+2) = [19/20; 1; sqrt(29/5)/4];
q(nls2+3) = 0;
q(nls2+4) = pi;

% sub-system 3
q(nls3:nls3+2) = [19/20; 2/5; (4/5)+(sqrt(29/5)/4)];
q(nls3+3) = 0; 
q(nls3+4) = pi;

% sub-system 4
q(nls4:nls4+2) = [19/20; 1; (4/5)+(sqrt(29/5)/4)];
q(nls4+3) = 0; 
q(nls4+4) = pi; 

% sub-system 5
q(nls5:nus5) = [1/10; 7/10; (2/5)+(sqrt(29/5)/4); 0; pi/2; 0];

% sub-system 6
q(nls6:nus6) = [19/20; 2/5; sqrt(29/5)/4; 0; pi/2; 0];

% sub-system 7
q(nls7:nus7) = [19/20; 1; sqrt(29/5)/4; 0; pi/2; 0];

% sub-system 8
q(nls8:nus8) = [19/20; 2/5; (4/5)+(sqrt(29/5)/4); 0; pi/2; 0];

% sub-system 9
q(nls9:nus9) = [19/20; 1; (4/5)+(sqrt(29/5)/4); 0; pi/2; 0];

% Prismatic joints

% sub-system 1
for ii = (nls1+3):3:nus1
	q(ii+2)=2*ali(1);
end

% sub-system 2
for jj = (nls2+3):3:nus2
	q(jj+2)=2*ali(2);
end

% sub-system 3
for kk = (nls3+3):3:nus3
	q(kk+2)=2*ali(3);
end

% sub-system 4
for ll = (nls4+3):3:nus4
	q(ll+2)=2*ali(4);
end

% Joint velocities
dq=zeros(n,1);

% Initial joints positions and velocities
y0=[q; dq];

% Time span
ti=0;
tf=6;
incr=.01; %Sampling time for adaptive solver and step size for fixed step solver

%Relative and absolute tolerance of the ODE solver
rtol=1e-4;
atol=1e-6;
int_type=1; %0 for ode45, 1 for ode15s
end

% Since, the values of the joint angles and cable lengths are dependent, change in one of these values demands an appropriate change in the values of the other variables which satisfies the loop-closure equations. 
