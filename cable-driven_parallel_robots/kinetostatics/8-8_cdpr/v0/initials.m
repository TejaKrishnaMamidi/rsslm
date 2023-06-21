% RSSLM-CDPR-Kinetostatics initials module. The initial conditions are defined here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

function[] = initials()

%System: 8-8 CDPR with each cable modelled by multiple modified rigid finite elements
 
% Global variables -- required
global n ali nls1 nls2 nls3 nls4 nls5 nls6 nls7 nls8 nls9 nus1 nus2 nus3 nus4 nus5 nus6 nus7 nus8 nus9;

% Global variables -- defined
global y0 ti tf incr rtol atol int_type q dq;

%% Initial conditions
q = -(pi/2)*ones(n,1);

% Joint values

% sub-system 1
q(nls1) = -2.64069;
q(nls1+1) = 1.90743853139933871444404580066;  

% sub-system 2
q(nls2) = -2.53177;
q(nls2+1) = 1.82504642259010141212440232789;

% sub-system 3
q(nls3) = 2.52691; 
q(nls3+1) = 1.91997797261789485715522531067;

% sub-system 4
q(nls4) = 0.628087; 
q(nls4+1) = 1.97029896903282138256176949092; 

% sub-system 5
q(nls5) = 2.61635; 
q(nls5+1) = 1.81430627761806817152247512394;

% sub-system 6
q(nls6) = 0.760992; 
q(nls6+1) = 1.87217485393891619715280331411;

% sub-system 7
q(nls7) = -0.750432; 
q(nls7+1) = 1.98537844937558431713538450558;

% sub-system 8
q(nls8) = -0.636145; 
q(nls8+1) = 1.86374846643852463456232958032;

% sub-system 9
q(nls9:nus9) = [2; 1; 0; 0; pi/2; 0];

% Prismatic joints

% sub-system 1
for ii = nls1:3:nus1
	q(ii+2)=2*ali(1);
end

% sub-system 2
for jj = nls2:3:nus2
	q(jj+2)=2*ali(2);
end

% sub-system 3
for kk = nls3:3:nus3
	q(kk+2)=2*ali(3);
end

% sub-system 4
for ll = nls4:3:nus4
	q(ll+2)=2*ali(4);
end

% sub-system 5
for mm = nls5:3:nus5
	q(mm+2)=2*ali(5);
end

% sub-system 6
for nn = nls6:3:nus6
	q(nn+2)=2*ali(6);
end

% sub-system 7
for oo = nls7:3:nus7
	q(oo+2)=2*ali(7);
end

% sub-system 8
for pp = nls8:3:nus8
	q(pp+2)=2*ali(8);
end

% Joint velocities
dq=zeros(n,1);

% Initial joints positions and velocities
y0=[q; dq];

% Time span
ti=0;
tf=30;
incr=.01; %Sampling time for adaptive solver and step size for fixed step solver

% Relative and absolute tolerance of the ODE solver
rtol=1e-4;
atol=1e-6;
int_type=1; %0 for ode45, 1 for ode15s

end

% Since, the values of the joint angles and cable lengths are dependent, change in one of these values demands an appropriate change in the values of the other variables which satisfies the loop-closure equations. 
