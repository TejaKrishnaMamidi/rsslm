% RSSLM-CDPR-Kinetostatics initials module. The initial conditions are defined here.

% Contibutors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

function[] = initials()

%System: 6-3 CDPR with each cable modelled by multiple modified rigid finite elements
 
% Global variables -- required
global n ali nls1 nls2 nls3 nls4 nls5 nls6 nus1 nus2 nus3 nus4 nus5 nus6 nls7 nus7;

% Global variables -- defined
global q dq y0 ti tf incr rtol atol int_type;

%% Initial conditions
q = -(pi/2)*ones(n,1);

% Joint values

% sub-system 1
q(nls1) = -0.758737815022850355972808533620;
q(nls1+1) = -0.695701499275881413378761914672;  

% sub-system 2
q(nls2) = -2.35619449019234492884698253746;
q(nls2+1) = -0.682263657597239855290615606900;

% sub-system 3
q(nls3) = -2.67990867739865984836507770828; 
q(nls3+1) = -1.54762612477102763222019671008;

% sub-system 4
q(nls4) = -2.35619449019234492884698253746; 
q(nls4+1) = -2.44686693428249145517174575741; 

% sub-system 5
q(nls5) = -0.821809407876077521971026089376; 
q(nls5+1) = -2.46043480451066546983746655070;

% sub-system 6
q(nls6) = -0.485660861797668725372977998267; 
q(nls6+1) = -1.57701323919068193344766582656;

% sub-system 7
q(nls7:nus7) = [150*sqrt(3); 150; -150; pi/4; pi/4; 0];

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

% Joint velocities
dq=zeros(n,1);

% Initial joints positions and velocities
y0=[q; dq];

% Time span
ti=0;
tf=200;
incr=.01; %Sampling time for adaptive solver and step size for fixed step solver

% Relative and absolute tolerance of the ODE solver
rtol=1e-4;
atol=1e-6;
int_type=1; %0 for ode45, 1 for ode15s

end

% Since, the values of the joint angles and cable lengths are dependent, change in one of these values demands an appropriate change in the values of the other variables which satisfies the loop-closure equations. 
