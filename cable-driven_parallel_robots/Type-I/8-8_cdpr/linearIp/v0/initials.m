% RSSLM-CDPR-Type-I initials module. The initial conditions are defined here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

% System: 8-8 CDPR with each cable modelled by multiple modified rigid finite elements and cable feed

function[] = initials()

% Importing the initial equilibrium configuration of the CDPR (from kinetostatic analyses)
load initConfig.dat initConfig;

% Global variables -- required
global n;

% Global variables -- defined
global y0 ti tf incr rtol atol int_type q dq tp;
global dmdt ddxdt ddydt ddzdt dIcxxdt dIcyydt dIczzdt dIcxydt dIcyzdt dIczxdt; 

%% Initial conditions
q = initConfig(1:n)';

% Joint velocities
dq= initConfig(n+1:2*n)';

% Initial joints positions and velocities
y0=[q; dq];

% Time span
ti=0;
tf=30;
incr=.01; %Sampling time for adaptive solver and step size for fixed step solver

%Relative and absolute tolerance of the ODE solver
rtol=1e-4;
atol=1e-6;
int_type=1; %0 for ode45, 1 for ode15s

%
tp=ti;

% Initialising (Modified in inputs.m)
dmdt=zeros(n,1);
ddxdt=zeros(n,1);
ddydt=zeros(n,1);
ddzdt=zeros(n,1);
dIcxxdt=zeros(n,1);
dIcyydt=zeros(n,1);
dIczzdt=zeros(n,1);
dIcxydt=zeros(n,1);
dIcyzdt=zeros(n,1);
dIczxdt=zeros(n,1);

end

% Since, the values of the joint angles and cable lengths are dependent, change in one of these values demands an appropriate change in the values of the other variables which satisfies the loop-closure equations. 
