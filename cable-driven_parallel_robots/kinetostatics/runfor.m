% RSSLM runfor module. This module perform the simulation.
% Contibutors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to sys_ode, jacob_acc, event_termination

function []=runfor()
disp('----------------------------------------------------------------------------------');
disp('RSSLM forward dynamics ');
disp('Contibutors: Dr. Teja Krishna Mamidi and Prof. Sandipan Bandyopadhyay @IIT Madras ');
disp('Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi ');
disp('----------------------------------------------------------------------------------');

% Initialisation
data;
inputs;
initials;
jacobian_init;

% Global variables
global n y0 ti tf incr rtol atol int_type startTime fid1 fid2 fid3 pos posIndex jtor nse cflam;

% Auxillary output files
fid1= fopen('posvelacc.dat','a');
fid2= fopen('jtorque.dat','a');
fid3= fopen('lambda.dat','a');
pos = zeros(3*n+1,1000);
jtor= zeros(2*n,1000);
cflam = zeros(18,1000);
posIndex = 1;
nse = 5;

% Set the options of ODE solver
disp('Simulation Started');
startTime=tic;
overTime = tic;
option = odeset('RelTol',rtol,'AbsTol',atol,'stats','on','Jacobian',@jacob_acc, 'Events',@event_termination);

% Choosing an ODE solver
if int_type==0
    [T,Y]=ode45(@sys_ode,ti:incr:tf,y0,option);
elseif int_type==1
    [T,Y]=ode15s(@sys_ode,ti:incr:tf,y0,option);
else
    error('Select integrator correctly');
end
toc(overTime);

% Output files
fomode='w';
fip1=fopen('timevar.dat',fomode);%time
fip2=fopen('statevar.dat',fomode);%all state variables

% Wrting the configurations (states) at every time instance to the output files
for j=1:length(T)
    tsim=T(j);
    fprintf(fip1,'%.16e\n',tsim);
    fprintf(fip2,'%.16e ',Y(j,:));
    fprintf(fip2,'\n');
end

fclose('all');

disp('----------------------------------------------------------------------------------');
disp('RSSLM forward dynamics ');
disp('Contibutors: Dr. Teja Krishna Mamidi and Prof. Sandipan Bandyopadhyay @IIT Madras ');
disp('Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi ');
disp('----------------------------------------------------------------------------------');

end
