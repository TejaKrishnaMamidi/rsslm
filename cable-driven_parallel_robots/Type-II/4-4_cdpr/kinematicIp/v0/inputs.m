% RSSLM-CDPR-Type-II-KinematicIp inputs module. The model parameters are defined here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Prof. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

%System: 4-4 CDPR with each cable modelled by multiple modified rigid finite elements and movements in the exit points of cables.

function []=inputs() 

% Global variables -- required
global n nc ali ni nei; 

% Global variables -- defined/updated
global dof type alp a b th bt r dx dy dz m g Icxx Icyy Iczz Icxy Icyz Iczx aj al ka ktr kgr ca ctr cgr fltr invmp;
global nls1 nus1 nls2 nus2 nls3 nus3 nls4 nus4 nls5 nus5 cdmp;

% Degrees of freedom of the system
dof = n-nc;

% System type
type = 1; %1 for Closed-loop and 0 for open-loop

%% Initialisation

% Link lengths 
%al=[0; 0; 0; 0; ali; ali; 0; ali; ali; 0; 0; 0; 0; ali; ali; 0; ali; ali; 0; 0; 0; 0; ali; ali; 0; ali; ali; 0; 0; 0; 0; ali; ali; 0; ali; ali; 0; 0; 0; 0; 0; 0];
al = [ali(1)*ones(ni(1),1); ali(2)*ones(ni(2),1); ali(3)*ones(ni(3),1); ali(4)*ones(ni(4),1); zeros(ni(5),1)];


% Joint values

% sub-system 1
th11 = 0;
th12 = pi;
  
% sub-system 2
th21 = 0;
th22 = pi;

% sub-system 3
th31 = 0; 
th32 = pi;

% sub-system 4
th41 = 0; 
th42 = pi; 

% sub-system 5
th51 = 0;
th52 = pi/2;
th53 = 0;
th5x = 7/10;
th5y = (2/5)+(sqrt(29/5)/4);
th5z = 1/10;

% Length, height and width of the moving platform
lmp=6e-1;
hmp=2e-1;
wmp=8e-1;

% Radius of the cylindrical links.
%rl=[0; 0; 0; 0; rl0; rl0; 0; rl0; rl0; 0; 0; 0; 0; rl0; rl0; 0; rl0; rl0; 0; 0; 0; 0; rl0; rl0; 0; rl0; rl0; 0; 0; 0; 0; rl0; rl0; 0; rl0; rl0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
rl0 = 2/1000;

% DH parameters for two MRFEs per cable
% sub-system 1 (Cable#1)
%alp= [0;    pi/2; pi/2; pi/2; pi/2; -pi/2; -pi/2;     -pi/2;      -pi/2...];
%a  = [0;    0;    0;    0;    0;     0;     0;         0;          0   ...];
%b  = [b1z;  b1x;  b1y;  0;    0;     th13;  0;         0;          th16...];
%th = [pi/2; pi/2; pi/2; th11; th12;  0;     th14-pi/2; th15-pi/2; -pi/2...];

% sub-system 2 (Cable#2)
%alp = [0;    pi/2; pi/2; pi/2; pi/2; -pi/2; -pi/2;      -pi/2;      -pi/2...]; 
%a   = [0;    0;    0;    0;    0;     0;     0;          0;          0   ...];     
%b   = [b2z;  b2x;  b2y;  0;    0;     th23;  0;          0;          th26...];      
%th  = [pi/2; pi/2; pi/2; th21; th22;  0;     th24-pi/2;  th25-pi/2; -pi/2...];       

% sub-system 3 (Cable#3)
%alp = [0;    pi/2; pi/2; pi/2; pi/2; -pi/2; -pi/2;      -pi/2;      -pi/2...]; 
%a   = [0;    0;    0;    0;    0;     0;     0;          0;          0   ...]; 
%b   = [b3z;  b3x;  b3y;  0;    0;     th33;  0;          0;          th36...]; 
%th  = [pi/2; pi/2; pi/2; th31; th32;  0;     th34-pi/2;  th35-pi/2; -pi/2...];

% sub-system 4 (Cable#4)
%alp = [0;    pi/2; pi/2; pi/2; pi/2; -pi/2; -pi/2;      -pi/2;      -pi/2...]; 
%a   = [0;    0;    0;    0;    0;     0;     0;          0;          0   ...];     
%b   = [b4z;  b4x;  b4y;  0;    0;     th43;  0;          0;          th46...];      
%th  = [pi/2; pi/2; pi/2; th41; th42;  0;     th44-pi/2;  th45-pi/2; -pi/2...];  

% sub-system 5 (Moving platform)
%alp = [0;    pi/2; pi/2; -pi/2; pi/2;      pi/2]; 
%a   = [0;    0;    0;    0;     0;         0;  ];     
%b   = [th5z; th5x; th5y; 0;     0;         0;  ];      
%th  = [pi/2; pi/2; 0;    th51;  th52+pi/2; th53]; 

piby2=pi/2;

alp = -piby2*ones(n,1);
a   =  zeros(n,1);
b   =  zeros(n,1);
th  = -piby2*ones(n,1);

% Joint type , r=1 if revolute and r=0 if prismatic
%r=[0; 0; 0; 1; 1; 0; 1; 1; 0; ... 0; 0; 0; 1; 1; 0; 1; 1; 0; ... 0; 0; 0; 1; 1; 0; 1; 1; 0; ... 0; 0; 0; 1; 1; 0; 1; 1; 0; ... 0; 0; 0; 1; 1; 1];
r = ones(n,1);

% Parent array
%bt=[0; 1; 2; 3; 4; 5; 6; 7; 8; 0; 10; 11; 12; 13; 14; 15; 16; 17; 0; 19; 20; 21; 22; 23; 24; 25; 26; 0; 28; 29; 30; 31; 32; 33; 34; 35; 0; 37; 38; 39; 40; 41];
bt = cumsum([0;ones(n-1,1)]);

%Actuated joints of open tree
%aj=[1 1 1 0 0 ... 1 1 1 0 0 ... 1 1 1 0 0 ... 1 1 1 0 0 ... 0 0 0 0 0 0]; %enter 1 for actuated joints and 0 otherwise
aj = zeros(1,n);

% d - vector from origin to center-of-gravity of every link 

% all X coordinates
%dx=[0; 0; 0; 0; 0;     0;     0; 0;     0;     0; 0; 0; 0; 0;     0;     0; 0;     0;     0; 0; 0; 0; 0;     0;     0; 0;     0;     ... 0];
% all Y coordinates
%dy=[0; 0; 0; 0; ali/2; 0;     0; ali/2; 0;     0; 0; 0; 0; ali/2; 0;     0; ali/2; 0;     0; 0; 0; 0; ali/2; 0;     0; ali/2; 0;     ... 0];
% all Z coordinates   
%dz=[0; 0; 0; 0; 0;    -ali/2; 0; 0;    -ali/2; 0; 0; 0; 0; 0;    -ali/2; 0; 0;    -ali/2; 0; 0; 0; 0; 0;    -ali/2; 0; 0;    -ali/2; ... 0];

dx = zeros(n,1);
dy = zeros(n,1);
dz = zeros(n,1);

% mass, moment of inertia, and acceleration due to gravity
%m=[0; 0; 0; 0; m01; m01; 0; m01; m01; ... 0; 0; 0; 0; m02; m02; 0; m02; m02; ... 0; 0; 0; 0; m03; m03; 0; m03; m03; ... 0; 0; 0; 0; m03; m03; 0; m03; m03; ... 0; 0; 0; 0; 0; mp];
g = [0; 0; -981/100];
% g=[0; 0; 0];

m = zeros(n,1);

%Inertia tensor of the kth link about the center-of-mass in ith frame
%which is rigidly attached to the link
Icxx = zeros(n,1);Icyy = zeros(n,1);Iczz = zeros(n,1); % Initialization 
Icxy = zeros(n,1);Icyz = zeros(n,1);Iczx = zeros(n,1); % Initialization 

%% Assignments

% sub-system 1

nls1= 1;
nus1= ni(1);

ms10 = (25/1000)/nei(1);
m(nls1:nus1,1)  = ms10*ones(ni(1),1);

Ics1 = (1/12)*[ms10*(3*rl0*rl0 + ali(1)*ali(1)); 6*ms10*rl0*rl0];

for ii=nls1:3:nus1
	iip1=ii+1;
	iip2=iip1+1;
	al(ii)=0;
	b(iip2)=2*ali(1);
	r(iip2)=0; 
	dy(iip1)=ali(1)/2;
	dz(iip2)=-ali(1)/2;
	m(ii)=0;
	Icxx(iip1)=Ics1(1);
	Icyy(iip1)=Ics1(2);
	Iczz(iip1)=Ics1(1);
	Icxx(iip2)=Ics1(1);
	Icyy(iip2)=Ics1(1);
	Iczz(iip2)=Ics1(2);
end

al(nls1+1:nls1+2) = zeros(2,1);

alp(nls1) = 0;
alp(nls1+1:nls1+4) = piby2*ones(4,1);

b(nls1:nls1+2) = [19/20; 2/5; sqrt(29/5)/4];

th(nls1:nls1+2) = piby2*ones(3,1);
th(nls1+3:nls1+5) = [th11; th12; 0];

r(nls1:nls1+1)=zeros(2,1);

bt(nus1+1) = 0;

dy(nls1+1) = 0;
dz(nls1+2) = 0;

m(nls1+1:nls1+2) = zeros(2,1);
Icxx(nls1+1:nls1+2)=zeros(2,1);
Icyy(nls1+1:nls1+2)=zeros(2,1);
Iczz(nls1+1:nls1+2)=zeros(2,1);

% sub-system 2

nls2= nus1+1;
nus2= nus1+ni(2);

ms20 = (25/1000)/nei(2);
m(nls2:nus2,1)  = ms20*ones(ni(2),1);

Ics2 = (1/12)*[ms20*(3*rl0*rl0 + ali(2)*ali(2)); 6*ms20*rl0*rl0];

for jj=nls2:3:nus2
	jjp1=jj+1;
	jjp2=jjp1+1;
	al(jj)=0;
	b(jjp2)=2*ali(2);
	r(jjp2)=0;
	dy(jjp1)=ali(2)/2;
	dz(jjp2)=-ali(2)/2;
	m(jj)=0;
	Icxx(jjp1)=Ics2(1);
	Icyy(jjp1)=Ics2(2);
	Iczz(jjp1)=Ics2(1);
	Icxx(jjp2)=Ics2(1);
	Icyy(jjp2)=Ics2(1);
	Iczz(jjp2)=Ics2(2);
end

al(nls2+1:nls2+2) = zeros(2,1);

alp(nls2) = 0;
alp(nls2+1:nls2+4) = piby2*ones(4,1);

b(nls2:nls2+2) = [19/20; 1; sqrt(29/5)/4];

th(nls2:nls2+2) = piby2*ones(3,1);
th(nls2+3:nls2+5) = [th21; th22; 0];

r(nls2:nls2+1)=zeros(2,1);

bt(nus2+1) = 0;

dy(nls2+1) = 0;
dz(nls2+2) = 0;

m(nls2+1:nls2+2) = zeros(2,1);
Icxx(nls2+1:nls2+2)=zeros(2,1);
Icyy(nls2+1:nls2+2)=zeros(2,1);
Iczz(nls2+1:nls2+2)=zeros(2,1);

% sub-system 3

nls3= nus2+1;
nus3= nus2+ni(3);

ms30 = (25/1000)/nei(3);
m(nls3:nus3,1)  = ms30*ones(ni(3),1);

Ics3 = (1/12)*[ms30*(3*rl0*rl0 + ali(3)*ali(3)); 6*ms30*rl0*rl0];

for kk=nls3:3:nus3
	kkp1=kk+1;
	kkp2=kkp1+1;
	al(kk)=0;
	b(kkp2)=2*ali(3);
	r(kkp2)=0;
	dy(kkp1)=ali(3)/2;
	dz(kkp2)=-ali(3)/2;
	m(kk)=0;
	Icxx(kkp1)=Ics3(1);
	Icyy(kkp1)=Ics3(2);
	Iczz(kkp1)=Ics3(1);
	Icxx(kkp2)=Ics3(1);
	Icyy(kkp2)=Ics3(1);
	Iczz(kkp2)=Ics3(2);
end

al(nls3+1:nls3+2) = zeros(2,1);

alp(nls3) = 0;
alp(nls3+1:nls3+4) = piby2*ones(4,1);

b(nls3:nls3+2) = [19/20; 2/5; (4/5)+(sqrt(29/5)/4)];

th(nls3:nls3+2) = piby2*ones(3,1);
th(nls3+3:nls3+5) = [th31; th32; 0];

r(nls3:nls3+1)=zeros(2,1);

bt(nus3+1) = 0;

dy(nls3+1) = 0;
dz(nls3+2) = 0;

m(nls3+1:nls3+2) = zeros(2,1);
Icxx(nls3+1:nls3+2)=zeros(2,1);
Icyy(nls3+1:nls3+2)=zeros(2,1);
Iczz(nls3+1:nls3+2)=zeros(2,1);

% sub-system 4
nls4= nus3+1;
nus4= nus3+ni(4);

ms40 = (25/1000)/nei(4);
m(nls4:nus4,1)  = ms40*ones(ni(4),1);

Ics4 = (1/12)*[ms40*(3*rl0*rl0 + ali(4)*ali(4)); 6*ms40*rl0*rl0];

for ll=nls4:3:nus4
	llp1=ll+1;
	llp2=llp1+1;
	al(ll)=0;
	b(llp2)=ali(4)*2;
	r(llp2)=0;
	dy(llp1)=ali(4)/2;
	dz(llp2)=-ali(4)/2;
	m(ll)=0;
	Icxx(llp1)=Ics4(1);
	Icyy(llp1)=Ics4(2);
	Iczz(llp1)=Ics4(1);
	Icxx(llp2)=Ics4(1);
	Icyy(llp2)=Ics4(1);
	Iczz(llp2)=Ics4(2);
end

al(nls4+1:nls4+2) = zeros(2,1);

alp(nls4) = 0;
alp(nls4+1:nls4+4) = piby2*ones(4,1);

b(nls4:nls4+2) = [19/20; 1; (4/5)+(sqrt(29/5)/4)];

th(nls4:nls4+2) = piby2*ones(3,1);
th(nls4+3:nls4+5) = [th41; th42; 0];

r(nls4:nls4+1)=zeros(2,1);

bt(nus4+1) = 0;

dy(nls4+1) = 0;
dz(nls4+2) = 0;

m(nls4+1:nls4+2) = zeros(2,1);
Icxx(nls4+1:nls4+2)=zeros(2,1);
Icyy(nls4+1:nls4+2)=zeros(2,1);
Iczz(nls4+1:nls4+2)=zeros(2,1);

% sub-system 5
nls5= nus4+1;
nus5= nus4+ni(5);

alp(nls5+1:nus5-3)=-alp(nls5+1:nus5-3);
alp(nus5-1:nus5)=-alp(nus5-1:nus5);
alp(nls5)     =0;

b(nls5:nls5+2)   = [th5z; th5x; th5y];

th(nls5:nls5+1)=-th(nls5:nls5+1);
th(nus5-3:nus5) = [0; th51; th52+(pi/2); th53];

r(nls5:nls5+2)=zeros(3,1);

bt(nus5+1) = 0;

m(nus5)=5e-1;
invmp = 2;

Icxx(nus5)= (m(nus5)/12)*(hmp^2+wmp^2);
Icyy(nus5)= (m(nus5)/12)*(lmp^2+hmp^2);  
Iczz(nus5)= (m(nus5)/12)*(lmp^2+wmp^2);


% Longitudinal, transverse and torsional spring constants 
ka  =(pi*556e3)./(2*ali);
ktr =(pi*2224e-3)./(8*ali);

% Longitudinal and transverse damping constants
ca=(pi*1148e1)./(2*ali);
ctr=(pi*4592e-5)./(8*ali);

% sub-system 5
% Natural: 1/2*rho*Cd*A
% Artificial
cdmp = 0;

% Ground properties
kgr = 1e3;
cgr = 1e2;

% Free lengths of the transverse and lateral springs
fltr = (-pi/2)*ones(n,1);

% sub-system 1
fltr(nls1+3) = th11;
fltr(nls1+4) = th12;

% sub-system 2
fltr(nls2+3) = th21;
fltr(nls2+4) = th22;

% sub-system 3
fltr(nls3+3) = th31; 
fltr(nls3+4) = th32;

% sub-system 4
fltr(nls4+3) = th41; 
fltr(nls4+4) = th42; 

% sub-system 5
fltr(nls5) = th5z;

end


