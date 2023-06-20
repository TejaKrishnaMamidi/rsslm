% RSSLM-CDPR-Type-I data module. This module contains the information of geometric, elastic, and inertia parameters of the model of the mechanical system under study.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

function []=data() 

% Global variables -- defined
global nen nei ni n nc ali dof type alp a b th bt r;
global dx dy dz m g Icxx Icyy Iczz Icxy Icyz Iczx aj al ka ktr ca ctr fltr;
global crmp invmp cdmp rl0;
global nls1 nus1 nls2 nus2 nls3 nus3 nls4 nus4 nls5 nus5 nls6 nus6 nls7 nus7 nls7p2 nls7p3;

% System: 6-3 CDPR with each cable modelled by multiple modified rigid finite elements (MRFEs) with cable feed

% Number of MRFEs associated with each sub-system (cable) of the 6-3 CDPR. 
nen = 20;
nei = nen*[1;1;1;1;1;1];

% Number of links for each of the sub-systems
ni=[3*nei(1:6);6];

% Total number of links
n=sum(ni);

% Cable lengths (Half of the actual lengths)
li=[165.548038578429090884062913499; 163.725847309351392376191313313; 166.768033917156933117799048395; 164.065734366576274820655872905; 168.552581236984964950104483129; 166.526851825357787870363936235];

% Link lengths
ali=li./nei;

% Number of constraints
nc = 18;

% Degrees of freedom of the system
dof = n-nc;

% System type
type = 1; %1 for Closed-loop and 0 for open-loop

%% Initialisation

% Link lengths 
%al=[0; ali; ali; 0; ali; ali; 0; ali; ali; 0; ali; ali; 0; ali; ali; 0; ali; ali; ...; 0; 0; 0; 0; 0; 0];
al = [ali(1)*ones(ni(1),1); ali(2)*ones(ni(2),1); ali(3)*ones(ni(3),1); ali(4)*ones(ni(4),1); ali(5)*ones(ni(5),1); ali(6)*ones(ni(6),1); zeros(ni(7),1)];

% Joint values

sq3 = sqrt(3);

% sub-system 1
th11 = -0.758737815022850355972808533620;
th12 = -0.695701499275881413378761914672;
  
% sub-system 2
th21 = -2.35619449019234492884698253746;
th22 = -0.682263657597239855290615606900;

% sub-system 3
th31 = -2.67990867739865984836507770828; 
th32 = -1.54762612477102763222019671008;

% sub-system 4
th41 = -2.35619449019234492884698253746; 
th42 = -2.44686693428249145517174575741; 

% sub-system 5
th51 = -0.821809407876077521971026089376; 
th52 = -2.46043480451066546983746655070;

% sub-system 6
th61 = -0.485660861797668725372977998267; 
th62 = -1.57701323919068193344766582656;

% sub-system 7
th71 = pi/4;
th72 = pi/4;
th73 = 0;
th7x = 150;
th7y = -150;
th7z = 150*sq3;

% Locations of the vertices of the base platform

b2x = 300;
b2z = 0;
b3x = 450;
b3z = 150*sq3; 
b4x = 300;
b4z = 300*sq3;
b5x = 0;
b5z = 300*sq3;
b6x = -150;
b6z = 150*sq3;

% Dimensions (side length and circumradius) of the moving platform
trb = 8*sq3;
crmp = 8;

% Radius of the cylindrical links.
%rl=[0; rl0; rl0; 0; rl0; rl0; 0; rl0; rl0; 0; rl0; rl0; 0; rl0; rl0; 0; rl0; rl0; 0; rl0; rl0; 0; rl0; rl0; 0; rl0; rl0; 0; rl0; rl0; 0; rl0; rl0; 0; rl0; rl0; 0; 0; 0; 0; 0; 0];
rl0 = 2/100;

% DH Parameters for two MRFEs per cable
% sub-system 1 (Cable#1)
%alp= [0;   pi/2; -pi/2; -pi/2;     -pi/2;      -pi/2...];
%a  = [0;   0;     0;     0;         0;          0   ...];
%b  = [0;   0;     th13;  0;         0;          th16...];
%th = [th11; th12; 0;     th14-pi/2; th15-pi/2; -pi/2...];

% sub-system 2 (Cable#2)
%alp = [0;    pi/2; -pi/2; -pi/2;      -pi/2;      -pi/2...]; 
%a   = [b2x;  0;     0;     0;          0;          0   ...];     
%b   = [b2z;  0;     th23;  0;          0;          th26...];      
%th  = [th21; th22;  0;     th24-pi/2;  th25-pi/2; -pi/2...];       

% sub-system 3 (Cable#3)
%alp = [0;    pi/2; -pi/2; -pi/2;      -pi/2;      -pi/2...]; 
%a   = [b3x;  0;     0;     0;          0;          0   ...]; 
%b   = [b3z;  0;     th33;  0;          0;          th36...]; 
%th  = [th31; th32;  0;     th34-pi/2;  th35-pi/2; -pi/2...];

% sub-system 4 (Cable#4)
%alp = [0;    pi/2; -pi/2; -pi/2;      -pi/2;      -pi/2...]; 
%a   = [b4x;  0;     0;     0;          0;          0   ...];     
%b   = [b4z;  0;     th43;  0;          0;          th46...];      
%th  = [th41; th42;  0;     th44-pi/2;  th45-pi/2; -pi/2...];  

% sub-system 5 (Cable#5)
%alp = [0;    pi/2; -pi/2; -pi/2;      -pi/2;      -pi/2...]; 
%a   = [b5x;  0;     0;     0;          0;          0   ...];     
%b   = [b5z;  0;     th53;  0;          0;          th56...];      
%th  = [th51; th52;  0;     th54-pi/2;  th55-pi/2; -pi/2...];  

% sub-system 6 (Cable#6)
%alp = [0;    pi/2; -pi/2; -pi/2;      -pi/2;      -pi/2...]; 
%a   = [b6x;  0;     0;     0;          0;          0   ...];     
%b   = [b6z;  0;     th63;  0;          0;          th66...];      
%th  = [th61; th62;  0;     th64-pi/2;  th65-pi/2; -pi/2...];  

% sub-system 7 (Moving platform)
%alp = [0;    pi/2; pi/2; pi/2; -pi/2;  pi/2]; 
%a   = [0;    0;    0;    0;     0;     0;  ];     
%b   = [th7z; th7x; th7y; 0;     0;     0;  ];      
%th  = [pi/2; pi/2; pi/2; th71;  th72;  th73]; 

piby2=pi/2;

alp = -piby2*ones(n,1);
a   =  zeros(n,1);
b   =  zeros(n,1);
th  = -piby2*ones(n,1);

% Joint type , r=1 if revolute and r=0 if prismatic
%r=[1; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1; 0; 1; 1; 0; ...; 0; 0; 0; 1; 1; 1];
r = ones(n,1);

% Parent array
%bt=[0; 1; 2; 3; 4; 5; 0; 7; 8; 9; 10; 11; 0; 13; 14; 15; 16; 17; 
%0; 19; 20; 21; 22; 23; 0; 25; 26; 27; 28; 29; 0; 31; 32; 33; 34; 35; 0; 37; 38; 39; 40; 41];
bt = cumsum([0;ones(n-1,1)]);

% Actuated joints of open tree
%aj=[1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0 0 0]; %enter 1 for actuated joints and 0 otherwise
aj = ones(1,n);

% d - vector from origin to center-of-gravity of every link  

% all X coordinates
%dx=[0; 0;     0;     0; 0;     0;     0; 0;     0;     0; 0;     0;     0; 0;     0;     0; 0;     0;     ... 0];
% all Y coordinates
%dy=[0; ali/2; 0;     0; ali/2; 0;     0; ali/2; 0;     0; ali/2; 0;     0; ali/2; 0;     0; ali/2; 0;     ... 0];
% all Z coordinates   
%dz=[0; 0;    -ali/2; 0; 0;    -ali/2; 0; 0;    -ali/2; 0; 0;    -ali/2; 0; 0;    -ali/2; 0; 0;    -ali/2; ... 0];

dx = zeros(n,1);
dy = zeros(n,1);
dz = zeros(n,1);

% mass, moment of inertia, and acceleration due to gravity
%m=[0; m0; m0; 0; m0; m0; 0; m0; m0; 0; m0; m0; 0; m0; m0; 0; m0; m0; ... 0; 0; 0; 0; 0; mp];
g = [0; -981/100; 0];
% g=[0; 0; 0];

m = ones(n,1);

%Inertia tensor of the kth link about the center-of-mass in ith frame
%which is rigidly attached to the link
Icxx = zeros(n,1);Icyy = zeros(n,1);Iczz = zeros(n,1); % Initialization 
Icxy = zeros(n,1);Icyz = zeros(n,1);Iczx = zeros(n,1); % Initialization 

%% Assignments

% sub-system 1

nls1= 1;
nus1= ni(1);

ms10 = 935.346417968124363494955461269/nei(1);
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

alp(nls1) = 0;
alp(nls1+1) = piby2;

bt(nus1+1) = 0;

th(nls1) = th11;
th(nls1+1) = th12;
th(nls1+2) = 0;

aj(nus1-2) = 0;
aj(nus1-1) = 0;

% sub-system 2

nls2= nus1+1;
nus2= nus1+ni(2);

ms20 = 925.051037297835366925480920220/nei(2);
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
alp(nls2)  = 0;
alp(nls2+1)  = piby2;

a(nls2)    = b2x;

b(nls2)    = b2z;

bt(nus2+1) = 0;

th(nls2)   = th21;
th(nls2+1) = th22;
th(nls2+2) = 0;

aj(nus2-1) = 0;
aj(nus2-2) = 0;

% sub-system 3

nls3= nus2+1;
nus3= nus2+ni(3);

ms30 = 942.239391631936672115564623433/nei(3);
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

alp(nls3) = 0;
alp(nls3+1) = piby2;

a(nls3)=b3x;

b(nls3)=b3z;

bt(nus3+1)=0;

th(nls3)=th31;
th(nls3+1)=th32;
th(nls3+2)=0;

aj(nus3-1) =0;
aj(nus3-2) =0;

% sub-system 4
nls4= nus3+1;
nus4= nus3+ni(4);

ms40 =926.971399171155952736705681912/nei(4);
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

alp(nls4) = 0;
alp(nls4+1) = piby2;

a(nls4)=b4x;

b(nls4)=b4z;

bt(nus4+1)=0;

th(nls4)=th41;
th(nls4+1)=th42;
th(nls4+2)=0;

aj(nus4-1) =0;
aj(nus4-2) =0;

% sub-system 5
nls5= nus4+1;
nus5= nus4+ni(5);

ms50 = 952.322083988965051968090329682/nei(5);
m(nls5:nus5,1)  = ms50*ones(ni(5),1);

Ics5 = (1/12)*[ms50*(3*rl0*rl0 + ali(5)*ali(5)); 6*ms50*rl0*rl0];

for mm=nls5:3:nus5
	mmp1=mm+1;
	mmp2=mmp1+1;
	al(mm)=0;
	b(mmp2)=ali(5)*2;
	r(mmp2)=0;
	dy(mmp1)=ali(5)/2;
	dz(mmp2)=-ali(5)/2;
	m(mm)=0;
	Icxx(mmp1)=Ics5(1);
	Icyy(mmp1)=Ics5(2);
	Iczz(mmp1)=Ics5(1);
	Icxx(mmp2)=Ics5(1);
	Icyy(mmp2)=Ics5(1);
	Iczz(mmp2)=Ics5(2);
end

alp(nls5) = 0;
alp(nls5+1) = piby2;

a(nls5)=b5x;

b(nls5)=b5z;

bt(nus5+1)=0;

th(nls5)=th51;
th(nls5+1)=th52;
th(nls5+2)=0;

aj(nus5-1) =0;
aj(nus5-2) =0;

% sub-system 6
nls6= nus5+1;
nus6= nus5+ni(6);

ms60 = 940.876712813271501467556239730/nei(6);
m(nls6:nus6,1)  = ms60*ones(ni(6),1);

Ics6 = (1/12)*[ms60*(3*rl0*rl0 + ali(6)*ali(6)); 6*ms60*rl0*rl0];

for nn=nls6:3:nus6
	nnp1=nn+1;
	nnp2=nnp1+1;
	al(nn)=0;
	b(nnp2)=ali(6)*2;
	r(nnp2)=0;
	dy(nnp1)=ali(6)/2;
	dz(nnp2)=-ali(6)/2;
	m(nn)=0;
	Icxx(nnp1)=Ics6(1);
	Icyy(nnp1)=Ics6(2);
	Iczz(nnp1)=Ics6(1);
	Icxx(nnp2)=Ics6(1);
	Icyy(nnp2)=Ics6(1);
	Iczz(nnp2)=Ics6(2);
end

alp(nls6) = 0;
alp(nls6+1) = piby2;

a(nls6)=b6x;

b(nls6)=b6z;

th(nls6)=th61;
th(nls6+1)=th62;
th(nls6+2)=0;

aj(nus6-1) =0;
aj(nus6-2) =0;

% sub-system 7
nls7= nus6+1;
nus7= nus6+ni(7);
nls7p2 = nls7+2;
nls7p3 = nls7p2+1;

alp(nls7:nus7)=-alp(nls7:nus7);
alp(nls7)     =0;
alp(nus7-1)   =-piby2;

b(nls7)   = th7z;
b(nls7+1) = th7x;
b(nls7+2) = th7y;

th(nls7:nus7)=-th(nls7:nus7);
th(nus7)     =th73;
th(nus7-1)   =th72;
th(nus7-2)   =th71;

r(nls7)  =0;
r(nls7+1)=0;
r(nls7+2)=0;

bt(nls7) =0;

aj(nls7:nus7) = zeros(1,ni(7));

m(nls7:nus7) = zeros(ni(7),1); 
m(nus7)=25e3;
invmp=4e-5;

Icxx(nus7)= m(nus7)*trb^2/24;  
Iczz(nus7)= m(nus7)*trb^2/24;
Icyy(nus7)= Icxx(nus7)+Iczz(nus7);

% Longitudinal, transverse and torsional spring constants 
%ka  =nen*6.06e05;
%ktr =nen*6.06e01;
ka  =(pi*64e6)./(2*ali);
ktr =(pi*256e2)./(8*ali);

% Longitudinal and transverse damping constants
%ca=nen*1e09;
%ctr=nen*1e05;
ca=(pi*1148e8)./(2*ali);
ctr=(pi*4592e4)./(8*ali);

% Free lengths of the transverse and lateral springs
fltr = (-pi/2)*ones(n,1);

% sub-system 1
fltr(nls1) = th11;
fltr(nls1+1) = th12;

% sub-system 2
fltr(nus1+1) = th21;
fltr(nus1+2) = th22;

% sub-system 3
fltr(nus2+1) = th31; 
fltr(nus2+2) = th32;

% sub-system 4
fltr(nus3+1) = th41; 
fltr(nus3+2) = th42; 

% sub-system 5
fltr(nus4+1) = th51; 
fltr(nus4+2) = th52;

% sub-system 6
fltr(nus5+1) = th61; 
fltr(nus5+2) = th62;

% sub-system 7
% Natural: 1/2*rho*Cd*A
% cdmp = (1/2)*1.225e0*2.1e0*16*8;
% Artificial
cdmp = 0;

end
