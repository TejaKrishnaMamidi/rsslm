% RSSLM-CDPR-Kinetostatics inputs module. The model parameters are defined here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Prof. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

function []=inputs() 

%System: 8-8 CDPR with each cable modelled by multiple modified rigid finite elements

% Global variables -- required
global n ali ni nei nen;

% Global variables -- defined
global dof type alp a b th bt r dx dy dz m g Icxx Icyy Iczz Icxy Icyz Iczx aj al ka ktr ca ctr fltr invmp cdmp;
global nls1 nus1 nls2 nus2 nls3 nus3 nls4 nus4 nls5 nus5 nls6 nus6 nls7 nus7 nls8 nus8 nls9 nus9 nls9p2 nls9p3;

% Degrees of freedom of the system
dof = n-24;

% System type
type = 1; %1 for Closed-loop and 0 for open-loop

%% Initialisation

% Link lengths 
%al=[0; ali; ali; 0; ali; ali; 0; ali; ali; 0; ali; ali; 0; ali; ali; 0; ali; ali; ...; 0; 0; 0; 0; 0; 0];
al = [ali(1)*ones(ni(1),1); ali(2)*ones(ni(2),1); ali(3)*ones(ni(3),1); ali(4)*ones(ni(4),1); ali(5)*ones(ni(5),1); ali(6)*ones(ni(6),1); ali(7)*ones(ni(7),1); ali(8)*ones(ni(8),1); zeros(ni(9),1)];

% Joint values

% sub-system 1
th11 = -2.64069;
th12 = 1.90743853139933871444404580066;
  
% sub-system 2
th21 = -2.53177;
th22 = 1.82504642259010141212440232789;

% sub-system 3
th31 = 2.52691; 
th32 = 1.91997797261789485715522531067;

% sub-system 4
th41 = 0.628087; 
th42 = 1.97029896903282138256176949092; 

% sub-system 5
th51 = 2.61635; 
th52 = 1.81430627761806817152247512394;

% sub-system 6
th61 = 0.760992; 
th62 = 1.87217485393891619715280331411;

% sub-system 7
th71 = -0.750432; 
th72 = 1.98537844937558431713538450558;

% sub-system 8
th81 = -0.636145; 
th82 = 1.86374846643852463456232958032;

% sub-system 9
th91 = 0;
th92 = 0;
th93 = 0;
th9x = 1;
th9y = 0;
th9z = 2;

% Locations of the vertices of the base platform

b1 = [-7.17512; -5.24398; 5.46246];
b2 = [-7.31591; -5.10296; 5.47222];
b3 = [-7.30285; 5.23598; 5.47615];
b4 = [7.18206; 5.3476; 5.4883];
b5 = [-7.16098; 5.37281; 5.48539];
b6 = [7.32331; 5.20584; 5.49903];
b7 = [7.30156; -5.13255; 5.489];
b8 = [7.16129; -5.26946; 5.49707];

% Radius of the cylindrical links.
%rl=[0; rl0; rl0; 0; rl0; rl0; ... 0; 0; 0; 0; 0; 0];
rl0 = 1/200;

%DH PARAMETERs for two MRFEs per cable
% sub-system 1 (Cable#1)
%alp= [0;   pi/2; -pi/2; -pi/2;     -pi/2;      -pi/2...];
%a  = [b1x;   0;     0;     0;         0;          0   ...];
%b  = [b1z;   0;     th13;  0;         0;          th16...];
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

% sub-system 7 (Cable#7)
%alp = [0;    pi/2; -pi/2; -pi/2;      -pi/2;      -pi/2...]; 
%a   = [b7x;  0;     0;     0;          0;          0   ...];     
%b   = [b7z;  0;     th73;  0;          0;          th76...];      
%th  = [th71; th72;  0;     th74-pi/2;  th75-pi/2; -pi/2...];  

% sub-system 8 (Cable#8)
%alp = [0;    pi/2; -pi/2; -pi/2;      -pi/2;      -pi/2...]; 
%a   = [b8x;  0;     0;     0;          0;          0   ...];     
%b   = [b8z;  0;     th83;  0;          0;          th86...];      
%th  = [th81; th82;  0;     th84-pi/2;  th85-pi/2; -pi/2...];  

% sub-system 9 (Moving platform)
%alp = [0;    pi/2; pi/2; -pi/2; pi/2;  pi/2]; 
%a   = [0;    0;    0;    0;     0;     0;  ];     
%b   = [th9z; th9x; th9y; 0;     0;     0;  ];      
%th  = [pi/2; pi/2; 0;    th91;  th92+pi/2;  th93]; 

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
%0; 19; 20; 21; 22; 23; 0; 25; 26; 27; 28; 29; 0; 31; 32; 33; 34; 35; 
%0; 37; 38; 39; 40; 41; 0; 43; 44; 45; 46; 47; 0; 49; 50; 51; 52; 53];
bt = cumsum([0;ones(n-1,1)]);

% Actuated joints of open tree
%aj=[1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0 0 0]; %enter 1 for actuated joints and 0 otherwise
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
g = [0; 0; -981/100];
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

ms10 = 1.81341/nei(1);
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

a(nls1)    = b1(1);

b(nls1)    = b1(3);

bt(nus1+1) = 0;

th(nls1) = th11;
th(nls1+1) = th12;
th(nls1+2) = 0;

aj(nus1-2) = 0;
aj(nus1-1) = 0;

% sub-system 2

nls2= nus1+1;
nus2= nus1+ni(2);

ms20 = 1.70214/nei(2);
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

a(nls2)    = b2(1);

b(nls2)    = b2(3);

bt(nus2+1) = 0;

th(nls2)   = th21;
th(nls2+1) = th22;
th(nls2+2) = 0;

aj(nus2-1) = 0;
aj(nus2-2) = 0;

% sub-system 3

nls3= nus2+1;
nus3= nus2+ni(3);

ms30 = 1.75774/nei(3);
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

a(nls3)=b3(1);

b(nls3)=b3(3);

bt(nus3+1)=0;

th(nls3)=th31;
th(nls3+1)=th32;
th(nls3+2)=0;

aj(nus3-1) =0;
aj(nus3-2) =0;

% sub-system 4
nls4= nus3+1;
nus4= nus3+ni(4);

ms40 =1.55151/nei(4);
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

a(nls4)=b4(1);

b(nls4)=b4(3);

bt(nus4+1)=0;

th(nls4)=th41;
th(nls4+1)=th42;
th(nls4+2)=0;

aj(nus4-1) =0;
aj(nus4-2) =0;

% sub-system 5
nls5= nus4+1;
nus5= nus4+ni(5);

ms50 = 1.78363/nei(5);
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

a(nls5)=b5(1);

b(nls5)=b5(3);

bt(nus5+1)=0;

th(nls5)=th51;
th(nls5+1)=th52;
th(nls5+2)=0;

aj(nus5-1) =0;
aj(nus5-2) =0;

% sub-system 6
nls6= nus5+1;
nus6= nus5+ni(6);

ms60 = 1.45694/nei(6);
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

a(nls6)=b6(1);

b(nls6)=b6(3);

bt(nus6+1)=0;

th(nls6)=th61;
th(nls6+1)=th62;
th(nls6+2)=0;

aj(nus6-1) =0;
aj(nus6-2) =0;

% sub-system 7
nls7= nus6+1;
nus7= nus6+ni(7);

ms70 = 1.49874/nei(7);
m(nls7:nus7,1)  = ms70*ones(ni(7),1);

Ics7 = (1/12)*[ms70*(3*rl0*rl0 + ali(7)*ali(7)); 6*ms70*rl0*rl0];

for oo=nls7:3:nus7
	oop1=oo+1;
	oop2=oop1+1;
	al(oo)=0;
	b(oop2)=ali(7)*2;
	r(oop2)=0;
	dy(oop1)=ali(7)/2;
	dz(oop2)=-ali(7)/2;
	m(oo)=0;
	Icxx(oop1)=Ics7(1);
	Icyy(oop1)=Ics7(2);
	Iczz(oop1)=Ics7(1);
	Icxx(oop2)=Ics7(1);
	Icyy(oop2)=Ics7(1);
	Iczz(oop2)=Ics7(2);
end

alp(nls7) = 0;
alp(nls7+1) = piby2;

a(nls7)=b7(1);

b(nls7)=b7(3);

bt(nus7+1)=0;

th(nls7)=th71;
th(nls7+1)=th72;
th(nls7+2)=0;

aj(nus7-1) =0;
aj(nus7-2) =0;

% sub-system 8
nls8= nus7+1;
nus8= nus7+ni(8);

ms80 = 1.49741/nei(8);
m(nls8:nus8,1)  = ms80*ones(ni(8),1);

Ics8 = (1/12)*[ms80*(3*rl0*rl0 + ali(8)*ali(8)); 6*ms80*rl0*rl0];

for pp=nls8:3:nus8
	ppp1=pp+1;
	ppp2=ppp1+1;
	al(pp)=0;
	b(ppp2)=ali(8)*2;
	r(ppp2)=0;
	dy(ppp1)=ali(8)/2;
	dz(ppp2)=-ali(8)/2;
	m(pp)=0;
	Icxx(ppp1)=Ics8(1);
	Icyy(ppp1)=Ics8(2);
	Iczz(ppp1)=Ics8(1);
	Icxx(ppp2)=Ics8(1);
	Icyy(ppp2)=Ics8(1);
	Iczz(ppp2)=Ics8(2);
end

alp(nls8) = 0;
alp(nls8+1) = piby2;

a(nls8)=b8(1);

b(nls8)=b8(3);

th(nls8)=th81;
th(nls8+1)=th82;
th(nls8+2)=0;

aj(nus8-1) =0;
aj(nus8-2) =0;

% sub-system 9
nls9= nus8+1;
nus9= nus8+ni(9);
nls9p2 = nls9+2;
nls9p3 = nls9p2+1;

alp(nls9:nus9)=-alp(nls9:nus9);
alp(nls9)     =0;
alp(nus9-2)   =-piby2;

b(nls9)   = th9z;
b(nls9+1) = th9x;
b(nls9+2) = th9y;

th(nls9:nus9)=-th(nls9:nus9);
th(nus9)     =th93;
th(nus9-1)   =th92+pi/2;
th(nus9-2)   =th91;
th(nus9-3)   =0;

r(nls9)  =0;
r(nls9+1)=0;
r(nls9+2)=0;

bt(nls9) =0;

aj(nls9:nus9) = zeros(1,ni(9));

m(nls9:nus9) = zeros(ni(9),1); 
m(nus9)=10;
invmp=1e-1;

Icxx(nus9)= m(nus9)*(1/12)+m(nus9)*(1/2)^2;  
Iczz(nus9)= m(nus9)*(1/12);
Icyy(nus9)= m(nus9)*(1/12)+m(nus9)*(1/2)^2;

% Data from Tempel et al.

% Icxx(nus9)= 6.032;  
% Iczz(nus9)= 3.089;
% Icyy(nus9)= 2.723;
% Icxy(nus9)= -0.215;  
% Icyz(nus9)= 0.503;
% Iczx(nus9)= -0.008;

% Longitudinal, transverse and torsional spring constants 
ka  =(nen*pi*1e5)./(8*ali.*nei);
ktr =(nen*pi*1e1)./(128*ali.*nei);

% Longitudinal and transverse damping constants
ca = (nen*pi*1e6)./(8*ali.*nei);
ctr= (nen*pi*1e2)./(128*ali.*nei);

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
fltr(nus6+1) = th71; 
fltr(nus6+2) = th72;

% sub-system 8
fltr(nus7+1) = th81; 
fltr(nus7+2) = th82;

% sub-system 9
% Natural: 1/2*rho*Cd*A
% cdmp = (1/2)*1.225e0*2.1e0;
% Artificial
cdmp = 0;

end


