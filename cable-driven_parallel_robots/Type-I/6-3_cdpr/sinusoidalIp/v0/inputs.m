% RSSLM-CDPR-Type-I inputs module. The model parameters affected by the feeding and retreiving of cables are updated here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Prof. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

%System: 6-3 CDPR with each cable modelled by multiple modified rigid finite elements and (sinusoidal) cable feed

function []=inputs(time) 

% Global variables -- required
global nei nls1 nls2 nls3 nls4 nls5 nls6 nus1 nus2 nus3 nus4 nus5 nus6 rl0;

% Global variables -- defined/updated
global ali dy dz m Icxx Icyy Iczz tp ka ktr ca ctr dvi;
global dmdt ddydt ddzdt dIcxxdt dIcyydt dIczzdt; 

% Change in the cable lengths due to feed/retreival

%dtime = time-tp;
freq=pi/5;
phase=pi/10;
dl = 5e0*[sin(freq*time)-sin(freq*tp); sin(freq*time + phase)-sin(freq*tp + phase); sin(freq*time + 2*phase)-sin(freq*tp + 2*phase); sin(freq*time + 3*phase)-sin(freq*tp + 3*phase); sin(freq*time + 4*phase)-sin(freq*tp + 4*phase); sin(freq*time + 5*phase)-sin(freq*tp + 5*phase)];
dvi = (5e0*freq)*[cos(freq*time); cos(freq*time + phase); cos(freq*time + 2*phase); cos(freq*time + 3*phase); cos(freq*time + 4*phase); cos(freq*time + 5*phase)]./nei;

% Updating the parameters affected by changes in cables lengths

dli = dl./nei;
ali = ali+ dli;

ms0 = 5.65.*ali;
dms0dt = 5.65.*dvi;

% sub-system 1

Ics1 = (1/12)*[ms0(1)*(3*rl0*rl0 + ali(1)*ali(1)); 6*ms0(1)*rl0*rl0];
dIcs1dt = (1/12)*[dms0dt(1)*(3*rl0*rl0 + ali(1)*ali(1))+2*ms0(1)*ali(1)*dvi(1); 6*dms0dt(1)*rl0*rl0];

for ii=nls1:3:nus1
	iip1=ii+1;
	iip2=iip1+1;
	dy(iip1)=ali(1)/2;
	dz(iip2)=-ali(1)/2;
	ddydt(iip1)=dvi(1)/2;
	ddzdt(iip2)=-dvi(1)/2;
	m(iip1)=ms0(1);
	m(iip2)=ms0(1);
	dmdt(iip1)=dms0dt(1);
	dmdt(iip2)=dms0dt(1);
	Icxx(iip1)=Ics1(1);
	Icyy(iip1)=Ics1(2);
	Iczz(iip1)=Ics1(1);
	Icxx(iip2)=Ics1(1);
	Icyy(iip2)=Ics1(1);
	Iczz(iip2)=Ics1(2);
	dIcxxdt(iip1)=dIcs1dt(1);
	dIcyydt(iip1)=dIcs1dt(2);
	dIczzdt(iip1)=dIcs1dt(1);
	dIcxxdt(iip2)=dIcs1dt(1);
	dIcyydt(iip2)=dIcs1dt(1);
	dIczzdt(iip2)=dIcs1dt(2);
end

% sub-system 2

Ics2 = (1/12)*[ms0(2)*(3*rl0*rl0 + ali(2)*ali(2)); 6*ms0(2)*rl0*rl0];
dIcs2dt = (1/12)*[dms0dt(2)*(3*rl0*rl0 + ali(2)*ali(2))+ 2*ms0(2)*ali(2)*dvi(2); 6*dms0dt(2)*rl0*rl0];

for jj=nls2:3:nus2
	jjp1=jj+1;
	jjp2=jjp1+1;
	dy(jjp1)=ali(2)/2;
	dz(jjp2)=-ali(2)/2;
	ddydt(jjp1)=dvi(2)/2;
	ddzdt(jjp2)=-dvi(2)/2;
	m(jjp1)=ms0(2);
	m(jjp2)=ms0(2);
	dmdt(jjp1)=dms0dt(2);
	dmdt(jjp2)=dms0dt(2);
	Icxx(jjp1)=Ics2(1);
	Icyy(jjp1)=Ics2(2);
	Iczz(jjp1)=Ics2(1);
	Icxx(jjp2)=Ics2(1);
	Icyy(jjp2)=Ics2(1);
	Iczz(jjp2)=Ics2(2);
	dIcxxdt(jjp1)=dIcs2dt(1);
	dIcyydt(jjp1)=dIcs2dt(2);
	dIczzdt(jjp1)=dIcs2dt(1);
	dIcxxdt(jjp2)=dIcs2dt(1);
	dIcyydt(jjp2)=dIcs2dt(1);
	dIczzdt(jjp2)=dIcs2dt(2);
end

% sub-system 3

Ics3 = (1/12)*[ms0(3)*(3*rl0*rl0 + ali(3)*ali(3)); 6*ms0(3)*rl0*rl0];
dIcs3dt = (1/12)*[dms0dt(3)*(3*rl0*rl0 + ali(3)*ali(3))+ 2*ms0(3)*ali(3)*dvi(3); 6*dms0dt(3)*rl0*rl0];

for kk=nls3:3:nus3
	kkp1=kk+1;
	kkp2=kkp1+1;
	dy(kkp1)=ali(3)/2;
	dz(kkp2)=-ali(3)/2;
	ddydt(kkp1)=dvi(3)/2;
	ddzdt(kkp2)=-dvi(3)/2;
	m(kkp1)=ms0(3);
	m(kkp2)=ms0(3);
	dmdt(kkp1)=dms0dt(3);
	dmdt(kkp2)=dms0dt(3);
	Icxx(kkp1)=Ics3(1);
	Icyy(kkp1)=Ics3(2);
	Iczz(kkp1)=Ics3(1);
	Icxx(kkp2)=Ics3(1);
	Icyy(kkp2)=Ics3(1);
	Iczz(kkp2)=Ics3(2);
	dIcxxdt(kkp1)=dIcs3dt(1);
	dIcyydt(kkp1)=dIcs3dt(2);
	dIczzdt(kkp1)=dIcs3dt(1);
	dIcxxdt(kkp2)=dIcs3dt(1);
	dIcyydt(kkp2)=dIcs3dt(1);
	dIczzdt(kkp2)=dIcs3dt(2);
end

% sub-system 4

Ics4 = (1/12)*[ms0(4)*(3*rl0*rl0 + ali(4)*ali(4)); 6*ms0(4)*rl0*rl0];
dIcs4dt = (1/12)*[dms0dt(4)*(3*rl0*rl0 + ali(4)*ali(4))+ 2*ms0(4)*ali(4)*dvi(4); 6*dms0dt(4)*rl0*rl0];

for ll=nls4:3:nus4
	llp1=ll+1;
	llp2=llp1+1;
	dy(llp1)=ali(4)/2;
	dz(llp2)=-ali(4)/2;
	ddydt(llp1)=dvi(4)/2;
	ddzdt(llp2)=-dvi(4)/2;
	m(llp1)=ms0(4);
	m(llp2)=ms0(4);
	dmdt(llp1)=dms0dt(4);
	dmdt(llp2)=dms0dt(4);
	Icxx(llp1)=Ics4(1);
	Icyy(llp1)=Ics4(2);
	Iczz(llp1)=Ics4(1);
	Icxx(llp2)=Ics4(1);
	Icyy(llp2)=Ics4(1);
	Iczz(llp2)=Ics4(2);
	dIcxxdt(llp1)=dIcs4dt(1);
	dIcyydt(llp1)=dIcs4dt(2);
	dIczzdt(llp1)=dIcs4dt(1);
	dIcxxdt(llp2)=dIcs4dt(1);
	dIcyydt(llp2)=dIcs4dt(1);
	dIczzdt(llp2)=dIcs4dt(2);
end

% sub-system 5

Ics5 = (1/12)*[ms0(5)*(3*rl0*rl0 + ali(5)*ali(5)); 6*ms0(5)*rl0*rl0];
dIcs5dt = (1/12)*[dms0dt(5)*(3*rl0*rl0 + ali(5)*ali(5))+ 2*ms0(5)*ali(5)*dvi(5); 6*dms0dt(5)*rl0*rl0];

for mm=nls5:3:nus5
	mmp1=mm+1;
	mmp2=mmp1+1;
	dy(mmp1)=ali(5)/2;
	dz(mmp2)=-ali(5)/2;
	ddydt(mmp1)=dvi(5)/2;
	ddzdt(mmp2)=-dvi(5)/2;
	m(mmp1)=ms0(5);
	m(mmp2)=ms0(5);
	dmdt(mmp1)=dms0dt(5);
	dmdt(mmp2)=dms0dt(5);
	Icxx(mmp1)=Ics5(1);
	Icyy(mmp1)=Ics5(2);
	Iczz(mmp1)=Ics5(1);
	Icxx(mmp2)=Ics5(1);
	Icyy(mmp2)=Ics5(1);
	Iczz(mmp2)=Ics5(2);
	dIcxxdt(mmp1)=dIcs5dt(1);
	dIcyydt(mmp1)=dIcs5dt(2);
	dIczzdt(mmp1)=dIcs5dt(1);
	dIcxxdt(mmp2)=dIcs5dt(1);
	dIcyydt(mmp2)=dIcs5dt(1);
	dIczzdt(mmp2)=dIcs5dt(2);
end

% sub-system 6

Ics6 = (1/12)*[ms0(6)*(3*rl0*rl0 + ali(6)*ali(6)); 6*ms0(6)*rl0*rl0];
dIcs6dt = (1/12)*[dms0dt(6)*(3*rl0*rl0 + ali(6)*ali(6))+ 2*ms0(6)*ali(6)*dvi(6); 6*dms0dt(6)*rl0*rl0];

for nn=nls6:3:nus6
	nnp1=nn+1;
	nnp2=nnp1+1;
	dy(nnp1)=ali(6)/2;
	dz(nnp2)=-ali(6)/2;
	ddydt(nnp1)=dvi(6)/2;
	ddzdt(nnp2)=-dvi(6)/2;
	m(nnp1)=ms0(6);
	m(nnp2)=ms0(6);
	dmdt(nnp1)=dms0dt(6);
	dmdt(nnp2)=dms0dt(6);
	Icxx(nnp1)=Ics6(1);
	Icyy(nnp1)=Ics6(2);
	Iczz(nnp1)=Ics6(1);
	Icxx(nnp2)=Ics6(1);
	Icyy(nnp2)=Ics6(1);
	Iczz(nnp2)=Ics6(2);
	dIcxxdt(nnp1)=dIcs6dt(1);
	dIcyydt(nnp1)=dIcs6dt(2);
	dIczzdt(nnp1)=dIcs6dt(1);
	dIcxxdt(nnp2)=dIcs6dt(1);
	dIcyydt(nnp2)=dIcs6dt(1);
	dIczzdt(nnp2)=dIcs6dt(2);
end

% Re-initialising the previous instance of configuration vector

tp = time;

% Coefficients of stiffness and damping

ka  =(pi*64e6)./(2*ali);
ktr =(pi*256e2)./(8*ali);

ca=(pi*1148e8)./(2*ali);
ctr=(pi*4592e4)./(8*ali);

end


