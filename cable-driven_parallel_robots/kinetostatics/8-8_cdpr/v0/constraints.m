% RSSLM-CDPR-Kinetostatics constraints module. The values of constraint functions and their derivatives are computed.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

function [eta, deta]=constraints(il, iu, q, dq, tt, tb, so, st, vt, Qf, p, b, th)

% Depends on correctionYcoordinate

% System: 8-8 CDPR

% Global variables -- required
global alp a bt r dx dy dz;
global nus1 nus2 nus3 nus4 nus5 nus6 nus7 nus8 n;
global mpa mpv;

% Computation of the changes in the constraint functions 
for ii=il:iu
	%p(ii)=1-r(ii);
	th(ii)=th(ii)*p(ii)+q(ii)*r(ii);
	b(ii)=b(ii)*r(ii)+q(ii)*p(ii);
	cth=cos(th(ii)); calp=cos(alp(ii));
	sth=sin(th(ii)); salp=sin(alp(ii));
	if bt(ii)==0 %When parent of the link is ground link	
    	Qi=[cth,      -sth,       0
       	calp*sth,  calp*cth, -salp
       	salp*sth,  salp*cth,  calp];
    	Qf(:,:,ii)=Qi;
        
    	%Positions
    	di=[dx(ii);dy(ii);dz(ii)]+[dx(ii);dy(ii);dz(ii)];
        aim=[a(ii)
             correctionYcoordinate(ii)-b(ii)*salp
             b(ii)*calp];
    	so(:,ii)=aim;
    	st(:,ii)=so(:,ii)+Qf(:,:,ii)*di;
       
    	%w angular velocity
     	tt(:,ii)=[0;0;r(ii)*dq(ii)];
     	tti=tt(:,ii);
        
    	%v linear velocity
     	tb(:,ii)=[0;0;p(ii)*dq(ii)];
     	vt(:,ii)=[0;-salp*p(ii)*dq(ii);calp*p(ii)*dq(ii)]-Qi(:,1)*di(2)*tti(3)+Qi(:,2)*di(1)*tti(3);
        
	else %Calculation for the links other than those attached with ground
     	Qi=[cth,      -sth,       0
         calp*sth,  calp*cth, -salp
         salp*sth,  salp*cth,  calp];
     	Qf(:,:,ii)=Qf(:,:,bt(ii))*Qi;
        
     	%position vector from origin of link to origin of next link
     	aim=[a(ii)
            -b(ii)*salp
             b(ii)*calp];
     	di=[dx(ii);dy(ii);dz(ii)]+[dx(ii);dy(ii);dz(ii)];
        
     	%Positions
     	so(:,ii)=so(:,bt(ii))+Qf(:,:,bt(ii))*aim;
     	st(:,ii)=so(:,ii)+Qf(:,:,ii)*di;
        
     	%w angular velocity
     	ttbi=tt(:,bt(ii));
     	tt(:,ii)=Qi'*ttbi+[0;0;r(ii)*dq(ii)];
     	tti=tt(:,ii);
     	%v  linear velocity
     	ttbixaim=[ttbi(2)*aim(3)-aim(2)*ttbi(3);-(ttbi(1)*aim(3)-aim(1)*ttbi(3));ttbi(1)*aim(2)-aim(1)*ttbi(2)];
     	tb(:,ii)=Qi.'*(tb(:,bt(ii))+ttbixaim)+[0;0;p(ii)*dq(ii)];
     	ttixdi=[tti(2)*di(3)-di(2)*tti(3);-(tti(1)*di(3)-di(1)*tti(3));tti(1)*di(2)-di(1)*tti(2)];
     	vt(:,ii)=Qf(:,:,ii)*(tb(:,ii)+ttixdi);
	end
end

% Only if the variables q(n-6) to q(n) are varied the following computations are done.  

if il>nus8
    c73 = cos(q(n)); s73=sin(q(n));
    c72 = cos(q(n-1)-pi/2); s72=sin(q(n-1)-pi/2);
    c71 = cos(q(n-2)); s71=sin(q(n-2));
    dth73 = dq(n); dth72 = dq(n-1); dth71 = dq(n-2);
    
    rotmp = [c72*c73,                -c72*s73,                  s72
	     c73*s71*s72 + c71*s73,  c71*c73 - s71*s72*s73,    -s71*c72
	     -(c71*c73*s72)+s71*s73, c73*s71+c71*s72*s73,       c71*c72];

    rotmpDot = [-c73*s72*dth72-c72*s73*dth73, s72*s73*dth72-c72*c73*dth73, c72*dth72
	         c71*c73*s72*dth71 - s71*s73*dth71 + c72*c73*s71*dth72 + c71*c73*dth73 - s71*s72*s73*dth73, -c73*s71*dth71 -  c71*s72*s73*dth71 - c72*s71*s73*dth72 -  c73*s71*s72*dth73 - c71*s73*dth73, -c71*c72*dth71 + s71*s72*dth72
	         c73*s71*s72*dth71 + c71*s73*dth71 - c71*c72*c73*dth72 + c73*s71*dth73 + c71*s72*s73*dth73, c71*c73*dth71 -  s71*s72*s73*dth71 + c71*c72*s73*dth72 + c71*c73*s72*dth73 - s71*s73*dth73, -c72*s71*dth71 - c71*s72*dth72];

    mpcom = [q(n-4);q(n-3);q(n-5)];
    mpvc  = [dq(n-4);dq(n-3);dq(n-5)];
    
    posa1 = mpcom + rotmp*[0.503210; -0.492830; 0];
    posa2 = mpcom + rotmp*[-0.509740; 0.350900; 0.997530];
    posa3 = mpcom + rotmp*[-0.503210; -0.269900; 0];
    posa4 = mpcom + rotmp*[-0.503210; 0.492830; 0];
    posa5 = mpcom + rotmp*[0.496070; 0.355620; 0.999540];
    posa6 = mpcom + rotmp*[0.499640; -0.340280; 0.999180];
    posa7 = mpcom + rotmp*[0.502090; 0.274900; -0.000620]; 
    posa8 = mpcom + rotmp*[-0.504540; -0.346290; 0.997520]; 

    vela1 = mpvc + rotmpDot*[0.503210; -0.492830; 0];
    vela2 = mpvc + rotmpDot*[-0.509740; 0.350900; 0.997530];
    vela3 = mpvc + rotmpDot*[-0.503210; -0.269900; 0];
    vela4 = mpvc + rotmpDot*[-0.503210; 0.492830; 0];
    vela5 = mpvc + rotmpDot*[0.496070; 0.355620; 0.999540];
    vela6 = mpvc + rotmpDot*[0.499640; -0.340280; 0.999180];
    vela7 = mpvc + rotmpDot*[0.502090; 0.274900; -0.000620];
    vela8 = mpvc + rotmpDot*[-0.504540; -0.346290; 0.997520];
else
    posa1 = mpa(1:3,1);
    posa2 = mpa(1:3,2);
    posa3 = mpa(1:3,3);
    posa4 = mpa(1:3,4);
    posa5 = mpa(1:3,5);
    posa6 = mpa(1:3,6);
    posa7 = mpa(1:3,7);
    posa8 = mpa(1:3,8);

    vela1 = mpv(1:3,1);
    vela2 = mpv(1:3,2);
    vela3 = mpv(1:3,3);
    vela4 = mpv(1:3,4);
    vela5 = mpv(1:3,5);
    vela6 = mpv(1:3,6);
    vela7 = mpv(1:3,7);
    vela8 = mpv(1:3,8);
end

% Velocities of the points at which the sub-systems are joined with the moving platform

vonus1 = Qf(:,:,nus1)*tb(:,nus1);
vonus2 = Qf(:,:,nus2)*tb(:,nus2);
vonus3 = Qf(:,:,nus3)*tb(:,nus3);
vonus4 = Qf(:,:,nus4)*tb(:,nus4);
vonus5 = Qf(:,:,nus5)*tb(:,nus5);
vonus6 = Qf(:,:,nus6)*tb(:,nus6);
vonus7 = Qf(:,:,nus7)*tb(:,nus7);
vonus8 = Qf(:,:,nus8)*tb(:,nus8);

% sub-systems 1 and 9
eta1 = so(:,nus1)-posa1;
deta1 = vonus1-vela1;

% sub-systems 2 and 9
eta2 = so(:,nus2)-posa2;
deta2 = vonus2-vela2;

% sub-systems 3 and 9
eta3 = so(:,nus3)-posa3;
deta3 = vonus3-vela3;

% sub-systems 4 and 9
eta4 = so(:,nus4)-posa4;
deta4 = vonus4-vela4;

% sub-systems 5 and 9
eta5 = so(:,nus5)-posa5;
deta5 = vonus5-vela5;

% sub-systems 6 and 9
eta6 = so(:,nus6)-posa6;
deta6 = vonus6-vela6;

% sub-systems 7 and 9
eta7 = so(:,nus7)-posa7;
deta7 = vonus7-vela7;

% sub-systems 8 and 9
eta8 = so(:,nus8)-posa8;
deta8 = vonus8-vela8;

% Position constraints
eta=[eta1;eta2;eta3;eta4;eta5;eta6;eta7;eta8];

% Velocity constraints
deta=[deta1;deta2;deta3;deta4;deta5;deta6;deta7;deta8];

end

