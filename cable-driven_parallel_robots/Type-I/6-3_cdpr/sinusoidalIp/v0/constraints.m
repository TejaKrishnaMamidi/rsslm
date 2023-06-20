% RSSLM-CDPR-Type-I constraints module. The values of constraint functions and their derivatives are computed.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

function [eta, deta]=constraints(il, iu, q, dq, tt, tb, so, Qf, p, b, th)

% No function calls

% System: 6-3 CDPR with (linear or sinusoidal) cable feed

% Global variables -- required
global alp a bt r;
global nus1 nus2 nus3 nus4 nus5 nus6 n;
global mpa1 mpa2 mpa3 crmp mpv1 mpv2 mpv3;

% Computation of the change in the constraint equations
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
    	%di=[dx(ii);dy(ii);dz(ii)]+[dx(ii);dy(ii);dz(ii)];
        aim=[a(ii)
            -b(ii)*salp
             b(ii)*calp];
    	so(:,ii)=aim;
    	%st(:,ii)=so(:,ii)+Qf(:,:,ii)*di;
       
    	%w angular velocity
     	tt(:,ii)=[0;0;r(ii)*dq(ii)];
     	%tti=tt(:,ii);
        
    	%v linear velocity
     	tb(:,ii)=[0;0;p(ii)*dq(ii)];
     	%vt(:,ii)=[0;-salp*p(ii)*dq(ii);calp*p(ii)*dq(ii)]-Qi(:,1)*di(2)*tti(3)+Qi(:,2)*di(1)*tti(3);
        
	else %Calculation for the links other than those attached with ground
     	Qi=[cth,      -sth,       0
         calp*sth,  calp*cth, -salp
         salp*sth,  salp*cth,  calp];
     	Qf(:,:,ii)=Qf(:,:,bt(ii))*Qi;
        
     	%position vector from origin of link to origin of next link
     	aim=[a(ii)
            -b(ii)*salp
             b(ii)*calp];
     	%di=[dx(ii);dy(ii);dz(ii)]+[dx(ii);dy(ii);dz(ii)];
        
     	%Positions
     	so(:,ii)=so(:,bt(ii))+Qf(:,:,bt(ii))*aim;
     	%st(:,ii)=so(:,ii)+Qf(:,:,ii)*di;
        
     	%w angular velocity
     	ttbi=tt(:,bt(ii));
     	tt(:,ii)=Qi'*ttbi+[0;0;r(ii)*dq(ii)];
     	%tti=tt(:,ii);
     	%v  linear velocity
     	ttbixaim=[ttbi(2)*aim(3)-aim(2)*ttbi(3);-(ttbi(1)*aim(3)-aim(1)*ttbi(3));ttbi(1)*aim(2)-aim(1)*ttbi(2)];
     	tb(:,ii)=Qi.'*(tb(:,bt(ii))+ttbixaim)+[0;0;p(ii)*dq(ii)];
     	%ttixdi=[tti(2)*di(3)-di(2)*tti(3);-(tti(1)*di(3)-di(1)*tti(3));tti(1)*di(2)-di(1)*tti(2)];
     	%vt(:,ii)=Qf(:,:,ii)*(tb(:,ii)+ttixdi);
	end
end

% Only if the variables q(n-6) to q(n) are varied the following computations are done.  

if il>nus6
    c73 = cos(q(n)); s73=sin(q(n));
    c72 = cos(q(n-1)); s72=sin(q(n-1));
    c71 = cos(q(n-2)); s71=sin(q(n-2));
    dth73 = dq(n); dth72 = dq(n-1); dth71 = dq(n-2);
    
    rotmpcol1 = [c71*c72*c73 - s71*s73; c72*c73*s71 + c71*s73; -(c73*s72)];
    rotmpcol3 = [c71*s72; s71*s72; c72];

    rotmpDotcol1 = [-(c73*(c72*dth71*s71 + dth73*s71 + c71*dth72*s72)) - c71*(dth71 + c72*dth73)*s73; c71*c73*(c72*dth71 + dth73) - s71*(c73*dth72*s72 + (dth71 + c72*dth73)*s73); -(c72*c73*dth72) + dth73*s72*s73];
    rotmpDotcol3 = [ c71*c72*dth72 - dth71*s71*s72; c72*dth72*s71 + c71*dth71*s72; -(dth72*s72)];
    mpcom = [q(n-4);q(n-3);q(n-5)];
    mpvc  = [dq(n-4);dq(n-3);dq(n-5)];
    mpaux1 = rotmpcol1*(crmp/2);
    mpaux2 = rotmpcol3*(sqrt(3)*crmp/2);
    mpvaux1 = rotmpDotcol1*(crmp/2);
    mpvaux2 = rotmpDotcol3*(sqrt(3)*crmp/2);

    posa1 = mpcom + 2*mpaux1;
    posa2 = mpcom - mpaux1 + mpaux2;
    posa3 = mpcom - mpaux1 - mpaux2;  

    vela1 = mpvc + 2*mpvaux1;
    vela2 = mpvc - mpvaux1 + mpvaux2;
    vela3 = mpvc - mpvaux1 - mpvaux2;
else
    posa1 = mpa1;
    posa2 = mpa2;
    posa3 = mpa3;

    vela1 = mpv1;
    vela2 = mpv2;
    vela3 = mpv3;
end

% Velocities of the points at which the sub-systems are joined with the moving platform

vonus1 = Qf(:,:,nus1)*tb(:,nus1);
vonus2 = Qf(:,:,nus2)*tb(:,nus2);
vonus3 = Qf(:,:,nus3)*tb(:,nus3);
vonus4 = Qf(:,:,nus4)*tb(:,nus4);
vonus5 = Qf(:,:,nus5)*tb(:,nus5);
vonus6 = Qf(:,:,nus6)*tb(:,nus6);

% sub-systems 1 and 7
eta1 = so(:,nus1)-posa1;
deta1 = vonus1-vela1;

% sub-systems 2 and 7
eta2 = so(:,nus2)-posa1;
deta2 = vonus2-vela1;

% sub-systems 3 and 7
eta3 = so(:,nus3)-posa2;
deta3 = vonus3-vela2;

% sub-systems 4 and 7
eta4 = so(:,nus4)-posa2;
deta4 = vonus4-vela2;

% sub-systems 5 and 7
eta5 = so(:,nus5)-posa3;
deta5 = vonus5-vela3;

% sub-systems 6 and 7
eta6 = so(:,nus6)-posa3;
deta6 = vonus6-vela3;

% Position constraints
eta=[eta1;eta2;eta3;eta4;eta5;eta6];

% Velocity constraints
deta=[deta1;deta2;deta3;deta4;deta5;deta6];

end

