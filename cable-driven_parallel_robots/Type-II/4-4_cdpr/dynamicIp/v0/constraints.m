% RSSLM-CDPR-Type-II-DynamicIp constraints module. The values of constraint functions and their derivatives are computed.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

function [eta, deta]=constraints(il, iu, q, dq, tt, tb, so, st, vt, Qf, p, b, th)

% No function calls

% System: 4-4 CDPR with cables attached to quadcopters

% Global variables -- required
global alp a bt r dx dy dz nus1 nus2 nus3 nus4 nus5 mpa mpv;

% Computation of the change in the constraint functions
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
            -b(ii)*salp
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

% Only if the variables q(nls5) to q(nus5) are varied the following computations are done.  

if il>nus4
    c73 = cos(q(nus5)); s73=sin(q(nus5));
    c72 = cos(q(nus5-1)-pi/2); s72=sin(q(nus5-1)-pi/2);
    c71 = cos(q(nus5-2)); s71=sin(q(nus5-2));
    dth73 = dq(nus5); dth72 = dq(nus5-1); dth71 = dq(nus5-2);
    
    rotmp = [c72*c73,              -c72*s73,                   s72
             c73*s71*s72 + c71*s73,  c71*c73 - s71*s72*s73,    -s71*c72
            -(c71*c73*s72)+s71*s73, c73*s71+c71*s72*s73,       c71*c72];

    rotmpDot = [-c73*s72*dth72-c72*s73*dth73, s72*s73*dth72-c72*c73*dth73, c72*dth72
                c71*c73*s72*dth71 - s71*s73*dth71 + c72*c73*s71*dth72 + c71*c73*dth73 - s71*s72*s73*dth73, -c73*s71*dth71 -  c71*s72*s73*dth71 - c72*s71*s73*dth72 -  c73*s71*s72*dth73 - c71*s73*dth73, -c71*c72*dth71 + s71*s72*dth72
                c73*s71*s72*dth71 + c71*s73*dth71 - c71*c72*c73*dth72 + c73*s71*dth73 + c71*s72*s73*dth73, c71*c73*dth71 -  s71*s72*s73*dth71 + c71*c72*s73*dth72 + c71*c73*s72*dth73 - s71*s73*dth73, -c72*s71*dth71 - c71*s72*dth72];

    mpcom = [q(nus5-4);q(nus5-3);q(nus5-5)];
    mpvc  = [dq(nus5-4);dq(nus5-3);dq(nus5-5)];

    posa1 = mpcom + rotmp*[-3/10; -2/5; 1/10];
    posa2 = mpcom + rotmp*[3/10; -2/5; 1/10];
    posa3 = mpcom + rotmp*[-3/10; 2/5; 1/10]; 
    posa4 = mpcom + rotmp*[3/10; 2/5; 1/10]; 

    vela1 = mpvc + rotmpDot*[-3/10; -2/5; 1/10];
    vela2 = mpvc + rotmpDot*[3/10; -2/5; 1/10];
    vela3 = mpvc + rotmpDot*[-3/10; 2/5; 1/10];
    vela4 = mpvc + rotmpDot*[3/10; 2/5; 1/10];
else
    posa1 = mpa(:,1);
    posa2 = mpa(:,2);
    posa3 = mpa(:,3);
    posa4 = mpa(:,4);

    vela1 = mpv(:,1);
    vela2 = mpv(:,2);
    vela3 = mpv(:,3);
    vela4 = mpv(:,4);
end

% Velocities of the points at which the sub-systems are joined with the moving platform

vonus1 = Qf(:,:,nus1)*tb(:,nus1);
vonus2 = Qf(:,:,nus2)*tb(:,nus2);
vonus3 = Qf(:,:,nus3)*tb(:,nus3);
vonus4 = Qf(:,:,nus4)*tb(:,nus4);

% sub-systems 1 and 5
eta1 = so(:,nus1)-posa1;
deta1 = vonus1-vela1;

% sub-systems 2 and 5
eta2 = so(:,nus2)-posa2;
deta2 = vonus2-vela2;

% sub-systems 3 and 5
eta3 = so(:,nus3)-posa3;
deta3 = vonus3-vela3;

% sub-systems 4 and 5
eta4 = so(:,nus4)-posa4;
deta4 = vonus4-vela4;

% Position constraints
eta=[eta1;eta2;eta3;eta4];

% Velocity constraints
deta=[deta1;deta2;deta3;deta4];

end

