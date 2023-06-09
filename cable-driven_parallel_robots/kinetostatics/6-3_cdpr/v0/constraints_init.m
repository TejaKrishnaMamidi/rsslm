% RSSLM-CDPR-Kinetostatics constraints module. The variables that do not change with respect to a configuration variable while computing the values of constraint functions and their derivatives are initialised here.

% Contibutors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

function [th, b]=constraints_init(th, b)

% No function calls

% System: 6-3 CDPR

% Global variables -- required
global n alp a bt r dx dy dz tt tb so st vt Qf p q dq crmp;

% Global variables -- defined
global mpa1 mpa2 mpa3 mpv1 mpv2 mpv3;

%% Initialisation
tt=zeros(3,n);
tb=zeros(3,n);
so=zeros(3,n);
st=zeros(3,n);
vt=zeros(3,n);
Qf=zeros(3,3,n);
p=zeros(n,1);

%% Positions and velocities of the distal points on the links of cables
for i=1:(n-6)
    p(i)=1-r(i);
    th(i)=th(i)*p(i)+q(i)*r(i);
    b(i)=b(i)*r(i)+q(i)*p(i);
    cth=cos(th(i)); calp=cos(alp(i));
    sth=sin(th(i)); salp=sin(alp(i));
    if bt(i)==0 %When parent of the link is ground link
        Qi=[cth,      -sth,       0
            calp*sth,  calp*cth, -salp
            salp*sth,  salp*cth,  calp];
        Qf(:,:,i)=Qi;
        
        %Positions
        di=[dx(i);dy(i);dz(i)]+[dx(i);dy(i);dz(i)];
        aim=[a(i)
            -b(i)*salp
            b(i)*calp];
        so(:,i)=aim;
        st(:,i)=so(:,i)+Qf(:,:,i)*di;
        
        %w angular velocity
        tt(:,i)=[0;0;r(i)*dq(i)];
        tti=tt(:,i);
        
        %v linear velocity
        tb(:,i)=[0;0;p(i)*dq(i)];
        vt(:,i)=[0;-salp*p(i)*dq(i);calp*p(i)*dq(i)]-Qi(:,1)*di(2)*tti(3)+Qi(:,2)*di(1)*tti(3);
        
    else %Calculation for the links other than those attached with ground
        Qi=[cth,      -sth,       0
            calp*sth,  calp*cth, -salp
            salp*sth,  salp*cth,  calp];
        Qf(:,:,i)=Qf(:,:,bt(i))*Qi;
        
        % position vector from origin of link to origin of next link
        aim=[a(i)
            -b(i)*salp
            b(i)*calp];
        di=[dx(i);dy(i);dz(i)]+[dx(i);dy(i);dz(i)];
        
        %Positions
        so(:,i)=so(:,bt(i))+Qf(:,:,bt(i))*aim;
        st(:,i)=so(:,i)+Qf(:,:,i)*di;
        
        %w angular velocity
        ttbi=tt(:,bt(i));
        tt(:,i)=Qi'*ttbi+[0;0;r(i)*dq(i)];
        tti=tt(:,i);
        %v  linear velocity
        ttbixaim=[ttbi(2)*aim(3)-aim(2)*ttbi(3);-(ttbi(1)*aim(3)-aim(1)*ttbi(3));ttbi(1)*aim(2)-aim(1)*ttbi(2)];
        tb(:,i)=Qi.'*(tb(:,bt(i))+ttbixaim)+[0;0;p(i)*dq(i)];
        ttixdi=[tti(2)*di(3)-di(2)*tti(3);-(tti(1)*di(3)-di(1)*tti(3));tti(1)*di(2)-di(1)*tti(2)];
        vt(:,i)=Qf(:,:,i)*(tb(:,i)+ttixdi);
    end
end

% Rotation matrix representing the orientation of the moving platform.
c73 = cos(q(n)); s73=sin(q(n));
c72 = cos(q(n-1)); s72=sin(q(n-1));
c71 = cos(q(n-2)); s71=sin(q(n-2));

dth73 = dq(n); dth72 = dq(n-1); dth71 = dq(n-2);

%rotmp = [c71*c72*c73 - s71*s73, -(c73*s71) - c71*c72*s73, c71*s72
%	c72*c73*s71 + c71*s73,  c71*c73 - c72*s71*s73,    s71*s72
%	-(c73*s72),             s72*s73,                  c72];
rotmpcol1 = [c71*c72*c73 - s71*s73; c72*c73*s71 + c71*s73; -(c73*s72)];
rotmpcol3 = [c71*s72; s71*s72; c72];

%rotmpDot = [-(c73*(c72*dth71*s71 + dth73*s71 + c71*dth72*s72)) - c71*(dth71 + c72*dth73)*s73, -(c71*c73*(dth71 + c72*dth73)) + (c72*dth71 + dth73)*s71*s73 + c71*dth72*s72*s73, c71*c72*dth72 - dth71*s71*s72
%	    c71*c73*(c72*dth71 + dth73) - s71*(c73*dth72*s72 + (dth71 + c72*dth73)*s73), -(c73*(dth71 + c72*dth73)*s71) - c71*(c72*dth71 + dth73)*s73 + dth72*s71*s72*s73, c72*dth72*s71 + c71*dth71*s72
%	    -(c72*c73*dth72) + dth73*s72*s73, c73*dth73*s72 + c72*dth72*s73,-(dth72*s72)];

rotmpDotcol1 = [-(c73*(c72*dth71*s71 + dth73*s71 + c71*dth72*s72)) - c71*(dth71 + c72*dth73)*s73; c71*c73*(c72*dth71 + dth73) - s71*(c73*dth72*s72 + (dth71 + c72*dth73)*s73); -(c72*c73*dth72) + dth73*s72*s73];
rotmpDotcol3 = [ c71*c72*dth72 - dth71*s71*s72; c72*dth72*s71 + c71*dth71*s72; -(dth72*s72)];

% Positions and velocities of the vertices of the moving platform

mpcom = [q(n-4);q(n-3);q(n-5)];
mpvc  = [dq(n-4);dq(n-3);dq(n-5)];
mpaux1 = rotmpcol1*(crmp/2);
mpaux2 = rotmpcol3*(sqrt(3)*crmp/2);
mpvaux1 = rotmpDotcol1*(crmp/2);
mpvaux2 = rotmpDotcol3*(sqrt(3)*crmp/2);

mpa1 = mpcom + 2*mpaux1;
mpa2 = mpcom - mpaux1 + mpaux2;
mpa3 = mpcom - mpaux1 - mpaux2;  

mpv1 = mpvc + 2*mpvaux1;
mpv2 = mpvc - mpvaux1 + mpvaux2;
mpv3 = mpvc - mpvaux1 - mpvaux2;

end

