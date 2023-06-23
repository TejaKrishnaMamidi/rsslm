% RSSLM-CDPR-Type-I (or Kinetostatics) for_kine module. It solves the forward kinematics of tree-type systems.

% This code was developed by Dr. Suril V. Shah
% Updated latest : 10 Aug 2013 -- minor modifications were done after that to suit the needs.
% This program calculate forward kinematics in the inertial frame

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function call to correctionYcoordinate

function [so, sc, vc, tt, st]=for_kine(q, dq, n, alp, a, b, th, bt, r, dx, dy, dz)

% FORWARD RECURSION _FINDING TWIST AND TWIST RATE
%Initialization
e=[0;0;1];
tt=zeros(3,n);
tb=zeros(3,n);
so=zeros(3,n);
sc=zeros(3,n);
st=zeros(3,n);
vc=zeros(3,n);
Qf=zeros(3,3,n);
p=zeros(n,1);

% FOR LOOP STARTS
for i=1:n
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
        di=[dx(i);dy(i);dz(i)];
        aim=[a(i)
            correctionYcoordinate(i)-b(i)*sin(alp(i))
            b(i)*cos(alp(i))];
        so(:,i)=aim;
        sc(:,i)=so(:,i)+Qf(:,:,i)*di;
        st(:,i)=so(:,i)+Qf(:,:,i)*(2*di);
        %w angular velocity
        edq=e*dq(i);
        tt(:,i)=r(i)*edq;
        tti=tt(:,i);
        
        %v linear velocity
        tb(:,i)=p(i)*edq;
        ttixdi=[tti(2)*di(3)-di(2)*tti(3);-(tti(1)*di(3)-di(1)*tti(3));tti(1)*di(2)-di(1)*tti(2)];
        vc(:,i)=Qf(:,:,i)*(tb(:,i)+ttixdi);
        
    else %Calculation for the links other than those attached with ground
        Qi=[cth,      -sth,       0
         calp*sth,  calp*cth, -salp
         salp*sth,  salp*cth,  calp];
        Qf(:,:,i)=Qf(:,:,bt(i))*Qi;
        
        %position vector from origin of link to origin of next link
        aim=[a(i)
            -b(i)*sin(alp(i))
             b(i)*cos(alp(i))];
        di=[dx(i);dy(i);dz(i)];
        
        %Positions
        so(:,i)=so(:,bt(i))+Qf(:,:,bt(i))*aim;
        sc(:,i)=so(:,i)+Qf(:,:,i)*di;
        st(:,i)=so(:,i)+Qf(:,:,i)*(2*di);
        
        %w angular velocity
        ttbi=tt(:,bt(i));
        edq=e*dq(i);
        tt(:,i)=Qi'*ttbi+r(i)*edq;
        tti=tt(:,i);
        
        %v  linear velocity
        ttbixaim=[ttbi(2)*aim(3)-aim(2)*ttbi(3);-(ttbi(1)*aim(3)-aim(1)*ttbi(3));ttbi(1)*aim(2)-aim(1)*ttbi(2)];
        tb(:,i)=Qi.'*(tb(:,bt(i))+ttbixaim)+p(i)*edq;
        ttixdi=[tti(2)*di(3)-di(2)*tti(3);-(tti(1)*di(3)-di(1)*tti(3));tti(1)*di(2)-di(1)*tti(2)];
        vc(:,i)=Qf(:,:,i)*(tb(:,i)+ttixdi);
    end
end
end
