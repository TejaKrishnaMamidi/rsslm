% RSSLM-CDPR-Type-II-DynamicIp acceleration module. This module computes the Lagrange multipliers and the joint accelerations of a CDPR

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to ddq_tree_eff

% System: 4-4 CDPR with cables attached to quadcopters

function [qdd] = acceleration(q, dq, tu, tau_d)

% Global variables -- required
global alp a b th r dx dy dz m Icxx Icyy Iczz Icxy Icyz Iczx J dJ ni; 
global nls1 nls2 nls3 nls4 nls5 nls6 nls7 nls8 nls9 nus1 nus2 nus3 nus4 nus5 nus6 nus7 nus8 nus9 invmp invmq1 invmq2 invmq3 invmq4;

% Global variables -- defined
global lam;

% Initialisation
phihat = zeros(nus9,1);
ihatTr = zeros(nus9,24);
qdd = zeros(nus9,1);
ibar = zeros(24,24);
psibar = zeros(24,1);
lam = zeros(24,1);

phil = tau_d - tu;

%% Generalised inertia matrices of fifth to ninth sub-systems and their inverses.

I5 = GIM_tree_sub_system(nls5,nus5,q,alp,a,b,th,r,dx,dy,dz,m,Icxx,Icyy,Iczz,Icxy,Icyz,Iczx);
I52 = inv3by3(I5(4:6,4:6));

I6 = GIM_tree_sub_system(nls6,nus6,q,alp,a,b,th,r,dx,dy,dz,m,Icxx,Icyy,Iczz,Icxy,Icyz,Iczx);
I62 = inv3by3(I6(4:6,4:6));

I7 = GIM_tree_sub_system(nls7,nus7,q,alp,a,b,th,r,dx,dy,dz,m,Icxx,Icyy,Iczz,Icxy,Icyz,Iczx);
I72 = inv3by3(I7(4:6,4:6));

I8 = GIM_tree_sub_system(nls8,nus8,q,alp,a,b,th,r,dx,dy,dz,m,Icxx,Icyy,Iczz,Icxy,Icyz,Iczx);
I82 = inv3by3(I8(4:6,4:6));

I9 = GIM_tree_sub_system(nls9,nus9,q,alp,a,b,th,r,dx,dy,dz,m,Icxx,Icyy,Iczz,Icxy,Icyz,Iczx);
I92 = inv3by3(I9(4:6,4:6));

%% Instead of directly inverting the sub-system inertia matrix, used reverse Gaussian elimination technique. 

% phihat(nls1:nus1) = I(nls1:nus1,nls1:nus1)\phil(nls1:nus1);
% phihat(nls2:nus2) = I(nls2:nus2,nls2:nus2)\phil(nls2:nus2);
% phihat(nls3:nus3) = I(nls3:nus3,nls3:nus3)\phil(nls3:nus3);
% phihat(nls4:nus4) = I(nls4:nus4,nls4:nus4)\phil(nls4:nus4);
% phihat(nls5:nus5) = I(nls5:nus5,nls5:nus5)\phil(nls5:nus5);
% phihat(nls6:nus6) = I(nls6:nus6,nls6:nus6)\phil(nls6:nus6);
% phihat(nls7:nus7) = I(nls7:nus7,nls7:nus7)\phil(nls7:nus7);
% phihat(nls8:nus8) = I(nls8:nus8,nls8:nus8)\phil(nls8:nus8);
% phihat(nls9:nus9) = I(nls9:nus9,nls9:nus9)\phil(nls9:nus9);

% Sub-systems 1 to 4

[ps1,ths1,bs1,hIs1,hFs1,hGs1,hetats1,hetabs1]=ddq_tree_eff_common(q(nls1:nus1), ni(1), b(nls1:nus1), th(nls1:nus1), r(nls1:nus1), dx(nls1:nus1), dy(nls1:nus1), dz(nls1:nus1), m(nls1:nus1), Icxx(nls1:nus1), Icyy(nls1:nus1), Iczz(nls1:nus1), Icxy(nls1:nus1), Icyz(nls1:nus1), Iczx(nls1:nus1));
[ps2,ths2,bs2,hIs2,hFs2,hGs2,hetats2,hetabs2]=ddq_tree_eff_common(q(nls2:nus2), ni(2), b(nls2:nus2), th(nls2:nus2), r(nls2:nus2), dx(nls2:nus2), dy(nls2:nus2), dz(nls2:nus2), m(nls2:nus2), Icxx(nls2:nus2), Icyy(nls2:nus2), Iczz(nls2:nus2), Icxy(nls2:nus2), Icyz(nls2:nus2), Iczx(nls2:nus2));
[ps3,ths3,bs3,hIs3,hFs3,hGs3,hetats3,hetabs3]=ddq_tree_eff_common(q(nls3:nus3), ni(3), b(nls3:nus3), th(nls3:nus3), r(nls3:nus3), dx(nls3:nus3), dy(nls3:nus3), dz(nls3:nus3), m(nls3:nus3), Icxx(nls3:nus3), Icyy(nls3:nus3), Iczz(nls3:nus3), Icxy(nls3:nus3), Icyz(nls3:nus3), Iczx(nls3:nus3));
[ps4,ths4,bs4,hIs4,hFs4,hGs4,hetats4,hetabs4]=ddq_tree_eff_common(q(nls4:nus4), ni(4), b(nls4:nus4), th(nls4:nus4), r(nls4:nus4), dx(nls4:nus4), dy(nls4:nus4), dz(nls4:nus4), m(nls4:nus4), Icxx(nls4:nus4), Icyy(nls4:nus4), Iczz(nls4:nus4), Icxy(nls4:nus4), Icyz(nls4:nus4), Iczx(nls4:nus4));
%
phihat(nls1:nus1) = ddq_tree_eff_variable(ni(1), alp(nls1:nus1), a(nls1:nus1),[0;cumsum(ones(ni(1)-1,1))], r(nls1:nus1), ps1, ths1, bs1, hIs1, hFs1, hGs1, hetats1, hetabs1, phil(nls1:nus1));
phihat(nls2:nus2) = ddq_tree_eff_variable(ni(2), alp(nls2:nus2), a(nls2:nus2),[0;cumsum(ones(ni(2)-1,1))], r(nls2:nus2), ps2, ths2, bs2, hIs2, hFs2, hGs2, hetats2, hetabs2, phil(nls2:nus2));
phihat(nls3:nus3) = ddq_tree_eff_variable(ni(3), alp(nls3:nus3), a(nls3:nus3),[0;cumsum(ones(ni(3)-1,1))], r(nls3:nus3), ps3, ths3, bs3, hIs3, hFs3, hGs3, hetats3, hetabs3, phil(nls3:nus3));
phihat(nls4:nus4) = ddq_tree_eff_variable(ni(4), alp(nls4:nus4), a(nls4:nus4),[0;cumsum(ones(ni(4)-1,1))], r(nls4:nus4), ps4, ths4, bs4, hIs4, hFs4, hGs4, hetats4, hetabs4, phil(nls4:nus4));

% Sub-systems 5 to 9

phihat(nls5:nls5+2) = invmp.*phil(nls5:nls5+2);
phihat(nls5+3:nus5) = I52*phil(nls5+3:nus5);

phihat(nls6:nls6+2) = invmq1.*phil(nls6:nls6+2);
phihat(nls6+3:nus6) = I62*phil(nls6+3:nus6);

phihat(nls7:nls7+2) = invmq2.*phil(nls7:nls7+2);
phihat(nls7+3:nus7) = I72*phil(nls7+3:nus7);

phihat(nls8:nls8+2) = invmq3.*phil(nls8:nls8+2);
phihat(nls8+3:nus8) = I82*phil(nls8+3:nus8);

phihat(nls9:nls9+2) = invmq4.*phil(nls9:nls9+2);
phihat(nls9+3:nus9) = I92*phil(nls9+3:nus9);

%% Optimised the code below considering the sparsity of the Jacobian matrix (4-4 CDPR).

% ihatTr(nls1:nus1,1:24) = I(nls1:nus1,nls1:nus1)\J(:,nls1:nus1)';
% ihatTr(nls2:nus2,1:24) = I(nls2:nus2,nls2:nus2)\J(:,nls2:nus2)';
% ihatTr(nls3:nus3,1:24) = I(nls3:nus3,nls3:nus3)\J(:,nls3:nus3)';
% ihatTr(nls4:nus4,1:24) = I(nls4:nus4,nls4:nus4)\J(:,nls4:nus4)';
% ihatTr(nls5:nus5,1:24) = I(nls5:nus5,nls5:nus5)\J(:,nls5:nus5)';
% ihatTr(nls6:nus6,1:24) = I(nls6:nus6,nls6:nus6)\J(:,nls6:nus6)';
% ihatTr(nls7:nus7,1:24) = I(nls7:nus7,nls7:nus7)\J(:,nls7:nus7)';
% ihatTr(nls8:nus8,1:24) = I(nls8:nus8,nls8:nus8)\J(:,nls8:nus8)';
% ihatTr(nls9:nus9,1:24) = I(nls8:nus9,nls9:nus9)\J(:,nls9:nus9)';

for ii = 1:1:3
    ihatTr(nls1:nus1,ii) = ddq_tree_eff_variable(ni(1), alp(nls1:nus1), a(nls1:nus1),[0;cumsum(ones(ni(1)-1,1))], r(nls1:nus1), ps1, ths1, bs1, hIs1, hFs1, hGs1, hetats1, hetabs1, J(ii,nls1:nus1)');
    ihatTr(nls5:nls5+2,ii) = invmp.*J(ii,nls5:nls5+2)';
    ihatTr(nls5+3:nus5,ii) = I52*J(ii,nls5+3:nus5)';
end

for ii = 4:1:6
    ihatTr(nls2:nus2,ii) = ddq_tree_eff_variable(ni(2), alp(nls2:nus2), a(nls2:nus2),[0;cumsum(ones(ni(2)-1,1))], r(nls2:nus2), ps2, ths2, bs2, hIs2, hFs2, hGs2, hetats2, hetabs2, J(ii,nls2:nus2)');
    ihatTr(nls5:nls5+2,ii) = invmp.*J(ii,nls5:nls5+2)';
    ihatTr(nls5+3:nus5,ii) = I52*J(ii,nls5+3:nus5)';
end

for ii = 7:1:9
    ihatTr(nls3:nus3,ii) = ddq_tree_eff_variable(ni(3), alp(nls3:nus3), a(nls3:nus3),[0;cumsum(ones(ni(3)-1,1))], r(nls3:nus3), ps3, ths3, bs3, hIs3, hFs3, hGs3, hetats3, hetabs3, J(ii,nls3:nus3)');
    ihatTr(nls5:nls5+2,ii) = invmp.*J(ii,nls5:nls5+2)';
    ihatTr(nls5+3:nus5,ii) = I52*J(ii,nls5+3:nus5)';
end

for ii = 10:1:12
    ihatTr(nls4:nus4,ii) = ddq_tree_eff_variable(ni(4), alp(nls4:nus4), a(nls4:nus4),[0;cumsum(ones(ni(4)-1,1))], r(nls4:nus4), ps4, ths4, bs4, hIs4, hFs4, hGs4, hetats4, hetabs4, J(ii,nls4:nus4)');
    ihatTr(nls5:nls5+2,ii) = invmp.*J(ii,nls5:nls5+2)';
    ihatTr(nls5+3:nus5,ii) = I52*J(ii,nls5+3:nus5)';
end

for ii = 13:1:15
    ihatTr(nls1:nus1,ii) = ddq_tree_eff_variable(ni(1), alp(nls1:nus1), a(nls1:nus1),[0;cumsum(ones(ni(1)-1,1))], r(nls1:nus1), ps1, ths1, bs1, hIs1, hFs1, hGs1, hetats1, hetabs1, J(ii,nls1:nus1)');
    ihatTr(nls6:nls6+2,ii) = invmq1.*J(ii,nls6:nls6+2)';
end

for ii = 16:1:18
    ihatTr(nls2:nus2,ii) = ddq_tree_eff_variable(ni(2), alp(nls2:nus2), a(nls2:nus2),[0;cumsum(ones(ni(2)-1,1))], r(nls2:nus2), ps2, ths2, bs2, hIs2, hFs2, hGs2, hetats2, hetabs2, J(ii,nls2:nus2)');
    ihatTr(nls7:nls7+2,ii) = invmq2.*J(ii,nls7:nls7+2)';
end

for ii = 19:1:21
    ihatTr(nls3:nus3,ii) = ddq_tree_eff_variable(ni(3), alp(nls3:nus3), a(nls3:nus3),[0;cumsum(ones(ni(3)-1,1))], r(nls3:nus3), ps3, ths3, bs3, hIs3, hFs3, hGs3, hetats3, hetabs3, J(ii,nls3:nus3)');
    ihatTr(nls8:nls8+2,ii) = invmq3.*J(ii,nls8:nls8+2)';
end

for ii = 22:1:24
    ihatTr(nls4:nus4,ii) = ddq_tree_eff_variable(ni(4), alp(nls4:nus4), a(nls4:nus4),[0;cumsum(ones(ni(4)-1,1))], r(nls4:nus4), ps4, ths4, bs4, hIs4, hFs4, hGs4, hetats4, hetabs4, J(ii,nls4:nus4)');
    ihatTr(nls9:nls9+2,ii) = invmq4.*J(ii,nls9:nls9+2)';
end

%% Optimised the code below considering the sparsity of the Jacobian matrix and the properties of the matrix ibar.

%ibar = J*ihatTr;

% Rows 1-3
ibar(1:3,1:3) = J(1:3,nls1:nus1)*ihatTr(nls1:nus1,1:3)+J(1:3,nls5:nus5)*ihatTr(nls5:nus5,1:3);
ibar(1:3,4:12) = J(1:3,nls5:nus5)*ihatTr(nls5:nus5,4:12);
ibar(1:3,13:15) = J(1:3,nls1:nus1)*ihatTr(nls1:nus1,13:15);

% Rows 4-6
ibar(4:6,1:3) = ibar(1:3,4:6)';
ibar(4:6,4:6) = J(4:6,nls2:nus2)*ihatTr(nls2:nus2,4:6)+ J(4:6,nls5:nus5)*ihatTr(nls5:nus5,4:6);
ibar(4:6,7:12) = J(4:6,nls5:nus5)*ihatTr(nls5:nus5,7:12);
ibar(4:6,16:18) = J(4:6,nls2:nus2)*ihatTr(nls2:nus2,16:18);

% Rows 7-9
ibar(7:9,1:6) = ibar(1:6,7:9)';
ibar(7:9,7:9) = J(7:9,nls3:nus3)*ihatTr(nls3:nus3,7:9)+ J(7:9,nls5:nus5)*ihatTr(nls5:nus5,7:9);
ibar(7:9,10:12) = J(7:9,nls5:nus5)*ihatTr(nls5:nus5,10:12);
ibar(7:9,19:21) = J(7:9,nls3:nus3)*ihatTr(nls3:nus3,19:21);

% Rows 10-12
ibar(10:12,1:9) = ibar(1:9,10:12)';
ibar(10:12,10:12) = J(10:12,nls4:nus4)*ihatTr(nls4:nus4,10:12) + J(10:12,nls5:nus5)*ihatTr(nls5:nus5, 10:12);
ibar(10:12,22:24) = J(10:12,nls4:nus4)*ihatTr(nls4:nus4,22:24);

% Rows 13-15
ibar(13:15,1:12) = ibar(1:12,13:15)';
ibar(13:15,13:15) = J(13:15,nls1:nls1+2)*ihatTr(nls1:nls1+2,13:15) + J(13:15, nls6:nls6+2)*ihatTr(nls6:nls6+2, 13:15);

% Rows 16-18
ibar(16:18,1:15) = ibar(1:15,16:18)';
ibar(16:18,16:18) = J(16:18,nls2:nls2+2)*ihatTr(nls2:nls2+2,16:18) + J(16:18, nls7:nls7+2)*ihatTr(nls7:nls7+2, 16:18);

% Rows 19-21
ibar(19:21,1:18) = ibar(1:18,19:21)';
ibar(19:21,19:21) = J(19:21,nls3:nls3+2)*ihatTr(nls3:nls3+2,19:21) + J(19:21, nls8:nls8+2)*ihatTr(nls8:nls8+2, 19:21);

% Rows 22-24
ibar(22:24,1:21) = ibar(1:21,22:24)';
ibar(22:24,22:24) = J(22:24,nls4:nls4+2)*ihatTr(nls4:nls4+2,22:24) + J(22:24, nls9:nls9+2)*ihatTr(nls9:nls9+2, 22:24);

%% Optimised the code below considering the sparsity of the Jacobian matrix and its derivative.

% psibar = -dJ*dq-J*phihat;

psibar(1:3) = -dJ(1:3,nls1:nus1)*dq(nls1:nus1)-dJ(1:3,nls5+3:nus5)*dq(nls5+3:nus5) - J(1:3,nls1:nus1)*phihat(nls1:nus1)-J(1:3,nls5:nus5)*phihat(nls5:nus5);
psibar(4:6) = -dJ(4:6,nls2:nus2)*dq(nls2:nus2)-dJ(4:6,nls5+3:nus5)*dq(nls5+3:nus5) - J(4:6,nls2:nus2)*phihat(nls2:nus2)-J(4:6,nls5:nus5)*phihat(nls5:nus5);
psibar(7:9) = -dJ(7:9,nls3:nus3)*dq(nls3:nus3)-dJ(7:9,nls5+3:nus5)*dq(nls5+3:nus5) - J(7:9,nls3:nus3)*phihat(nls3:nus3)-J(7:9,nls5:nus5)*phihat(nls5:nus5);
psibar(10:12) = -dJ(10:12,nls4:nus4)*dq(nls4:nus4)-dJ(10:12,nls5+3:nus5)*dq(nls5+3:nus5) - J(10:12,nls4:nus4)*phihat(nls4:nus4)-J(10:12,nls5:nus5)*phihat(nls5:nus5);
psibar(13:15) = -J(13:15,nls1:nls1+2)*phihat(nls1:nls1+2)-J(13:15,nls6:nls6+2)*phihat(nls6:nls6+2);
psibar(16:18) = -J(16:18,nls2:nls2+2)*phihat(nls2:nls2+2)-J(16:18,nls7:nls7+2)*phihat(nls7:nls7+2);
psibar(19:21) = -J(19:21,nls3:nls3+2)*phihat(nls3:nls3+2)-J(19:21,nls8:nls8+2)*phihat(nls8:nls8+2);
psibar(22:24) = -J(22:24,nls4:nls4+2)*phihat(nls4:nls4+2)-J(22:24,nls9:nls9+2)*phihat(nls9:nls9+2);

%% Optimised the code below by performing block Gaussian elimination.
% lam = ibar\psibar;

auxa  = inv3by3(ibar(1:3,1:3));
auxa1 = ibar(4:6,1:3)*auxa;
auxa2 = ibar(7:9,1:3)*auxa;
auxa3 = ibar(10:12,1:3)*auxa;
auxa4 = ibar(13:15,1:3)*auxa;
%
auxb  = inv3by3(ibar(4:6,4:6)-auxa1*ibar(1:3,4:6));
auxb1 = (ibar(7:9,4:6)-auxa2*ibar(1:3,4:6))*auxb;
auxb2 = (ibar(10:12,4:6)-auxa3*ibar(1:3,4:6))*auxb;
auxb3 = (-auxa4*ibar(1:3,4:6))*auxb;
auxb4 = ibar(16:18,4:6)*auxb;
%
auxc0 = ibar(4:6,7:9)-auxa1*ibar(1:3,7:9);
auxc  = inv3by3(ibar(7:9,7:9)-auxa2*ibar(1:3,7:9)-auxb1*auxc0);
auxc1 = (ibar(10:12,7:9)-auxa3*ibar(1:3,7:9)-auxb2*auxc0)*auxc;
auxc2 = (-auxa4*ibar(1:3,7:9)-auxb3*auxc0)*auxc;
auxc3 = (-auxb4*auxc0)*auxc;
auxc4 = ibar(19:21,7:9)*auxc;
%
auxd01 = ibar(4:6,10:12)-auxa1*ibar(1:3,10:12);
auxd02 = ibar(7:9,10:12)-auxa2*ibar(1:3,10:12)-auxb1*auxd01; 
auxd = inv3by3(ibar(10:12,10:12)-auxa3*ibar(1:3,10:12)-auxb2*auxd01-auxc1*auxd02);
auxd1 = (-auxa4*ibar(1:3,10:12)-auxb3*auxd01-auxc2*auxd02)*auxd;
auxd2 = (-auxb4*auxd01-auxc3*auxd02)*auxd;
auxd3 = (-auxc4*auxd02)*auxd;
auxd4 = ibar(22:24,10:12)*auxd;
%
auxe01 = -auxa1*ibar(1:3,13:15);
auxe02 = -auxa2*ibar(1:3,13:15)-auxb1*auxe01;
auxe03 = -auxa3*ibar(1:3,13:15)-auxb2*auxe01-auxc1*auxe02;
auxe = inv3by3(ibar(13:15,13:15)-auxa4*ibar(1:3,13:15)-auxb3*auxe01-auxc2*auxe02-auxd1*auxe03);
auxe1 = (-auxb4*auxe01-auxc3*auxe02-auxd2*auxe03)*auxe;
auxe2 = (-auxc4*auxe02-auxd3*auxe03)*auxe;
auxe3 = (-auxd4*auxe03)*auxe;
%
auxf01 = ibar(4:6,16:18);
auxf02 = -auxb1*auxf01;
auxf03 = -auxb2*auxf01-auxc1*auxf02;
auxf04 = -auxb3*auxf01-auxc2*auxf02-auxd1*auxf03;
auxf = inv3by3(ibar(16:18,16:18)-auxb4*auxf01-auxc3*auxf02-auxd2*auxf03-auxe1*auxf04);
auxf1 = (-auxc4*auxf02-auxd3*auxf03-auxe2*auxf04)*auxf;
auxf2 = (-auxd4*auxf03-auxe3*auxf04)*auxf;
%
auxg01 = ibar(7:9,19:21);
auxg02 = -auxc1*auxg01;
auxg03 = -auxc2*auxg01-auxd1*auxg02;
auxg04 = -auxc3*auxg01-auxd2*auxg02-auxe1*auxg03;
auxg = inv3by3(ibar(19:21,19:21)-auxc4*auxg01-auxd3*auxg02-auxe2*auxg03-auxf1*auxg04);
auxg1 = (-auxd4*auxg02-auxe3*auxg03-auxf2*auxg04)*auxg;
%
auxh01 = ibar(10:12,22:24);
auxh02 = -auxd1*auxh01;
auxh03 = -auxd2*auxh01-auxe1*auxh02;
auxh04 = -auxd3*auxh01-auxe2*auxh02-auxf1*auxh03;
auxh = inv3by3(ibar(22:24,22:24)-auxd4*auxh01-auxe3*auxh02-auxf2*auxh03-auxg1*auxh04);
%
auxp1 = psibar(4:6)-auxa1*psibar(1:3);
auxp2 = psibar(7:9)-auxa2*psibar(1:3)-auxb1*auxp1;
auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
auxp5 = psibar(16:18) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
auxp6 = psibar(19:21) - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
auxp7 = psibar(22:24) - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;

% Back substitution

lam(22:24) = auxh*auxp7;
lam(19:21) = auxg*(auxp6-auxh04*lam(22:24));
lam(16:18) = auxf*(auxp5-auxh03*lam(22:24)-auxg04*lam(19:21));
lam(13:15) = auxe*(auxp4-auxh02*lam(22:24)-auxg03*lam(19:21)-auxf04*lam(16:18));
lam(10:12) = auxd*(auxp3-auxh01*lam(22:24)-auxg02*lam(19:21)-auxf03*lam(16:18)-auxe03*lam(13:15));
lam(7:9)   = auxc*(auxp2-auxg01*lam(19:21)-auxf02*lam(16:18)-auxe02*lam(13:15)-auxd02*lam(10:12));
lam(4:6)   = auxb*(auxp1-auxf01*lam(16:18)-auxe01*lam(13:15)-auxd01*lam(10:12)-auxc0*lam(7:9));
lam(1:3)   = auxa*(psibar(1:3)-ibar(1:3,4:6)*lam(4:6)-ibar(1:3,7:9)*lam(7:9)-ibar(1:3,10:12)*lam(10:12)-ibar(1:3,13:15)*lam(13:15));

%% Joint accelerations

% Optimised the code below using the sparsity of ihatTr
qdd(nls1:nus1) = phihat(nls1:nus1)+ihatTr(nls1:nus1,1:3)*lam(1:3)+ihatTr(nls1:nus1,13:15)*lam(13:15);
qdd(nls2:nus2) = phihat(nls2:nus2)+ihatTr(nls2:nus2,4:6)*lam(4:6)+ihatTr(nls2:nus2,16:18)*lam(16:18);
qdd(nls3:nus3) = phihat(nls3:nus3)+ihatTr(nls3:nus3,7:9)*lam(7:9)+ihatTr(nls3:nus3,19:21)*lam(19:21);
qdd(nls4:nus4) = phihat(nls4:nus4)+ihatTr(nls4:nus4,10:12)*lam(10:12)+ihatTr(nls4:nus4,22:24)*lam(22:24);
qdd(nls5:nus5) = phihat(nls5:nus5)+ihatTr(nls5:nus5,1:12)*lam(1:12);
qdd(nls6:nus6) = phihat(nls6:nus6)+[ihatTr(nls6:nls6+2,13:15)*lam(13:15); zeros(nus6-(nls6+2),1)];
qdd(nls7:nus7) = phihat(nls7:nus7)+[ihatTr(nls7:nls7+2,16:18)*lam(16:18); zeros(nus7-(nls7+2),1)];
qdd(nls8:nus8) = phihat(nls8:nus8)+[ihatTr(nls8:nls8+2,19:21)*lam(19:21); zeros(nus8-(nls8+2),1)];
qdd(nls9:nus9) = phihat(nls9:nus9)+[ihatTr(nls9:nls9+2,22:24)*lam(22:24); zeros(nus9-(nls9+2),1)];

%qdd(nls1:nus1) = phihat(nls1:nus1)+ihatTr(nls1:nus1,1:24)*lam;
%qdd(nls2:nus2) = phihat(nls2:nus2)+ihatTr(nls2:nus2,1:24)*lam;
%qdd(nls3:nus3) = phihat(nls3:nus3)+ihatTr(nls3:nus3,1:24)*lam;
%qdd(nls4:nus4) = phihat(nls4:nus4)+ihatTr(nls4:nus4,1:24)*lam;
%qdd(nls5:nus5) = phihat(nls5:nus5)+ihatTr(nls5:nus5,1:24)*lam;
%qdd(nls6:nus6) = phihat(nls6:nus6)+ihatTr(nls6:nus6,1:24)*lam;
%qdd(nls7:nus7) = phihat(nls7:nus7)+ihatTr(nls7:nus7,1:24)*lam;
%qdd(nls8:nus8) = phihat(nls8:nus8)+ihatTr(nls8:nus8,1:24)*lam;
%qdd(nls9:nus9) = phihat(nls9:nus9)+ihatTr(nls9:nus9,1:24)*lam;

end
