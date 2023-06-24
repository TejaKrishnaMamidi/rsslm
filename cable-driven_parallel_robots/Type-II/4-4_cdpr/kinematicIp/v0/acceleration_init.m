% RSSLM-CDPR-Type-II-KinematicIp acceleration_init module. The common portions in determining the changes in joint accelerations with perturbations in the joint positions and velocities are computed here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to ddq_tree_eff

% System: 4-4 CDPR with movements of the cables' exit points.

function [phihat, ihatTr, ibar, psibar, psibarAux1, psibarAux2, psibarAux3, auxa, auxa1, auxa2, auxa3, auxa4, auxb, auxb1, auxb2, auxb3, auxb4, auxc0, auxc, auxc1, auxc2, auxc3, auxc4, auxd01, auxd02, auxd, auxd1, auxd2, auxd3, auxd4, auxe01, auxe02, auxe03, auxe, auxe1, auxe2, auxe3, auxf01, auxf02, auxf03, auxf04, auxf, auxf1, auxf2, auxg01, auxg02, auxg03, auxg04, auxg, auxg1, auxh01, auxh02, auxh03, auxh04, auxh, auxp1, auxp2, auxp3, auxp4, auxp5, auxp6, auxp7] = acceleration_init(q, dq, tu, tau_d)

% Global variables -- required
global alp a b th r dx dy dz m Icxx Icyy Iczz Icxy Icyz Iczx J dJ ni;
global nls1 nls2 nls3 nls4 nls5 nus1 nus2 nus3 nus4 nus5 invmp n nc bitrj;

% Initialisation
phihat = zeros(n,1);
ihatTr = zeros(n,nc);
ibar = zeros(nc,nc);
psibarAux1 = zeros(nc,1);
psibarAux2 = zeros(nc,1);
psibarAux3 = zeros(nc,1);

phil = tau_d - tu;

%% Generalised inertia matrices of the fifth sub-system and its inverse.

I5 = GIM_tree_sub_system(nls5,nus5,q,alp,a,b,th,r,dx,dy,dz,m,Icxx,Icyy,Iczz,Icxy,Icyz,Iczx);
I52 = inv3by3(I5(4:6,4:6));

% Instead of directly inverting the sub-system inertia matrix, used reverse Gaussian elimination technique. 

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

% Sub-system 5

phihat(nls5:nls5+2) = invmp.*phil(nls5:nls5+2);
phihat(nls5+3:nus5) = I52*phil(nls5+3:nus5);

%% Optimised the code below considering the sparsity of the Jacobian matrix (4-4 CDPR).

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
end

for ii = 16:1:18
    ihatTr(nls2:nus2,ii) = ddq_tree_eff_variable(ni(2), alp(nls2:nus2), a(nls2:nus2),[0;cumsum(ones(ni(2)-1,1))], r(nls2:nus2), ps2, ths2, bs2, hIs2, hFs2, hGs2, hetats2, hetabs2, J(ii,nls2:nus2)');
end

for ii = 19:1:21
    ihatTr(nls3:nus3,ii) = ddq_tree_eff_variable(ni(3), alp(nls3:nus3), a(nls3:nus3),[0;cumsum(ones(ni(3)-1,1))], r(nls3:nus3), ps3, ths3, bs3, hIs3, hFs3, hGs3, hetats3, hetabs3, J(ii,nls3:nus3)');
end

for ii = 22:1:24
    ihatTr(nls4:nus4,ii) = ddq_tree_eff_variable(ni(4), alp(nls4:nus4), a(nls4:nus4),[0;cumsum(ones(ni(4)-1,1))], r(nls4:nus4), ps4, ths4, bs4, hIs4, hFs4, hGs4, hetats4, hetabs4, J(ii,nls4:nus4)');
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
ibar(7:9,7:9) = J(7:9,nls3:nus3)*ihatTr(nls3:nus3,7:9)+ J(7:9, nls5:nus5)*ihatTr(nls5:nus5,7:9);
ibar(7:9,10:12) = J(7:9,nls5:nus5)*ihatTr(nls5:nus5,10:12);
ibar(7:9,19:21) = J(7:9,nls3:nus3)*ihatTr(nls3:nus3,19:21);

% Rows 10-12
ibar(10:12,1:9) = ibar(1:9,10:12)';
ibar(10:12,10:12) = J(10:12,nls4:nus4)*ihatTr(nls4:nus4,10:12) + J(10:12,nls5:nus5)*ihatTr(nls5:nus5, 10:12);
ibar(10:12,22:24) = J(10:12,nls4:nus4)*ihatTr(nls4:nus4,22:24);

% Rows 13-15
ibar(13:15,1:12) = ibar(1:12,13:15)';
ibar(13:15,13:15) = J(13:15,nls1:nls1+2)*ihatTr(nls1:nls1+2,13:15);

% Rows 16-18
ibar(16:18,1:15) = ibar(1:15,16:18)';
ibar(16:18,16:18) = J(16:18,nls2:nls2+2)*ihatTr(nls2:nls2+2,16:18);

% Rows 19-21
ibar(19:21,1:18) = ibar(1:18,19:21)';
ibar(19:21,19:21) = J(19:21,nls3:nls3+2)*ihatTr(nls3:nls3+2,19:21);

% Rows 22-24
ibar(22:24,1:21) = ibar(1:21,22:24)';
ibar(22:24,22:24) = J(22:24,nls4:nls4+2)*ihatTr(nls4:nls4+2,22:24);

%% Optimised the code below considering the sparsity of the Jacobian matrix and its derivative.

%psibarAux1 = -J*phihat;
%psibarAux2 = -dJ*dq;
%psibarAux3 = -ddeta;

psibarAux1(1:3) = - J(1:3,nls1:nus1)*phihat(nls1:nus1)-J(1:3,nls5:nus5)*phihat(nls5:nus5);
psibarAux1(4:6) = - J(4:6,nls2:nus2)*phihat(nls2:nus2)-J(4:6,nls5:nus5)*phihat(nls5:nus5);
psibarAux1(7:9) = - J(7:9,nls3:nus3)*phihat(nls3:nus3)-J(7:9,nls5:nus5)*phihat(nls5:nus5);
psibarAux1(10:12) =  - J(10:12,nls4:nus4)*phihat(nls4:nus4)-J(10:12,nls5:nus5)*phihat(nls5:nus5);
psibarAux1(13:15) = -J(13:15,nls1:nls1+2)*phihat(nls1:nls1+2);
psibarAux1(16:18) = -J(16:18,nls2:nls2+2)*phihat(nls2:nls2+2);
psibarAux1(19:21) = -J(19:21,nls3:nls3+2)*phihat(nls3:nls3+2); 
psibarAux1(22:24) = -J(22:24,nls4:nls4+2)*phihat(nls4:nls4+2);

psibarAux2(1:3) = -dJ(1:3,nls1:nus1)*dq(nls1:nus1)-dJ(1:3,nls5+3:nus5)*dq(nls5+3:nus5);
psibarAux2(4:6) = -dJ(4:6,nls2:nus2)*dq(nls2:nus2)-dJ(4:6,nls5+3:nus5)*dq(nls5+3:nus5);
psibarAux2(7:9) = -dJ(7:9,nls3:nus3)*dq(nls3:nus3)-dJ(7:9,nls5+3:nus5)*dq(nls5+3:nus5);
psibarAux2(10:12) =-dJ(10:12,nls4:nus4)*dq(nls4:nus4)-dJ(10:12,nls5+3:nus5)*dq(nls5+3:nus5);

psibarAux3(13:24) = -bitrj;

psibar = psibarAux1 + psibarAux2 + psibarAux3;

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

end
