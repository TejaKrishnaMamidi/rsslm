% RSSLM-CDPR-Type-I acceleration module. This module computes the Lagrange multipliers and the joint accelerations of a CDPR

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to ddq_tree_eff

% System: 6-3 CDPR with (sinusoidal or linear) cable feed

function [qdd] = acceleration(q, dq, tu, tm, tau_d)

% Global variables -- required
global alp a b th r dx dy dz m Icxx Icyy Iczz Icxy Icyz Iczx J dJ ni;
global nls1 nls2 nls3 nls4 nls5 nls6 nls7 nus1 nus2 nus3 nus4 nus5 nus6 nus7 nls7p2 nls7p3 n nc invmp;

% Global variables -- defined
global lam;

% Initialisation
phihat = zeros(n,1);
ihatTr = zeros(n,nc);
qdd = zeros(n,1);
ibar = zeros(nc,nc);
psibar = zeros(nc,1);
lam = zeros(nc,1);

phil = tau_d - tu -tm;

% Generalised inertia matrix of seventh sub-system (Based on whether I is available or not).
I7 = GIM_tree_sub_system(nls7,nus7,q,alp,a,b,th,r,dx,dy,dz,m,Icxx,Icyy,Iczz,Icxy,Icyz,Iczx);
I72 = inv3by3(I7(4:6,4:6));
%I72 = inv3by3(I(nls7p3:nus7,nls7p3:nus7));

% Instead of directly inverting the sub-system inertia matrix, used reverse Gaussian elimination technique. 

% phihat(nls1:nus1) = I(nls1:nus1,nls1:nus1)\phil(nls1:nus1);
% phihat(nls2:nus2) = I(nls2:nus2,nls2:nus2)\phil(nls2:nus2);
% phihat(nls3:nus3) = I(nls3:nus3,nls3:nus3)\phil(nls3:nus3);
% phihat(nls4:nus4) = I(nls4:nus4,nls4:nus4)\phil(nls4:nus4);
% phihat(nls5:nus5) = I(nls5:nus5,nls5:nus5)\phil(nls5:nus5);
% phihat(nls6:nus6) = I(nls6:nus6,nls6:nus6)\phil(nls6:nus6);
% phihat(nls7:nus7) = I(nls7:nus7,nls7:nus7)\phil(nls7:nus7);
[ps1,ths1,bs1,hIs1,hFs1,hGs1,hetats1,hetabs1]=ddq_tree_eff_common(q(nls1:nus1), ni(1), b(nls1:nus1), th(nls1:nus1), r(nls1:nus1), dx(nls1:nus1), dy(nls1:nus1), dz(nls1:nus1), m(nls1:nus1), Icxx(nls1:nus1), Icyy(nls1:nus1), Iczz(nls1:nus1), Icxy(nls1:nus1), Icyz(nls1:nus1), Iczx(nls1:nus1));
[ps2,ths2,bs2,hIs2,hFs2,hGs2,hetats2,hetabs2]=ddq_tree_eff_common(q(nls2:nus2), ni(2), b(nls2:nus2), th(nls2:nus2), r(nls2:nus2), dx(nls2:nus2), dy(nls2:nus2), dz(nls2:nus2), m(nls2:nus2), Icxx(nls2:nus2), Icyy(nls2:nus2), Iczz(nls2:nus2), Icxy(nls2:nus2), Icyz(nls2:nus2), Iczx(nls2:nus2));
[ps3,ths3,bs3,hIs3,hFs3,hGs3,hetats3,hetabs3]=ddq_tree_eff_common(q(nls3:nus3), ni(3), b(nls3:nus3), th(nls3:nus3), r(nls3:nus3), dx(nls3:nus3), dy(nls3:nus3), dz(nls3:nus3), m(nls3:nus3), Icxx(nls3:nus3), Icyy(nls3:nus3), Iczz(nls3:nus3), Icxy(nls3:nus3), Icyz(nls3:nus3), Iczx(nls3:nus3));
[ps4,ths4,bs4,hIs4,hFs4,hGs4,hetats4,hetabs4]=ddq_tree_eff_common(q(nls4:nus4), ni(4), b(nls4:nus4), th(nls4:nus4), r(nls4:nus4), dx(nls4:nus4), dy(nls4:nus4), dz(nls4:nus4), m(nls4:nus4), Icxx(nls4:nus4), Icyy(nls4:nus4), Iczz(nls4:nus4), Icxy(nls4:nus4), Icyz(nls4:nus4), Iczx(nls4:nus4));
[ps5,ths5,bs5,hIs5,hFs5,hGs5,hetats5,hetabs5]=ddq_tree_eff_common(q(nls5:nus5), ni(5), b(nls5:nus5), th(nls5:nus5), r(nls5:nus5), dx(nls5:nus5), dy(nls5:nus5), dz(nls5:nus5), m(nls5:nus5), Icxx(nls5:nus5), Icyy(nls5:nus5), Iczz(nls5:nus5), Icxy(nls5:nus5), Icyz(nls5:nus5), Iczx(nls5:nus5));
[ps6,ths6,bs6,hIs6,hFs6,hGs6,hetats6,hetabs6]=ddq_tree_eff_common(q(nls6:nus6), ni(6), b(nls6:nus6), th(nls6:nus6), r(nls6:nus6), dx(nls6:nus6), dy(nls6:nus6), dz(nls6:nus6), m(nls6:nus6), Icxx(nls6:nus6), Icyy(nls6:nus6), Iczz(nls6:nus6), Icxy(nls6:nus6), Icyz(nls6:nus6), Iczx(nls6:nus6));
%
phihat(nls1:nus1) = ddq_tree_eff_variable(ni(1), alp(nls1:nus1), a(nls1:nus1),[0;cumsum(ones(ni(1)-1,1))], r(nls1:nus1), ps1, ths1, bs1, hIs1, hFs1, hGs1, hetats1, hetabs1, phil(nls1:nus1));
phihat(nls2:nus2) = ddq_tree_eff_variable(ni(2), alp(nls2:nus2), a(nls2:nus2),[0;cumsum(ones(ni(2)-1,1))], r(nls2:nus2), ps2, ths2, bs2, hIs2, hFs2, hGs2, hetats2, hetabs2, phil(nls2:nus2));
phihat(nls3:nus3) = ddq_tree_eff_variable(ni(3), alp(nls3:nus3), a(nls3:nus3),[0;cumsum(ones(ni(3)-1,1))], r(nls3:nus3), ps3, ths3, bs3, hIs3, hFs3, hGs3, hetats3, hetabs3, phil(nls3:nus3));
phihat(nls4:nus4) = ddq_tree_eff_variable(ni(4), alp(nls4:nus4), a(nls4:nus4),[0;cumsum(ones(ni(4)-1,1))], r(nls4:nus4), ps4, ths4, bs4, hIs4, hFs4, hGs4, hetats4, hetabs4, phil(nls4:nus4));
phihat(nls5:nus5) = ddq_tree_eff_variable(ni(5), alp(nls5:nus5), a(nls5:nus5),[0;cumsum(ones(ni(5)-1,1))], r(nls5:nus5), ps5, ths5, bs5, hIs5, hFs5, hGs5, hetats5, hetabs5, phil(nls5:nus5));
phihat(nls6:nus6) = ddq_tree_eff_variable(ni(6), alp(nls6:nus6), a(nls6:nus6),[0;cumsum(ones(ni(6)-1,1))], r(nls6:nus6), ps6, ths6, bs6, hIs6, hFs6, hGs6, hetats6, hetabs6, phil(nls6:nus6));
%phihat(nls7:nus7) = ddq_tree_eff(q(nls7:nus7), ni(7), alp(nls7:nus7), a(nls7:nus7), b(nls7:nus7), th(nls7:nus7), [0;cumsum(ones(ni(7)-1,1))], r(nls7:nus7), dx(nls7:nus7), dy(nls7:nus7), dz(nls7:nus7), m(nls7:nus7), Icxx(nls7:nus7), Icyy(nls7:nus7), Iczz(nls7:nus7), Icxy(nls7:nus7), Icyz(nls7:nus7), Iczx(nls7:nus7), phil(nls7:nus7));
phihat(nls7:nls7p2) = invmp.*phil(nls7:nls7p2);
phihat(nls7p3:nus7) = I72*phil(nls7p3:nus7);

%disp([cond(I),cond(I(nls1:nus1,nls1:nus1)),cond(I(nls2:nus2,nls2:nus2)),cond(I(nls3:nus3,nls3:nus3)),cond(I(nls4:nus4,nls4:nus4)),cond(I(nls5:nus5,nls5:nus5)),cond(I(nls6:nus6,nls6:nus6))]);

%% Optimised the code below considering the sparsity of the Jacobian matrix (6-3 CDPR).

% ihatTr(nls1:nus1,1:18) = I(nls1:nus1,nls1:nus1)\J(:,nls1:nus1)';
% ihatTr(nls2:nus2,1:18) = I(nls2:nus2,nls2:nus2)\J(:,nls2:nus2)';
% ihatTr(nls3:nus3,1:18) = I(nls3:nus3,nls3:nus3)\J(:,nls3:nus3)';
% ihatTr(nls4:nus4,1:18) = I(nls4:nus4,nls4:nus4)\J(:,nls4:nus4)';
% ihatTr(nls5:nus5,1:18) = I(nls5:nus5,nls5:nus5)\J(:,nls5:nus5)';
% ihatTr(nls6:nus6,1:18) = I(nls6:nus6,nls6:nus6)\J(:,nls6:nus6)';
% ihatTr(nls7:nus7,1:18) = I(nls7:nus7,nls7:nus7)\J(:,nls7:nus7)';

for ii = 1:1:3
    ihatTr(nls1:nus1,ii) = ddq_tree_eff_variable(ni(1), alp(nls1:nus1), a(nls1:nus1),[0;cumsum(ones(ni(1)-1,1))], r(nls1:nus1), ps1, ths1, bs1, hIs1, hFs1, hGs1, hetats1, hetabs1, J(ii,nls1:nus1)');
    %ihatTr(nls7:nus7,ii) = ddq_tree_eff(q(nls7:nus7), ni(7), alp(nls7:nus7), a(nls7:nus7), b(nls7:nus7), th(nls7:nus7), [0;cumsum(ones(ni(7)-1,1))], r(nls7:nus7), dx(nls7:nus7), dy(nls7:nus7), dz(nls7:nus7), m(nls7:nus7), Icxx(nls7:nus7), Icyy(nls7:nus7), Iczz(nls7:nus7), Icxy(nls7:nus7), Icyz(nls7:nus7), Iczx(nls7:nus7), J(ii,nls7:nus7)');
    ihatTr(nls7:nls7p2,ii) = invmp.*J(ii,nls7:nls7p2)';
    ihatTr(nls7p3:nus7,ii) = I72*J(ii,nls7p3:nus7)';
end

for ii = 4:1:6
    ihatTr(nls2:nus2,ii) = ddq_tree_eff_variable(ni(2), alp(nls2:nus2), a(nls2:nus2),[0;cumsum(ones(ni(2)-1,1))], r(nls2:nus2), ps2, ths2, bs2, hIs2, hFs2, hGs2, hetats2, hetabs2, J(ii,nls2:nus2)');
    %ihatTr(nls7:nus7,ii) = ddq_tree_eff(q(nls7:nus7), ni(7), alp(nls7:nus7), a(nls7:nus7), b(nls7:nus7), th(nls7:nus7), [0;cumsum(ones(ni(7)-1,1))], r(nls7:nus7), dx(nls7:nus7), dy(nls7:nus7), dz(nls7:nus7), m(nls7:nus7), Icxx(nls7:nus7), Icyy(nls7:nus7), Iczz(nls7:nus7), Icxy(nls7:nus7), Icyz(nls7:nus7), Iczx(nls7:nus7), J(ii,nls7:nus7)');
    ihatTr(nls7:nls7p2,ii) = invmp.*J(ii,nls7:nls7p2)';
    ihatTr(nls7p3:nus7,ii) = I72*J(ii,nls7p3:nus7)';
end

for ii = 7:1:9
    ihatTr(nls3:nus3,ii) = ddq_tree_eff_variable(ni(3), alp(nls3:nus3), a(nls3:nus3),[0;cumsum(ones(ni(3)-1,1))], r(nls3:nus3), ps3, ths3, bs3, hIs3, hFs3, hGs3, hetats3, hetabs3, J(ii,nls3:nus3)');
    %ihatTr(nls7:nus7,ii) = ddq_tree_eff(q(nls7:nus7), ni(7), alp(nls7:nus7), a(nls7:nus7), b(nls7:nus7), th(nls7:nus7), [0;cumsum(ones(ni(7)-1,1))], r(nls7:nus7), dx(nls7:nus7), dy(nls7:nus7), dz(nls7:nus7), m(nls7:nus7), Icxx(nls7:nus7), Icyy(nls7:nus7), Iczz(nls7:nus7), Icxy(nls7:nus7), Icyz(nls7:nus7), Iczx(nls7:nus7), J(ii,nls7:nus7)');
    ihatTr(nls7:nls7p2,ii) = invmp.*J(ii,nls7:nls7p2)';
    ihatTr(nls7p3:nus7,ii) = I72*J(ii,nls7p3:nus7)';
end

for ii = 10:1:12
    ihatTr(nls4:nus4,ii) = ddq_tree_eff_variable(ni(4), alp(nls4:nus4), a(nls4:nus4),[0;cumsum(ones(ni(4)-1,1))], r(nls4:nus4), ps4, ths4, bs4, hIs4, hFs4, hGs4, hetats4, hetabs4, J(ii,nls4:nus4)');
    %ihatTr(nls7:nus7,ii) = ddq_tree_eff(q(nls7:nus7), ni(7), alp(nls7:nus7), a(nls7:nus7), b(nls7:nus7), th(nls7:nus7), [0;cumsum(ones(ni(7)-1,1))], r(nls7:nus7), dx(nls7:nus7), dy(nls7:nus7), dz(nls7:nus7), m(nls7:nus7), Icxx(nls7:nus7), Icyy(nls7:nus7), Iczz(nls7:nus7), Icxy(nls7:nus7), Icyz(nls7:nus7), Iczx(nls7:nus7), J(ii,nls7:nus7)');
    ihatTr(nls7:nls7p2,ii) = invmp.*J(ii,nls7:nls7p2)';
    ihatTr(nls7p3:nus7,ii) = I72*J(ii,nls7p3:nus7)';
end

for ii = 13:1:15
    ihatTr(nls5:nus5,ii) = ddq_tree_eff_variable(ni(5), alp(nls5:nus5), a(nls5:nus5),[0;cumsum(ones(ni(5)-1,1))], r(nls5:nus5), ps5, ths5, bs5, hIs5, hFs5, hGs5, hetats5, hetabs5, J(ii,nls5:nus5)');
    %ihatTr(nls7:nus7,ii) = ddq_tree_eff(q(nls7:nus7), ni(7), alp(nls7:nus7), a(nls7:nus7), b(nls7:nus7), th(nls7:nus7), [0;cumsum(ones(ni(7)-1,1))], r(nls7:nus7), dx(nls7:nus7), dy(nls7:nus7), dz(nls7:nus7), m(nls7:nus7), Icxx(nls7:nus7), Icyy(nls7:nus7), Iczz(nls7:nus7), Icxy(nls7:nus7), Icyz(nls7:nus7), Iczx(nls7:nus7), J(ii,nls7:nus7)');
    ihatTr(nls7:nls7p2,ii) = invmp.*J(ii,nls7:nls7p2)';
    ihatTr(nls7p3:nus7,ii) = I72*J(ii,nls7p3:nus7)';
end

for ii = 16:1:18
    ihatTr(nls6:nus6,ii) = ddq_tree_eff_variable(ni(6), alp(nls6:nus6), a(nls6:nus6),[0;cumsum(ones(ni(6)-1,1))], r(nls6:nus6), ps6, ths6, bs6, hIs6, hFs6, hGs6, hetats6, hetabs6, J(ii,nls6:nus6)');
    %ihatTr(nls7:nus7,ii) = ddq_tree_eff(q(nls7:nus7), ni(7), alp(nls7:nus7), a(nls7:nus7), b(nls7:nus7), th(nls7:nus7), [0;cumsum(ones(ni(7)-1,1))], r(nls7:nus7), dx(nls7:nus7), dy(nls7:nus7), dz(nls7:nus7), m(nls7:nus7), Icxx(nls7:nus7), Icyy(nls7:nus7), Iczz(nls7:nus7), Icxy(nls7:nus7), Icyz(nls7:nus7), Iczx(nls7:nus7), J(ii,nls7:nus7)');
    ihatTr(nls7:nls7p2,ii) = invmp.*J(ii,nls7:nls7p2)';
    ihatTr(nls7p3:nus7,ii) = I72*J(ii,nls7p3:nus7)';
end

%% Optimised the code below considering the sparsity of the Jacobian matrix and the properties of the matrix ibar.

%ibar = J*ihatTr;

% Rows 1-3
ibar(1:3,1:3) = J(1:3,nls1:nus1)*ihatTr(nls1:nus1,1:3)+J(1:3,nls7:nus7)*ihatTr(nls7:nus7,1:3);
ibar(1:3,4:18) = J(1:3,nls7:nus7)*ihatTr(nls7:nus7,4:18);

% Rows 4-6
ibar(4:6,1:3) = ibar(1:3,4:6)';
ibar(4:6,4:6) = J(4:6,nls2:nus2)*ihatTr(nls2:nus2,4:6)+ J(4:6,nls7:nus7)*ihatTr(nls7:nus7,4:6);
ibar(4:6,7:18) = J(4:6,nls7:nus7)*ihatTr(nls7:nus7,7:18);

% Rows 7-9
ibar(7:9,1:6) = ibar(1:6,7:9)';
ibar(7:9,7:9) = J(7:9,nls3:nus3)*ihatTr(nls3:nus3,7:9)+ J(7:9, nls7:nus7)*ihatTr(nls7:nus7,7:9);
ibar(7:9,10:18) = J(7:9,nls7:nus7)*ihatTr(nls7:nus7,10:18);

% Rows 10-12
ibar(10:12,1:9) = ibar(1:9,10:12)';
ibar(10:12,10:12) = J(10:12,nls4:nus4)*ihatTr(nls4:nus4,10:12) + J(10:12,nls7:nus7)*ihatTr(nls7:nus7, 10:12);
ibar(10:12,13:18) = J(10:12,nls7:nus7)*ihatTr(nls7:nus7,13:18);

% Rows 13-15
ibar(13:15,1:12) = ibar(1:12,13:15)';
ibar(13:15,13:15) = J(13:15,nls5:nus5)*ihatTr(nls5:nus5,13:15) + J(13:15, nls7:nus7)*ihatTr(nls7:nus7, 13:15);
ibar(13:15,16:18) = J(13:15,nls7:nus7)*ihatTr(nls7:nus7,16:18);

% Rows 16-18
ibar(16:18,1:15) = ibar(1:15,16:18)';
ibar(16:18,16:18) = J(16:18,nls6:nus6)*ihatTr(nls6:nus6,16:18) + J(16:18, nls7:nus7)*ihatTr(nls7:nus7, 16:18);

%% Optimised the code below considering the sparsity of the Jacobian matrix and its derivative.

% psibar = -dJ*dq-J*phihat;

psibar(1:3) = -dJ(1:3,nls1:nus1)*dq(nls1:nus1)-dJ(1:3,nls7+3:nus7)*dq(nls7+3:nus7) - J(1:3,nls1:nus1)*phihat(nls1:nus1)-J(1:3,nls7:nus7)*phihat(nls7:nus7);
psibar(4:6) = -dJ(4:6,nls2:nus2)*dq(nls2:nus2)-dJ(4:6,nls7+3:nus7)*dq(nls7+3:nus7) - J(4:6,nls2:nus2)*phihat(nls2:nus2)-J(4:6,nls7:nus7)*phihat(nls7:nus7);
psibar(7:9) = -dJ(7:9,nls3:nus3)*dq(nls3:nus3)-dJ(7:9,nls7+3:nus7)*dq(nls7+3:nus7) - J(7:9,nls3:nus3)*phihat(nls3:nus3)-J(7:9,nls7:nus7)*phihat(nls7:nus7);
psibar(10:12) = -dJ(10:12,nls4:nus4)*dq(nls4:nus4)-dJ(10:12,nls7+3:nus7)*dq(nls7+3:nus7) - J(10:12,nls4:nus4)*phihat(nls4:nus4)-J(10:12,nls7:nus7)*phihat(nls7:nus7);
psibar(13:15) = -dJ(13:15,nls5:nus5)*dq(nls5:nus5)-dJ(13:15,nls7+3:nus7)*dq(nls7+3:nus7) - J(13:15,nls5:nus5)*phihat(nls5:nus5)-J(13:15,nls7:nus7)*phihat(nls7:nus7);
psibar(16:18) = -dJ(16:18,nls6:nus6)*dq(nls6:nus6)-dJ(16:18,nls7+3:nus7)*dq(nls7+3:nus7) - J(16:18,nls6:nus6)*phihat(nls6:nus6)-J(16:18,nls7:nus7)*phihat(nls7:nus7);

%% Optimised the code below by performing block Gaussian elimination.
% lam = ibar\psibar;

auxa  = inv3by3(ibar(1:3,1:3));
auxa1 = ibar(4:6,1:3)*auxa;
auxa2 = ibar(7:9,1:3)*auxa;
auxa3 = ibar(10:12,1:3)*auxa;
auxa4 = ibar(13:15,1:3)*auxa;
auxa5 = ibar(16:18,1:3)*auxa;
%
auxb  = inv3by3(ibar(4:6,4:6)-auxa1*ibar(1:3,4:6));
auxb1 = (ibar(7:9,4:6)-auxa2*ibar(1:3,4:6))*auxb;
auxb2 = (ibar(10:12,4:6)-auxa3*ibar(1:3,4:6))*auxb;
auxb3 = (ibar(13:15,4:6)-auxa4*ibar(1:3,4:6))*auxb;
auxb4 = (ibar(16:18,4:6)-auxa5*ibar(1:3,4:6))*auxb;
%
auxc0 = ibar(4:6,7:9)-auxa1*ibar(1:3,7:9);
auxc  = inv3by3(ibar(7:9,7:9)-auxa2*ibar(1:3,7:9)-auxb1*auxc0);
auxc1 = (ibar(10:12,7:9)-auxa3*ibar(1:3,7:9)-auxb2*auxc0)*auxc;
auxc2 = (ibar(13:15,7:9)-auxa4*ibar(1:3,7:9)-auxb3*auxc0)*auxc;
auxc3 = (ibar(16:18,7:9)-auxa5*ibar(1:3,7:9)-auxb4*auxc0)*auxc;
%
auxd01 = ibar(4:6,10:12)-auxa1*ibar(1:3,10:12);
auxd02 = ibar(7:9,10:12)-auxa2*ibar(1:3,10:12)-auxb1*auxd01; 
auxd = inv3by3(ibar(10:12,10:12)-auxa3*ibar(1:3,10:12)-auxb2*auxd01-auxc1*auxd02);
auxd1 = (ibar(13:15,10:12)-auxa4*ibar(1:3,10:12)-auxb3*auxd01-auxc2*auxd02)*auxd;
auxd2 = (ibar(16:18,10:12)-auxa5*ibar(1:3,10:12)-auxb4*auxd01-auxc3*auxd02)*auxd;
%
auxe01 = ibar(4:6,13:15)-auxa1*ibar(1:3,13:15);
auxe02 = ibar(7:9,13:15)-auxa2*ibar(1:3,13:15)-auxb1*auxe01;
auxe03 = ibar(10:12,13:15)-auxa3*ibar(1:3,13:15)-auxb2*auxe01-auxc1*auxe02;
auxe = inv3by3(ibar(13:15,13:15)-auxa4*ibar(1:3,13:15)-auxb3*auxe01-auxc2*auxe02-auxd1*auxe03);
auxe1 = (ibar(16:18,13:15)-auxa5*ibar(1:3,13:15)-auxb4*auxe01-auxc3*auxe02-auxd2*auxe03)*auxe;
%
auxf01 = ibar(4:6,16:18)-auxa1*ibar(1:3,16:18);
auxf02 = ibar(7:9,16:18)-auxa2*ibar(1:3,16:18)-auxb1*auxf01;
auxf03 = ibar(10:12,16:18)-auxa3*ibar(1:3,16:18)-auxb2*auxf01-auxc1*auxf02;
auxf04 = ibar(13:15,16:18)-auxa4*ibar(1:3,16:18)-auxb3*auxf01-auxc2*auxf02-auxd1*auxf03;
%
auxp1 = psibar(4:6)-auxa1*psibar(1:3);
auxp2 = psibar(7:9)-auxa2*psibar(1:3)-auxb1*auxp1;
auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;

% Back substitution

lam(16:18) = inv3by3(ibar(16:18,16:18)-auxa5*ibar(1:3,16:18)-auxb4*auxf01-auxc3*auxf02-auxd2*auxf03-auxe1*auxf04)*auxp5;
lam(13:15) = auxe*(auxp4-auxf04*lam(16:18));
lam(10:12) = auxd*(auxp3-auxe03*lam(13:15)-auxf03*lam(16:18));
lam(7:9)   = auxc*(auxp2-auxd02*lam(10:12)-auxe02*lam(13:15)-auxf02*lam(16:18));
lam(4:6)   = auxb*(auxp1-auxc0*lam(7:9)-auxd01*lam(10:12)-auxe01*lam(13:15)-auxf01*lam(16:18));
lam(1:3)   = auxa*(psibar(1:3)-ibar(1:3,4:6)*lam(4:6)-ibar(1:3,7:9)*lam(7:9)-ibar(1:3,10:12)*lam(10:12)-ibar(1:3,13:15)*lam(13:15)-ibar(1:3,16:18)*lam(16:18));

%% Joint accelerations

% Optimised the code below using the sparsity of ihatTr
qdd(nls1:nus1) = phihat(nls1:nus1)+ihatTr(nls1:nus1,1:3)*lam(1:3);
qdd(nls2:nus2) = phihat(nls2:nus2)+ihatTr(nls2:nus2,4:6)*lam(4:6);
qdd(nls3:nus3) = phihat(nls3:nus3)+ihatTr(nls3:nus3,7:9)*lam(7:9);
qdd(nls4:nus4) = phihat(nls4:nus4)+ihatTr(nls4:nus4,10:12)*lam(10:12);
qdd(nls5:nus5) = phihat(nls5:nus5)+ihatTr(nls5:nus5,13:15)*lam(13:15);
qdd(nls6:nus6) = phihat(nls6:nus6)+ihatTr(nls6:nus6,16:18)*lam(16:18);
qdd(nls7:nus7) = phihat(nls7:nus7)+ihatTr(nls7:nus7,:)*lam;

%qdd(nls1:nus1) = phihat(nls1:nus1)+ihatTr(nls1:nus1,1:18)*lam;
%qdd(nls2:nus2) = phihat(nls2:nus2)+ihatTr(nls2:nus2,1:18)*lam;
%qdd(nls3:nus3) = phihat(nls3:nus3)+ihatTr(nls3:nus3,1:18)*lam;
%qdd(nls4:nus4) = phihat(nls4:nus4)+ihatTr(nls4:nus4,1:18)*lam;
%qdd(nls5:nus5) = phihat(nls5:nus5)+ihatTr(nls5:nus5,1:18)*lam;
%qdd(nls6:nus6) = phihat(nls6:nus6)+ihatTr(nls6:nus6,1:18)*lam;
%qdd(nls7:nus7) = phihat(nls7:nus7)+ihatTr(nls7:nus7,1:18)*lam;

end
