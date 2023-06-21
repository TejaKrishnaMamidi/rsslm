% RSSLM-CDPR-Kinetostatics acceleration_var module. The varying portions in determining the changes in joint accelerations with perturbations in the joint positions and velocities are computed here. 

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to ddq_tree_eff

% System: 8-8 CDPR

function [qdd] = acceleration_var(jj, fl, q, dq, tu, tau_d, phihat, ihatTr, ibar, psibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxa6, auxa7, auxb, auxb1, auxb2, auxb3, auxb4, auxb5, auxb6, auxc0, auxc, auxc1, auxc2, auxc3, auxc4, auxc5, auxd01, auxd02, auxd, auxd1, auxd2, auxd3, auxd4, auxe01, auxe02, auxe03, auxe, auxe1, auxe2, auxe3, auxf01, auxf02, auxf03, auxf04, auxf, auxf1, auxf2, auxg01, auxg02, auxg03, auxg04, auxg05, auxg, auxg1, auxh01, auxh02, auxh03, auxh04, auxh05, auxh06, auxh, auxp1, auxp2, auxp3, auxp4, auxp5, auxp6, auxp7)

% Global variables -- required
global alp a b th r dx dy dz m Icxx Icyy Iczz Icxy Icyz Iczx J dJ ni;
% Limits for the sub-systems
global nls1 nls2 nls3 nls4 nls5 nls6 nls7 nls8 nls9 nus1 nus2 nus3 nus4 nus5 nus6 nus7 nus8 nus9 nls9p2 nls9p3 invmp;

% Initialisation
qdd = zeros(nus9,1);
lam = zeros(24,1);

phil = tau_d - tu;

% Generalised inertia matrix of ninth sub-system
I9 = GIM_tree_sub_system(nls9,nus9,q,alp,a,b,th,r,dx,dy,dz,m,Icxx,Icyy,Iczz,Icxy,Icyz,Iczx);
I92 = inv3by3(I9(4:6,4:6));

% Instead of directly inverting the sub-system inertia matrix, used reverse Gaussian elimination technique powered by DeNoC. 
%phihat(lsi:lui) = I(lsi:lui,lsi:lui)\phil(lsi:lui);

% Optimised the code below considering the sparsity of the Jacobian matrix (8-8 CDPR).
% ihatTr(lsi:lui,1:15) = I(lsi:lui,lsi:lui)\J(:,lsi:lui)'

% Optimised the code below considering the sparsity of the Jacobian matrix and the properties of the matrix ibar.
% ibar = J*ihatTr;

%psibarAux1 = -J*phihat;
%psibarAux2 = -dJ*dq;

if fl == 0
    % sub-system 1
    if jj<=nus1
        [ps1,ths1,bs1,hIs1,hFs1,hGs1,hetats1,hetabs1]=ddq_tree_eff_common(q(nls1:nus1), ni(1), b(nls1:nus1), th(nls1:nus1), r(nls1:nus1), dx(nls1:nus1), dy(nls1:nus1), dz(nls1:nus1), m(nls1:nus1), Icxx(nls1:nus1), Icyy(nls1:nus1), Iczz(nls1:nus1), Icxy(nls1:nus1), Icyz(nls1:nus1), Iczx(nls1:nus1));
        phihat(nls1:nus1) = ddq_tree_eff_variable(ni(1), alp(nls1:nus1), a(nls1:nus1),[0;cumsum(ones(ni(1)-1,1))], r(nls1:nus1), ps1, ths1, bs1, hIs1, hFs1, hGs1, hetats1, hetabs1, phil(nls1:nus1));
        for ii = 1:1:3
            ihatTr(nls1:nus1,ii) = ddq_tree_eff_variable(ni(1), alp(nls1:nus1), a(nls1:nus1),[0;cumsum(ones(ni(1)-1,1))], r(nls1:nus1), ps1, ths1, bs1, hIs1, hFs1, hGs1, hetats1, hetabs1, J(ii,nls1:nus1)');
        end
        % Ibar
    	ibar(1:3,1:3) = J(1:3,nls1:nus1)*ihatTr(nls1:nus1,1:3)+J(1:3,nls9:nus9)*ihatTr(nls9:nus9,1:3);
        % Psibar
        psibarAux1(1:3) = - J(1:3,nls1:nus1)*phihat(nls1:nus1)-J(1:3,nls9:nus9)*phihat(nls9:nus9);
        psibar(1:3) = psibarAux1(1:3)+psibarAux2(1:3);
        % Computation of lambda
        auxa  = inv3by3(ibar(1:3,1:3));
        auxa1 = ibar(4:6,1:3)*auxa;
    	auxa2 = ibar(7:9,1:3)*auxa;
    	auxa3 = ibar(10:12,1:3)*auxa;
        auxa4 = ibar(13:15,1:3)*auxa;
        auxa5 = ibar(16:18,1:3)*auxa;
        auxa6 = ibar(19:21,1:3)*auxa;
        auxa7 = ibar(22:24,1:3)*auxa;
        %
        auxb  = inv3by3(ibar(4:6,4:6)-auxa1*ibar(1:3,4:6));
        auxb1 = (ibar(7:9,4:6)-auxa2*ibar(1:3,4:6))*auxb;
        auxb2 = (ibar(10:12,4:6)-auxa3*ibar(1:3,4:6))*auxb;
        auxb3 = (ibar(13:15,4:6)-auxa4*ibar(1:3,4:6))*auxb;
        auxb4 = (ibar(16:18,4:6)-auxa5*ibar(1:3,4:6))*auxb;
        auxb5 = (ibar(19:21,4:6)-auxa6*ibar(1:3,4:6))*auxb;
        auxb6 = (ibar(22:24,4:6)-auxa7*ibar(1:3,4:6))*auxb;
        %
        auxc0 = ibar(4:6,7:9)-auxa1*ibar(1:3,7:9);
        auxc  = inv3by3(ibar(7:9,7:9)-auxa2*ibar(1:3,7:9)-auxb1*auxc0);
        auxc1 = (ibar(10:12,7:9)-auxa3*ibar(1:3,7:9)-auxb2*auxc0)*auxc;
        auxc2 = (ibar(13:15,7:9)-auxa4*ibar(1:3,7:9)-auxb3*auxc0)*auxc;
        auxc3 = (ibar(16:18,7:9)-auxa5*ibar(1:3,7:9)-auxb4*auxc0)*auxc;
        auxc4 = (ibar(19:21,7:9)-auxa6*ibar(1:3,7:9)-auxb5*auxc0)*auxc;
        auxc5 = (ibar(22:24,7:9)-auxa7*ibar(1:3,7:9)-auxb6*auxc0)*auxc;
        %
        auxd01 = ibar(4:6,10:12)-auxa1*ibar(1:3,10:12);
        auxd02 = ibar(7:9,10:12)-auxa2*ibar(1:3,10:12)-auxb1*auxd01; 
        auxd = inv3by3(ibar(10:12,10:12)-auxa3*ibar(1:3,10:12)-auxb2*auxd01-auxc1*auxd02);
    	auxd1 = (ibar(13:15,10:12)-auxa4*ibar(1:3,10:12)-auxb3*auxd01-auxc2*auxd02)*auxd;
    	auxd2 = (ibar(16:18,10:12)-auxa5*ibar(1:3,10:12)-auxb4*auxd01-auxc3*auxd02)*auxd;
    	auxd3 = (ibar(19:21,10:12)-auxa6*ibar(1:3,10:12)-auxb5*auxd01-auxc4*auxd02)*auxd;
    	auxd4 = (ibar(22:24,10:12)-auxa7*ibar(1:3,10:12)-auxb6*auxd01-auxc5*auxd02)*auxd;
        %
        auxe01 = ibar(4:6,13:15)-auxa1*ibar(1:3,13:15);
        auxe02 = ibar(7:9,13:15)-auxa2*ibar(1:3,13:15)-auxb1*auxe01;
    	auxe03 = ibar(10:12,13:15)-auxa3*ibar(1:3,13:15)-auxb2*auxe01-auxc1*auxe02;
        auxe = inv3by3(ibar(13:15,13:15)-auxa4*ibar(1:3,13:15)-auxb3*auxe01-auxc2*auxe02-auxd1*auxe03);
        auxe1 = (ibar(16:18,13:15)-auxa5*ibar(1:3,13:15)-auxb4*auxe01-auxc3*auxe02-	auxd2*auxe03)*auxe;
        auxe2 = (ibar(19:21,13:15)-auxa6*ibar(1:3,13:15)-auxb5*auxe01-auxc4*auxe02-auxd3*auxe03)*auxe;
        auxe3 = (ibar(22:24,13:15)-auxa7*ibar(1:3,13:15)-auxb6*auxe01-auxc5*auxe02-auxd4*auxe03)*auxe;
        %
        auxf01 = ibar(4:6,16:18)-auxa1*ibar(1:3,16:18);
        auxf02 = ibar(7:9,16:18)-auxa2*ibar(1:3,16:18)-auxb1*auxf01;
        auxf03 = ibar(10:12,16:18)-auxa3*ibar(1:3,16:18)-auxb2*auxf01-auxc1*auxf02;
        auxf04 = ibar(13:15,16:18)-auxa4*ibar(1:3,16:18)-auxb3*auxf01-auxc2*auxf02-auxd1*auxf03;
        auxf = inv3by3(ibar(16:18,16:18)-auxa5*ibar(1:3,16:18)-auxb4*auxf01-auxc3*auxf02-auxd2*auxf03-auxe1*auxf04);
        auxf1 = (ibar(19:21,16:18)-auxa6*ibar(1:3,16:18)-auxb5*auxf01-auxc4*auxf02-auxd3*auxf03-auxe2*auxf04)*auxf;
        auxf2 = (ibar(22:24,16:18)-auxa7*ibar(1:3,16:18)-auxb6*auxf01-auxc5*auxf02-auxd4*auxf03-auxe3*auxf04)*auxf;
        %
        auxg01 = ibar(4:6,19:21)-auxa1*ibar(1:3,19:21);
        auxg02 = ibar(7:9,19:21)-auxa2*ibar(1:3,19:21)-auxb1*auxg01;
        auxg03 = ibar(10:12,19:21)-auxa3*ibar(1:3,19:21)-auxb2*auxg01-auxc1*auxg02;
        auxg04 = ibar(13:15,19:21)-auxa4*ibar(1:3,19:21)-auxb3*auxg01-auxc2*auxg02-auxd1*auxg03;
        auxg05 = ibar(16:18,19:21)-auxa5*ibar(1:3,19:21)-auxb4*auxg01-auxc3*auxg02-auxd2*auxg03-auxe1*auxg04;
        auxg = inv3by3(ibar(19:21,19:21)-auxa6*ibar(1:3,19:21)-auxb5*auxg01-auxc4*auxg02-auxd3*auxg03-auxe2*auxg04-auxf1*auxg05);
        auxg1 = (ibar(22:24,19:21)-auxa7*ibar(1:3,19:21)-auxb6*auxg01-auxc5*auxg02-auxd4*auxg03-auxe3*auxg04-auxf2*auxg05)*auxg;
        %
        auxh01 = ibar(4:6,22:24)-auxa1*ibar(1:3,22:24);
        auxh02 = ibar(7:9,22:24)-auxa2*ibar(1:3,22:24)-auxb1*auxh01;
        auxh03 = ibar(10:12,22:24)-auxa3*ibar(1:3,22:24)-auxb2*auxh01-auxc1*auxh02;
        auxh04 = ibar(13:15,22:24)-auxa4*ibar(1:3,22:24)-auxb3*auxh01-auxc2*auxh02-auxd1*auxh03;
        auxh05 = ibar(16:18,22:24)-auxa5*ibar(1:3,22:24)-auxb4*auxh01-auxc3*auxh02-auxd2*auxh03-auxe1*auxh04;
        auxh06 = ibar(19:21,22:24)-auxa6*ibar(1:3,22:24)-auxb5*auxh01-auxc4*auxh02-auxd3*auxh03-auxe2*auxh04-auxf1*auxh05;
        auxh = inv3by3(ibar(22:24,22:24)-auxa7*ibar(1:3,22:24)-auxb6*auxh01-auxc5*auxh02-auxd4*auxh03-auxe3*auxh04-auxf2*auxh05-auxg1*auxh06);
        %
        auxp1 = psibar(4:6)-auxa1*psibar(1:3);
        auxp2 = psibar(7:9)-auxa2*psibar(1:3)-auxb1*auxp1;
        auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
        auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 2
    elseif jj<=nus2
        [ps2,ths2,bs2,hIs2,hFs2,hGs2,hetats2,hetabs2]=ddq_tree_eff_common(q(nls2:nus2), ni(2), b(nls2:nus2), th(nls2:nus2), r(nls2:nus2), dx(nls2:nus2), dy(nls2:nus2), dz(nls2:nus2), m(nls2:nus2), Icxx(nls2:nus2), Icyy(nls2:nus2), Iczz(nls2:nus2), Icxy(nls2:nus2), Icyz(nls2:nus2), Iczx(nls2:nus2));
        phihat(nls2:nus2) = ddq_tree_eff_variable(ni(2), alp(nls2:nus2), a(nls2:nus2),[0;cumsum(ones(ni(2)-1,1))], r(nls2:nus2), ps2, ths2, bs2, hIs2, hFs2, hGs2, hetats2, hetabs2, phil(nls2:nus2));
        for ii = 4:1:6
            ihatTr(nls2:nus2,ii) = ddq_tree_eff_variable(ni(2), alp(nls2:nus2), a(nls2:nus2),[0;cumsum(ones(ni(2)-1,1))], r(nls2:nus2), ps2, ths2, bs2, hIs2, hFs2, hGs2, hetats2, hetabs2, J(ii,nls2:nus2)');
        end
    	% Ibar
    	ibar(4:6,4:6) = J(4:6,nls2:nus2)*ihatTr(nls2:nus2,4:6)+ J(4:6,nls9:nus9)*ihatTr(nls9:nus9,4:6);
        % Psibar
        psibarAux1(4:6) = - J(4:6,nls2:nus2)*phihat(nls2:nus2)-J(4:6,nls9:nus9)*phihat(nls9:nus9);
    	psibar(4:6) = psibarAux1(4:6)+psibarAux2(4:6);
        % Computation of lambda
        auxb  = inv3by3(ibar(4:6,4:6)-auxa1*ibar(1:3,4:6));
        auxb1 = (ibar(7:9,4:6)-auxa2*ibar(1:3,4:6))*auxb;
    	auxb2 = (ibar(10:12,4:6)-auxa3*ibar(1:3,4:6))*auxb;
    	auxb3 = (ibar(13:15,4:6)-auxa4*ibar(1:3,4:6))*auxb;
    	auxb4 = (ibar(16:18,4:6)-auxa5*ibar(1:3,4:6))*auxb;
    	auxb5 = (ibar(19:21,4:6)-auxa6*ibar(1:3,4:6))*auxb;
        auxb6 = (ibar(22:24,4:6)-auxa7*ibar(1:3,4:6))*auxb;
        %
        auxc  = inv3by3(ibar(7:9,7:9)-auxa2*ibar(1:3,7:9)-auxb1*auxc0);
    	auxc1 = (ibar(10:12,7:9)-auxa3*ibar(1:3,7:9)-auxb2*auxc0)*auxc;	
        auxc2 = (ibar(13:15,7:9)-auxa4*ibar(1:3,7:9)-auxb3*auxc0)*auxc;
        auxc3 = (ibar(16:18,7:9)-auxa5*ibar(1:3,7:9)-auxb4*auxc0)*auxc;
        auxc4 = (ibar(19:21,7:9)-auxa6*ibar(1:3,7:9)-auxb5*auxc0)*auxc;
        auxc5 = (ibar(22:24,7:9)-auxa7*ibar(1:3,7:9)-auxb6*auxc0)*auxc;
        %
        auxd02 = ibar(7:9,10:12)-auxa2*ibar(1:3,10:12)-auxb1*auxd01; 
        auxd = inv3by3(ibar(10:12,10:12)-auxa3*ibar(1:3,10:12)-auxb2*auxd01-auxc1*auxd02);
        auxd1 = (ibar(13:15,10:12)-auxa4*ibar(1:3,10:12)-auxb3*auxd01-auxc2*auxd02)*auxd;
        auxd2 = (ibar(16:18,10:12)-auxa5*ibar(1:3,10:12)-auxb4*auxd01-auxc3*auxd02)*auxd;
        auxd3 = (ibar(19:21,10:12)-auxa6*ibar(1:3,10:12)-auxb5*auxd01-auxc4*auxd02)*auxd;
        auxd4 = (ibar(22:24,10:12)-auxa7*ibar(1:3,10:12)-auxb6*auxd01-auxc5*auxd02)*auxd;
        %
        auxe02 = ibar(7:9,13:15)-auxa2*ibar(1:3,13:15)-auxb1*auxe01;
        auxe03 = ibar(10:12,13:15)-auxa3*ibar(1:3,13:15)-auxb2*auxe01-auxc1*auxe02;
        auxe = inv3by3(ibar(13:15,13:15)-auxa4*ibar(1:3,13:15)-auxb3*auxe01-auxc2*auxe02-	auxd1*auxe03);
        auxe1 = (ibar(16:18,13:15)-auxa5*ibar(1:3,13:15)-auxb4*auxe01-auxc3*auxe02-auxd2*auxe03)*auxe;
        auxe2 = (ibar(19:21,13:15)-auxa6*ibar(1:3,13:15)-auxb5*auxe01-auxc4*auxe02-auxd3*auxe03)*auxe;
        auxe3 = (ibar(22:24,13:15)-auxa7*ibar(1:3,13:15)-auxb6*auxe01-auxc5*auxe02-auxd4*auxe03)*auxe;
        %
        auxf02 = ibar(7:9,16:18)-auxa2*ibar(1:3,16:18)-auxb1*auxf01;
        auxf03 = ibar(10:12,16:18)-auxa3*ibar(1:3,16:18)-auxb2*auxf01-auxc1*auxf02;
        auxf04 = ibar(13:15,16:18)-auxa4*ibar(1:3,16:18)-auxb3*auxf01-auxc2*auxf02-auxd1*auxf03;
        auxf = inv3by3(ibar(16:18,16:18)-auxa5*ibar(1:3,16:18)-auxb4*auxf01-auxc3*auxf02-auxd2*auxf03-auxe1*auxf04);
        auxf1 = (ibar(19:21,16:18)-auxa6*ibar(1:3,16:18)-auxb5*auxf01-auxc4*auxf02-auxd3*auxf03-auxe2*auxf04)*auxf;
        auxf2 = (ibar(22:24,16:18)-auxa7*ibar(1:3,16:18)-auxb6*auxf01-auxc5*auxf02-auxd4*auxf03-auxe3*auxf04)*auxf;
        %
        auxg02 = ibar(7:9,19:21)-auxa2*ibar(1:3,19:21)-auxb1*auxg01;
        auxg03 = ibar(10:12,19:21)-auxa3*ibar(1:3,19:21)-auxb2*auxg01-auxc1*auxg02;
        auxg04 = ibar(13:15,19:21)-auxa4*ibar(1:3,19:21)-auxb3*auxg01-auxc2*auxg02-auxd1*auxg03;
        auxg05 = ibar(16:18,19:21)-auxa5*ibar(1:3,19:21)-auxb4*auxg01-auxc3*auxg02-auxd2*auxg03-auxe1*auxg04;
        auxg = inv3by3(ibar(19:21,19:21)-auxa6*ibar(1:3,19:21)-auxb5*auxg01-auxc4*auxg02-auxd3*auxg03-auxe2*auxg04-auxf1*auxg05);
        auxg1 = (ibar(22:24,19:21)-auxa7*ibar(1:3,19:21)-auxb6*auxg01-auxc5*auxg02-auxd4*auxg03-auxe3*auxg04-auxf2*auxg05)*auxg;
        %
        auxh02 = ibar(7:9,22:24)-auxa2*ibar(1:3,22:24)-auxb1*auxh01;
        auxh03 = ibar(10:12,22:24)-auxa3*ibar(1:3,22:24)-auxb2*auxh01-auxc1*auxh02;
        auxh04 = ibar(13:15,22:24)-auxa4*ibar(1:3,22:24)-auxb3*auxh01-auxc2*auxh02-auxd1*auxh03;
        auxh05 = ibar(16:18,22:24)-auxa5*ibar(1:3,22:24)-auxb4*auxh01-auxc3*auxh02-auxd2*auxh03-auxe1*auxh04;
        auxh06 = ibar(19:21,22:24)-auxa6*ibar(1:3,22:24)-auxb5*auxh01-auxc4*auxh02-auxd3*auxh03-auxe2*auxh04-auxf1*auxh05;
        auxh = inv3by3(ibar(22:24,22:24)-auxa7*ibar(1:3,22:24)-auxb6*auxh01-auxc5*auxh02-auxd4*auxh03-auxe3*auxh04-auxf2*auxh05-auxg1*auxh06);
        %
        auxp1 = psibar(4:6)-auxa1*psibar(1:3);
        auxp2 = psibar(7:9)-auxa2*psibar(1:3)-auxb1*auxp1;
        auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
        auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 3
    elseif jj<=nus3
        [ps3,ths3,bs3,hIs3,hFs3,hGs3,hetats3,hetabs3]=ddq_tree_eff_common(q(nls3:nus3), ni(3), b(nls3:nus3), th(nls3:nus3), r(nls3:nus3), dx(nls3:nus3), dy(nls3:nus3), dz(nls3:nus3), m(nls3:nus3), Icxx(nls3:nus3), Icyy(nls3:nus3), Iczz(nls3:nus3), Icxy(nls3:nus3), Icyz(nls3:nus3), Iczx(nls3:nus3));
        phihat(nls3:nus3) = ddq_tree_eff_variable(ni(3), alp(nls3:nus3), a(nls3:nus3),[0;cumsum(ones(ni(3)-1,1))], r(nls3:nus3), ps3, ths3, bs3, hIs3, hFs3, hGs3, hetats3, hetabs3, phil(nls3:nus3));
        for ii = 7:1:9
    	    ihatTr(nls3:nus3,ii) = ddq_tree_eff_variable(ni(3), alp(nls3:nus3), a(nls3:nus3),[0;cumsum(ones(ni(3)-1,1))], r(nls3:nus3), ps3, ths3, bs3, hIs3, hFs3, hGs3, hetats3, hetabs3, J(ii,nls3:nus3)');
        end
    	% Ibar
    	ibar(7:9,7:9) = J(7:9,nls3:nus3)*ihatTr(nls3:nus3,7:9)+ J(7:9, nls9:nus9)*ihatTr(nls9:nus9,7:9);
        % Psibar
    	psibarAux1(7:9) = - J(7:9,nls3:nus3)*phihat(nls3:nus3)-J(7:9,nls9:nus9)*phihat(nls9:nus9);
    	psibar(7:9) = psibarAux1(7:9)+psibarAux2(7:9);
    	% Computation of lambda
        auxc  = inv3by3(ibar(7:9,7:9)-auxa2*ibar(1:3,7:9)-auxb1*auxc0);
        auxc1 = (ibar(10:12,7:9)-auxa3*ibar(1:3,7:9)-auxb2*auxc0)*auxc;
        auxc2 = (ibar(13:15,7:9)-auxa4*ibar(1:3,7:9)-auxb3*auxc0)*auxc;
        auxc3 = (ibar(16:18,7:9)-auxa5*ibar(1:3,7:9)-auxb4*auxc0)*auxc;
        auxc4 = (ibar(19:21,7:9)-auxa6*ibar(1:3,7:9)-auxb5*auxc0)*auxc;
        auxc5 = (ibar(22:24,7:9)-auxa7*ibar(1:3,7:9)-auxb6*auxc0)*auxc;
        %
        auxd = inv3by3(ibar(10:12,10:12)-auxa3*ibar(1:3,10:12)-auxb2*auxd01-auxc1*auxd02);
        auxd1 = (ibar(13:15,10:12)-auxa4*ibar(1:3,10:12)-auxb3*auxd01-auxc2*auxd02)*auxd;
        auxd2 = (ibar(16:18,10:12)-auxa5*ibar(1:3,10:12)-auxb4*auxd01-auxc3*auxd02)*auxd;
        auxd3 = (ibar(19:21,10:12)-auxa6*ibar(1:3,10:12)-auxb5*auxd01-auxc4*auxd02)*auxd;
        auxd4 = (ibar(22:24,10:12)-auxa7*ibar(1:3,10:12)-auxb6*auxd01-auxc5*auxd02)*auxd;
        %
        auxe03 = ibar(10:12,13:15)-auxa3*ibar(1:3,13:15)-auxb2*auxe01-auxc1*auxe02;
        auxe = inv3by3(ibar(13:15,13:15)-auxa4*ibar(1:3,13:15)-auxb3*auxe01-auxc2*auxe02-	auxd1*auxe03);
        auxe1 = (ibar(16:18,13:15)-auxa5*ibar(1:3,13:15)-auxb4*auxe01-auxc3*auxe02-auxd2*auxe03)*auxe;
        auxe2 = (ibar(19:21,13:15)-auxa6*ibar(1:3,13:15)-auxb5*auxe01-auxc4*auxe02-auxd3*auxe03)*auxe;
        auxe3 = (ibar(22:24,13:15)-auxa7*ibar(1:3,13:15)-auxb6*auxe01-auxc5*auxe02-auxd4*auxe03)*auxe;
        %
        auxf03 = ibar(10:12,16:18)-auxa3*ibar(1:3,16:18)-auxb2*auxf01-auxc1*auxf02;
        auxf04 = ibar(13:15,16:18)-auxa4*ibar(1:3,16:18)-auxb3*auxf01-auxc2*auxf02-auxd1*auxf03;
        auxf = inv3by3(ibar(16:18,16:18)-auxa5*ibar(1:3,16:18)-auxb4*auxf01-auxc3*auxf02-auxd2*auxf03-auxe1*auxf04);
        auxf1 = (ibar(19:21,16:18)-auxa6*ibar(1:3,16:18)-auxb5*auxf01-auxc4*auxf02-auxd3*auxf03-auxe2*auxf04)*auxf;
        auxf2 = (ibar(22:24,16:18)-auxa7*ibar(1:3,16:18)-auxb6*auxf01-auxc5*auxf02-auxd4*auxf03-auxe3*auxf04)*auxf;
        %
        auxg03 = ibar(10:12,19:21)-auxa3*ibar(1:3,19:21)-auxb2*auxg01-auxc1*auxg02;
        auxg04 = ibar(13:15,19:21)-auxa4*ibar(1:3,19:21)-auxb3*auxg01-auxc2*auxg02-auxd1*auxg03;
        auxg05 = ibar(16:18,19:21)-auxa5*ibar(1:3,19:21)-auxb4*auxg01-auxc3*auxg02-auxd2*auxg03-auxe1*auxg04;
        auxg = inv3by3(ibar(19:21,19:21)-auxa6*ibar(1:3,19:21)-auxb5*auxg01-auxc4*auxg02-auxd3*auxg03-auxe2*auxg04-auxf1*auxg05);
        auxg1 = (ibar(22:24,19:21)-auxa7*ibar(1:3,19:21)-auxb6*auxg01-auxc5*auxg02-auxd4*auxg03-auxe3*auxg04-auxf2*auxg05)*auxg;
        %
        auxh03 = ibar(10:12,22:24)-auxa3*ibar(1:3,22:24)-auxb2*auxh01-auxc1*auxh02;
        auxh04 = ibar(13:15,22:24)-auxa4*ibar(1:3,22:24)-auxb3*auxh01-auxc2*auxh02-auxd1*auxh03;
        auxh05 = ibar(16:18,22:24)-auxa5*ibar(1:3,22:24)-auxb4*auxh01-auxc3*auxh02-auxd2*auxh03-auxe1*auxh04;
        auxh06 = ibar(19:21,22:24)-auxa6*ibar(1:3,22:24)-auxb5*auxh01-auxc4*auxh02-auxd3*auxh03-auxe2*auxh04-auxf1*auxh05;
        auxh = inv3by3(ibar(22:24,22:24)-auxa7*ibar(1:3,22:24)-auxb6*auxh01-auxc5*auxh02-auxd4*auxh03-auxe3*auxh04-auxf2*auxh05-auxg1*auxh06);
        %
        auxp2 = psibar(7:9)-auxa2*psibar(1:3)-auxb1*auxp1;
        auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
        auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 4
    elseif jj<=nus4
        [ps4,ths4,bs4,hIs4,hFs4,hGs4,hetats4,hetabs4]=ddq_tree_eff_common(q(nls4:nus4), ni(4), b(nls4:nus4), th(nls4:nus4), r(nls4:nus4), dx(nls4:nus4), dy(nls4:nus4), dz(nls4:nus4), m(nls4:nus4), Icxx(nls4:nus4), Icyy(nls4:nus4), Iczz(nls4:nus4), Icxy(nls4:nus4), Icyz(nls4:nus4), Iczx(nls4:nus4));
    	phihat(nls4:nus4) = ddq_tree_eff_variable(ni(4), alp(nls4:nus4), a(nls4:nus4),[0;cumsum(ones(ni(4)-1,1))], r(nls4:nus4), ps4, ths4, bs4, hIs4, hFs4, hGs4, hetats4, hetabs4, phil(nls4:nus4));
        for ii = 10:1:12
    	    ihatTr(nls4:nus4,ii) = ddq_tree_eff_variable(ni(4), alp(nls4:nus4), a(nls4:nus4),[0;cumsum(ones(ni(4)-1,1))], r(nls4:nus4), ps4, ths4, bs4, hIs4, hFs4, hGs4, hetats4, hetabs4, J(ii,nls4:nus4)');
        end
    	% Ibar
    	ibar(10:12,10:12) = J(10:12,nls4:nus4)*ihatTr(nls4:nus4,10:12) + J(10:12,nls9:nus9)*ihatTr(nls9:nus9, 10:12);
        % Psibar
        psibarAux1(10:12) = -J(10:12,nls4:nus4)*phihat(nls4:nus4)-J(10:12,nls9:nus9)*phihat(nls9:nus9);
        psibar(10:12) = psibarAux1(10:12)+psibarAux2(10:12);
        % Computation of lambda
        auxd = inv3by3(ibar(10:12,10:12)-auxa3*ibar(1:3,10:12)-auxb2*auxd01-auxc1*auxd02);
        auxd1 = (ibar(13:15,10:12)-auxa4*ibar(1:3,10:12)-auxb3*auxd01-auxc2*auxd02)*auxd;
        auxd2 = (ibar(16:18,10:12)-auxa5*ibar(1:3,10:12)-auxb4*auxd01-auxc3*auxd02)*auxd;
        auxd3 = (ibar(19:21,10:12)-auxa6*ibar(1:3,10:12)-auxb5*auxd01-auxc4*auxd02)*auxd;
        auxd4 = (ibar(22:24,10:12)-auxa7*ibar(1:3,10:12)-auxb6*auxd01-auxc5*auxd02)*auxd;
        %
        auxe = inv3by3(ibar(13:15,13:15)-auxa4*ibar(1:3,13:15)-auxb3*auxe01-auxc2*auxe02-auxd1*auxe03);
        auxe1 = (ibar(16:18,13:15)-auxa5*ibar(1:3,13:15)-auxb4*auxe01-auxc3*auxe02-auxd2*auxe03)*auxe;
        auxe2 = (ibar(19:21,13:15)-auxa6*ibar(1:3,13:15)-auxb5*auxe01-auxc4*auxe02-auxd3*auxe03)*auxe;
        auxe3 = (ibar(22:24,13:15)-auxa7*ibar(1:3,13:15)-auxb6*auxe01-auxc5*auxe02-auxd4*auxe03)*auxe;
        %
        auxf04 = ibar(13:15,16:18)-auxa4*ibar(1:3,16:18)-auxb3*auxf01-auxc2*auxf02-auxd1*auxf03;
        auxf = inv3by3(ibar(16:18,16:18)-auxa5*ibar(1:3,16:18)-auxb4*auxf01-auxc3*auxf02-auxd2*auxf03-auxe1*auxf04);
        auxf1 = (ibar(19:21,16:18)-auxa6*ibar(1:3,16:18)-auxb5*auxf01-auxc4*auxf02-auxd3*auxf03-auxe2*auxf04)*auxf;
        auxf2 = (ibar(22:24,16:18)-auxa7*ibar(1:3,16:18)-auxb6*auxf01-auxc5*auxf02-auxd4*auxf03-auxe3*auxf04)*auxf;
        %
        auxg04 = ibar(13:15,19:21)-auxa4*ibar(1:3,19:21)-auxb3*auxg01-auxc2*auxg02-auxd1*auxg03;
        auxg05 = ibar(16:18,19:21)-auxa5*ibar(1:3,19:21)-auxb4*auxg01-auxc3*auxg02-auxd2*auxg03-auxe1*auxg04;
        auxg = inv3by3(ibar(19:21,19:21)-auxa6*ibar(1:3,19:21)-auxb5*auxg01-auxc4*auxg02-auxd3*auxg03-auxe2*auxg04-auxf1*auxg05);
        auxg1 = (ibar(22:24,19:21)-auxa7*ibar(1:3,19:21)-auxb6*auxg01-auxc5*auxg02-auxd4*auxg03-auxe3*auxg04-auxf2*auxg05)*auxg;
        %
        auxh04 = ibar(13:15,22:24)-auxa4*ibar(1:3,22:24)-auxb3*auxh01-auxc2*auxh02-auxd1*auxh03;
        auxh05 = ibar(16:18,22:24)-auxa5*ibar(1:3,22:24)-auxb4*auxh01-auxc3*auxh02-auxd2*auxh03-auxe1*auxh04;
        auxh06 = ibar(19:21,22:24)-auxa6*ibar(1:3,22:24)-auxb5*auxh01-auxc4*auxh02-auxd3*auxh03-auxe2*auxh04-auxf1*auxh05;
        auxh = inv3by3(ibar(22:24,22:24)-auxa7*ibar(1:3,22:24)-auxb6*auxh01-auxc5*auxh02-auxd4*auxh03-auxe3*auxh04-auxf2*auxh05-auxg1*auxh06);
        %
        auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
        auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 5
    elseif jj<=nus5
        [ps5,ths5,bs5,hIs5,hFs5,hGs5,hetats5,hetabs5]=ddq_tree_eff_common(q(nls5:nus5), ni(5), b(nls5:nus5), th(nls5:nus5), r(nls5:nus5), dx(nls5:nus5), dy(nls5:nus5), dz(nls5:nus5), m(nls5:nus5), Icxx(nls5:nus5), Icyy(nls5:nus5), Iczz(nls5:nus5), Icxy(nls5:nus5), Icyz(nls5:nus5), Iczx(nls5:nus5));
    	phihat(nls5:nus5) = ddq_tree_eff_variable(ni(5), alp(nls5:nus5), a(nls5:nus5),[0;cumsum(ones(ni(5)-1,1))], r(nls5:nus5), ps5, ths5, bs5, hIs5, hFs5, hGs5, hetats5, hetabs5, phil(nls5:nus5));
        for ii = 13:1:15
    	    ihatTr(nls5:nus5,ii) = ddq_tree_eff_variable(ni(5), alp(nls5:nus5), a(nls5:nus5),[0;cumsum(ones(ni(5)-1,1))], r(nls5:nus5), ps5, ths5, bs5, hIs5, hFs5, hGs5, hetats5, hetabs5, J(ii,nls5:nus5)');
        end
    	% Ibar
    	ibar(13:15,13:15) = J(13:15,nls5:nus5)*ihatTr(nls5:nus5,13:15) + J(13:15, nls9:nus9)*ihatTr(nls9:nus9, 13:15);
        % Psibar
        psibarAux1(13:15) = -J(13:15,nls5:nus5)*phihat(nls5:nus5)-J(13:15,nls9:nus9)*phihat(nls9:nus9);
    	psibar(13:15) = psibarAux1(13:15)+psibarAux2(13:15);
    	% Computation of lambda
    	auxe = inv3by3(ibar(13:15,13:15)-auxa4*ibar(1:3,13:15)-auxb3*auxe01-auxc2*auxe02-auxd1*auxe03);
        auxe1 = (ibar(16:18,13:15)-auxa5*ibar(1:3,13:15)-auxb4*auxe01-auxc3*auxe02-auxd2*auxe03)*auxe;
        auxe2 = (ibar(19:21,13:15)-auxa6*ibar(1:3,13:15)-auxb5*auxe01-auxc4*auxe02-auxd3*auxe03)*auxe;
        auxe3 = (ibar(22:24,13:15)-auxa7*ibar(1:3,13:15)-auxb6*auxe01-auxc5*auxe02-auxd4*auxe03)*auxe;
        %
        auxf = inv3by3(ibar(16:18,16:18)-auxa5*ibar(1:3,16:18)-auxb4*auxf01-auxc3*auxf02-auxd2*auxf03-auxe1*auxf04);
        auxf1 = (ibar(19:21,16:18)-auxa6*ibar(1:3,16:18)-auxb5*auxf01-auxc4*auxf02-auxd3*auxf03-auxe2*auxf04)*auxf;
        auxf2 = (ibar(22:24,16:18)-auxa7*ibar(1:3,16:18)-auxb6*auxf01-auxc5*auxf02-auxd4*auxf03-auxe3*auxf04)*auxf;
        %
        auxg05 = ibar(16:18,19:21)-auxa5*ibar(1:3,19:21)-auxb4*auxg01-auxc3*auxg02-auxd2*auxg03-auxe1*auxg04;
        auxg = inv3by3(ibar(19:21,19:21)-auxa6*ibar(1:3,19:21)-auxb5*auxg01-auxc4*auxg02-auxd3*auxg03-auxe2*auxg04-auxf1*auxg05);
        auxg1 = (ibar(22:24,19:21)-auxa7*ibar(1:3,19:21)-auxb6*auxg01-auxc5*auxg02-auxd4*auxg03-auxe3*auxg04-auxf2*auxg05)*auxg;
        %
        auxh05 = ibar(16:18,22:24)-auxa5*ibar(1:3,22:24)-auxb4*auxh01-auxc3*auxh02-auxd2*auxh03-auxe1*auxh04;
        auxh06 = ibar(19:21,22:24)-auxa6*ibar(1:3,22:24)-auxb5*auxh01-auxc4*auxh02-auxd3*auxh03-auxe2*auxh04-auxf1*auxh05;
        auxh = inv3by3(ibar(22:24,22:24)-auxa7*ibar(1:3,22:24)-auxb6*auxh01-auxc5*auxh02-auxd4*auxh03-auxe3*auxh04-auxf2*auxh05-auxg1*auxh06);
        %
        auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 6
    elseif jj<=nus6
        [ps6,ths6,bs6,hIs6,hFs6,hGs6,hetats6,hetabs6]=ddq_tree_eff_common(q(nls6:nus6), ni(6), b(nls6:nus6), th(nls6:nus6), r(nls6:nus6), dx(nls6:nus6), dy(nls6:nus6), dz(nls6:nus6), m(nls6:nus6), Icxx(nls6:nus6), Icyy(nls6:nus6), Iczz(nls6:nus6), Icxy(nls6:nus6), Icyz(nls6:nus6), Iczx(nls6:nus6));
        phihat(nls6:nus6) = ddq_tree_eff_variable(ni(6), alp(nls6:nus6), a(nls6:nus6),[0;cumsum(ones(ni(6)-1,1))], r(nls6:nus6), ps6, ths6, bs6, hIs6, hFs6, hGs6, hetats6, hetabs6, phil(nls6:nus6));
        for ii = 16:1:18
            ihatTr(nls6:nus6,ii) = ddq_tree_eff_variable(ni(6), alp(nls6:nus6), a(nls6:nus6),[0;cumsum(ones(ni(6)-1,1))], r(nls6:nus6), ps6, ths6, bs6, hIs6, hFs6, hGs6, hetats6, hetabs6, J(ii,nls6:nus6)');
        end
        % Ibar
    	ibar(16:18,16:18) = J(16:18,nls6:nus6)*ihatTr(nls6:nus6,16:18) + J(16:18, nls9:nus9)*ihatTr(nls9:nus9, 16:18);
        % Psibar
    	psibarAux1(16:18) = -J(16:18,nls6:nus6)*phihat(nls6:nus6)-J(16:18,nls9:nus9)*phihat(nls9:nus9);
    	psibar(16:18) = psibarAux1(16:18)+psibarAux2(16:18);
    	% Computation of lambda
    	auxf = inv3by3(ibar(16:18,16:18)-auxa5*ibar(1:3,16:18)-auxb4*auxf01-auxc3*auxf02-auxd2*auxf03-auxe1*auxf04);
        auxf1 = (ibar(19:21,16:18)-auxa6*ibar(1:3,16:18)-auxb5*auxf01-auxc4*auxf02-auxd3*auxf03-auxe2*auxf04)*auxf;
        auxf2 = (ibar(22:24,16:18)-auxa7*ibar(1:3,16:18)-auxb6*auxf01-auxc5*auxf02-auxd4*auxf03-auxe3*auxf04)*auxf;
        %
        auxg = inv3by3(ibar(19:21,19:21)-auxa6*ibar(1:3,19:21)-auxb5*auxg01-auxc4*auxg02-auxd3*auxg03-auxe2*auxg04-auxf1*auxg05);
        auxg1 = (ibar(22:24,19:21)-auxa7*ibar(1:3,19:21)-auxb6*auxg01-auxc5*auxg02-auxd4*auxg03-auxe3*auxg04-auxf2*auxg05)*auxg;
        %
        auxh06 = ibar(19:21,22:24)-auxa6*ibar(1:3,22:24)-auxb5*auxh01-auxc4*auxh02-auxd3*auxh03-auxe2*auxh04-auxf1*auxh05;
        auxh = inv3by3(ibar(22:24,22:24)-auxa7*ibar(1:3,22:24)-auxb6*auxh01-auxc5*auxh02-auxd4*auxh03-auxe3*auxh04-auxf2*auxh05-auxg1*auxh06);
        %
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 7
    elseif jj<=nus7
        [ps7,ths7,bs7,hIs7,hFs7,hGs7,hetats7,hetabs7]=ddq_tree_eff_common(q(nls7:nus7), ni(7), b(nls7:nus7), th(nls7:nus7), r(nls7:nus7), dx(nls7:nus7), dy(nls7:nus7), dz(nls7:nus7), m(nls7:nus7), Icxx(nls7:nus7), Icyy(nls7:nus7), Iczz(nls7:nus7), Icxy(nls7:nus7), Icyz(nls7:nus7), Iczx(nls7:nus7));
        phihat(nls7:nus7) = ddq_tree_eff_variable(ni(7), alp(nls7:nus7), a(nls7:nus7),[0;cumsum(ones(ni(7)-1,1))], r(nls7:nus7), ps7, ths7, bs7, hIs7, hFs7, hGs7, hetats7, hetabs7, phil(nls7:nus7));
        for ii = 19:1:21
            ihatTr(nls7:nus7,ii) = ddq_tree_eff_variable(ni(7), alp(nls7:nus7), a(nls7:nus7),[0;cumsum(ones(ni(7)-1,1))], r(nls7:nus7), ps7, ths7, bs7, hIs7, hFs7, hGs7, hetats7, hetabs7, J(ii,nls7:nus7)');
        end
        % Ibar
    	ibar(19:21,19:21) = J(19:21,nls7:nus7)*ihatTr(nls7:nus7,19:21) + J(19:21, nls9:nus9)*ihatTr(nls9:nus9, 19:21);
        % Psibar
    	psibarAux1(19:21) = -J(19:21,nls7:nus7)*phihat(nls7:nus7)-J(19:21,nls9:nus9)*phihat(nls9:nus9);
    	psibar(19:21) = psibarAux1(19:21)+psibarAux2(19:21);
    	% Computation of lambda
    	auxg = inv3by3(ibar(19:21,19:21)-auxa6*ibar(1:3,19:21)-auxb5*auxg01-auxc4*auxg02-auxd3*auxg03-auxe2*auxg04-auxf1*auxg05);
        auxg1 = (ibar(22:24,19:21)-auxa7*ibar(1:3,19:21)-auxb6*auxg01-auxc5*auxg02-auxd4*auxg03-auxe3*auxg04-auxf2*auxg05)*auxg;
        %
        auxh = inv3by3(ibar(22:24,22:24)-auxa7*ibar(1:3,22:24)-auxb6*auxh01-auxc5*auxh02-auxd4*auxh03-auxe3*auxh04-auxf2*auxh05-auxg1*auxh06);
        %
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 8
    elseif jj<=nus8
        [ps8,ths8,bs8,hIs8,hFs8,hGs8,hetats8,hetabs8]=ddq_tree_eff_common(q(nls8:nus8), ni(8), b(nls8:nus8), th(nls8:nus8), r(nls8:nus8), dx(nls8:nus8), dy(nls8:nus8), dz(nls8:nus8), m(nls8:nus8), Icxx(nls8:nus8), Icyy(nls8:nus8), Iczz(nls8:nus8), Icxy(nls8:nus8), Icyz(nls8:nus8), Iczx(nls8:nus8));
        phihat(nls8:nus8) = ddq_tree_eff_variable(ni(8), alp(nls8:nus8), a(nls8:nus8),[0;cumsum(ones(ni(8)-1,1))], r(nls8:nus8), ps8, ths8, bs8, hIs8, hFs8, hGs8, hetats8, hetabs8, phil(nls8:nus8));
        for ii = 22:1:24
            ihatTr(nls8:nus8,ii) = ddq_tree_eff_variable(ni(8), alp(nls8:nus8), a(nls8:nus8),[0;cumsum(ones(ni(8)-1,1))], r(nls8:nus8), ps8, ths8, bs8, hIs8, hFs8, hGs8, hetats8, hetabs8, J(ii,nls8:nus8)');
        end
        % Ibar
    	ibar(22:24,22:24) = J(22:24,nls8:nus8)*ihatTr(nls8:nus8,22:24) + J(22:24, nls9:nus9)*ihatTr(nls9:nus9, 22:24);
        % Psibar
    	psibarAux1(22:24) = -J(22:24,nls8:nus8)*phihat(nls8:nus8)-J(22:24,nls9:nus9)*phihat(nls9:nus9);
    	psibar(22:24) = psibarAux1(22:24)+psibarAux2(22:24);
    	% Computation of lambda
    	auxh = inv3by3(ibar(22:24,22:24)-auxa7*ibar(1:3,22:24)-auxb6*auxh01-auxc5*auxh02-auxd4*auxh03-auxe3*auxh04-auxf2*auxh05-auxg1*auxh06);
    	auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    %sub-system 9
    else
        phihat(nls9:nls9p2) = invmp.*phil(nls9:nls9p2);
        phihat(nls9p3:nus9) = I92*phil(nls9p3:nus9);
        for ii = 1:1:24
            ihatTr(nls9:nls9p2,ii) = invmp.*J(ii,nls9:nls9p2)';
            ihatTr(nls9p3:nus9,ii) = I92*J(ii,nls9p3:nus9)';
        end
        % Ibar
        % Rows 1-3
        ibar(1:3,1:3) = J(1:3,nls1:nus1)*ihatTr(nls1:nus1,1:3)+J(1:3,nls9:nus9)*ihatTr(nls9:nus9,1:3);
        ibar(1:3,4:24) = J(1:3,nls9:nus9)*ihatTr(nls9:nus9,4:24);
        % Rows 4-6
        ibar(4:6,1:3) = ibar(1:3,4:6)';
        ibar(4:6,4:6) = J(4:6,nls2:nus2)*ihatTr(nls2:nus2,4:6)+ J(4:6,nls9:nus9)*ihatTr(nls9:nus9,4:6);
        ibar(4:6,7:24) = J(4:6,nls9:nus9)*ihatTr(nls9:nus9,7:24);
        % Rows 7-9
        ibar(7:9,1:6) = ibar(1:6,7:9)';
        ibar(7:9,7:9) = J(7:9,nls3:nus3)*ihatTr(nls3:nus3,7:9)+ J(7:9, nls9:nus9)*ihatTr(nls9:nus9,7:9);
        ibar(7:9,10:24) = J(7:9,nls9:nus9)*ihatTr(nls9:nus9,10:24);
        % Rows 10-12
        ibar(10:12,1:9) = ibar(1:9,10:12)';
        ibar(10:12,10:12) = J(10:12,nls4:nus4)*ihatTr(nls4:nus4,10:12) + J(10:12,nls9:nus9)*ihatTr(nls9:nus9, 10:12);
        ibar(10:12,13:24) = J(10:12,nls9:nus9)*ihatTr(nls9:nus9,13:24);
        % Rows 13-15
        ibar(13:15,1:12) = ibar(1:12,13:15)';
        ibar(13:15,13:15) = J(13:15,nls5:nus5)*ihatTr(nls5:nus5,13:15) + J(13:15, nls9:nus9)*ihatTr(nls9:nus9, 13:15);
        ibar(13:15,16:24) = J(13:15,nls9:nus9)*ihatTr(nls9:nus9,16:24);
        % Rows 16-18
        ibar(16:18,1:15) = ibar(1:15,16:18)';
        ibar(16:18,16:18) = J(16:18,nls6:nus6)*ihatTr(nls6:nus6,16:18) + J(16:18, nls9:nus9)*ihatTr(nls9:nus9, 16:18);
        ibar(16:18,19:24) = J(16:18,nls9:nus9)*ihatTr(nls9:nus9,19:24);
        % Rows 19-21
        ibar(19:21,1:18) = ibar(1:18,19:21)';
        ibar(19:21,19:21) = J(19:21,nls7:nus7)*ihatTr(nls7:nus7,19:21) + J(19:21, nls9:nus9)*ihatTr(nls9:nus9, 19:21);
        ibar(19:21,22:24) = J(19:21,nls9:nus9)*ihatTr(nls9:nus9,22:24);
        % Rows 22-24
        ibar(22:24,1:21) = ibar(1:21,22:24)';
        ibar(22:24,22:24) = J(22:24,nls8:nus8)*ihatTr(nls8:nus8,22:24) + J(22:24, nls9:nus9)*ihatTr(nls9:nus9, 22:24);
      	% Psibar
    	psibarAux1(1:3) = -J(1:3,nls1:nus1)*phihat(nls1:nus1)-J(1:3,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(4:6) = -J(4:6,nls2:nus2)*phihat(nls2:nus2)-J(4:6,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(7:9) = -J(7:9,nls3:nus3)*phihat(nls3:nus3)-J(7:9,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(10:12) = -J(10:12,nls4:nus4)*phihat(nls4:nus4)-J(10:12,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(13:15) = -J(13:15,nls5:nus5)*phihat(nls5:nus5)-J(13:15,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(16:18) = -J(16:18,nls6:nus6)*phihat(nls6:nus6)-J(16:18,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(19:21) = -J(19:21,nls7:nus7)*phihat(nls7:nus7)-J(19:21,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(22:24) = -J(22:24,nls8:nus8)*phihat(nls8:nus8)-J(22:24,nls9:nus9)*phihat(nls9:nus9);
    	psibar = psibarAux1+psibarAux2;
    	% Computation of lambda
        auxa  = inv3by3(ibar(1:3,1:3));
    	auxa1 = ibar(4:6,1:3)*auxa;
    	auxa2 = ibar(7:9,1:3)*auxa;
    	auxa3 = ibar(10:12,1:3)*auxa;
    	auxa4 = ibar(13:15,1:3)*auxa;
    	auxa5 = ibar(16:18,1:3)*auxa;
    	auxa6 = ibar(19:21,1:3)*auxa;
    	auxa7 = ibar(22:24,1:3)*auxa;
    	%
    	auxb  = inv3by3(ibar(4:6,4:6)-auxa1*ibar(1:3,4:6));
    	auxb1 = (ibar(7:9,4:6)-auxa2*ibar(1:3,4:6))*auxb;
    	auxb2 = (ibar(10:12,4:6)-auxa3*ibar(1:3,4:6))*auxb;
    	auxb3 = (ibar(13:15,4:6)-auxa4*ibar(1:3,4:6))*auxb;
    	auxb4 = (ibar(16:18,4:6)-auxa5*ibar(1:3,4:6))*auxb;
    	auxb5 = (ibar(19:21,4:6)-auxa6*ibar(1:3,4:6))*auxb;
    	auxb6 = (ibar(22:24,4:6)-auxa7*ibar(1:3,4:6))*auxb;
    	%
    	auxc0 = ibar(4:6,7:9)-auxa1*ibar(1:3,7:9);
    	auxc  = inv3by3(ibar(7:9,7:9)-auxa2*ibar(1:3,7:9)-auxb1*auxc0);
    	auxc1 = (ibar(10:12,7:9)-auxa3*ibar(1:3,7:9)-auxb2*auxc0)*auxc;
    	auxc2 = (ibar(13:15,7:9)-auxa4*ibar(1:3,7:9)-auxb3*auxc0)*auxc;
    	auxc3 = (ibar(16:18,7:9)-auxa5*ibar(1:3,7:9)-auxb4*auxc0)*auxc;
    	auxc4 = (ibar(19:21,7:9)-auxa6*ibar(1:3,7:9)-auxb5*auxc0)*auxc;
    	auxc5 = (ibar(22:24,7:9)-auxa7*ibar(1:3,7:9)-auxb6*auxc0)*auxc;
    	%
    	auxd01 = ibar(4:6,10:12)-auxa1*ibar(1:3,10:12);
    	auxd02 = ibar(7:9,10:12)-auxa2*ibar(1:3,10:12)-auxb1*auxd01; 
    	auxd = inv3by3(ibar(10:12,10:12)-auxa3*ibar(1:3,10:12)-auxb2*auxd01-auxc1*auxd02);
    	auxd1 = (ibar(13:15,10:12)-auxa4*ibar(1:3,10:12)-auxb3*auxd01-auxc2*auxd02)*auxd;
    	auxd2 = (ibar(16:18,10:12)-auxa5*ibar(1:3,10:12)-auxb4*auxd01-auxc3*auxd02)*auxd;
    	auxd3 = (ibar(19:21,10:12)-auxa6*ibar(1:3,10:12)-auxb5*auxd01-auxc4*auxd02)*auxd;
    	auxd4 = (ibar(22:24,10:12)-auxa7*ibar(1:3,10:12)-auxb6*auxd01-auxc5*auxd02)*auxd;
    	%
    	auxe01 = ibar(4:6,13:15)-auxa1*ibar(1:3,13:15);
    	auxe02 = ibar(7:9,13:15)-auxa2*ibar(1:3,13:15)-auxb1*auxe01;
    	auxe03 = ibar(10:12,13:15)-auxa3*ibar(1:3,13:15)-auxb2*auxe01-auxc1*auxe02;
    	auxe = inv3by3(ibar(13:15,13:15)-auxa4*ibar(1:3,13:15)-auxb3*auxe01-auxc2*auxe02-auxd1*auxe03);
    	auxe1 = (ibar(16:18,13:15)-auxa5*ibar(1:3,13:15)-auxb4*auxe01-auxc3*auxe02-auxd2*auxe03)*auxe;
    	auxe2 = (ibar(19:21,13:15)-auxa6*ibar(1:3,13:15)-auxb5*auxe01-auxc4*auxe02-auxd3*auxe03)*auxe;
    	auxe3 = (ibar(22:24,13:15)-auxa7*ibar(1:3,13:15)-auxb6*auxe01-auxc5*auxe02-auxd4*auxe03)*auxe;
    	%
    	auxf01 = ibar(4:6,16:18)-auxa1*ibar(1:3,16:18);
    	auxf02 = ibar(7:9,16:18)-auxa2*ibar(1:3,16:18)-auxb1*auxf01;
    	auxf03 = ibar(10:12,16:18)-auxa3*ibar(1:3,16:18)-auxb2*auxf01-auxc1*auxf02;
    	auxf04 = ibar(13:15,16:18)-auxa4*ibar(1:3,16:18)-auxb3*auxf01-auxc2*auxf02-auxd1*auxf03;
    	auxf = inv3by3(ibar(16:18,16:18)-auxa5*ibar(1:3,16:18)-auxb4*auxf01-auxc3*auxf02-auxd2*auxf03-auxe1*auxf04);
    	auxf1 = (ibar(19:21,16:18)-auxa6*ibar(1:3,16:18)-auxb5*auxf01-auxc4*auxf02-auxd3*auxf03-auxe2*auxf04)*auxf;
    	auxf2 = (ibar(22:24,16:18)-auxa7*ibar(1:3,16:18)-auxb6*auxf01-auxc5*auxf02-auxd4*auxf03-auxe3*auxf04)*auxf;
    	%
    	auxg01 = ibar(4:6,19:21)-auxa1*ibar(1:3,19:21);
    	auxg02 = ibar(7:9,19:21)-auxa2*ibar(1:3,19:21)-auxb1*auxg01;
    	auxg03 = ibar(10:12,19:21)-auxa3*ibar(1:3,19:21)-auxb2*auxg01-auxc1*auxg02;
    	auxg04 = ibar(13:15,19:21)-auxa4*ibar(1:3,19:21)-auxb3*auxg01-auxc2*auxg02-auxd1*auxg03;
    	auxg05 = ibar(16:18,19:21)-auxa5*ibar(1:3,19:21)-auxb4*auxg01-auxc3*auxg02-auxd2*auxg03-auxe1*auxg04;
    	auxg = inv3by3(ibar(19:21,19:21)-auxa6*ibar(1:3,19:21)-auxb5*auxg01-auxc4*auxg02-auxd3*auxg03-auxe2*auxg04-auxf1*auxg05);
    	auxg1 = (ibar(22:24,19:21)-auxa7*ibar(1:3,19:21)-auxb6*auxg01-auxc5*auxg02-auxd4*auxg03-auxe3*auxg04-auxf2*auxg05)*auxg;
    	%
    	auxh01 = ibar(4:6,22:24)-auxa1*ibar(1:3,22:24);
    	auxh02 = ibar(7:9,22:24)-auxa2*ibar(1:3,22:24)-auxb1*auxh01;
    	auxh03 = ibar(10:12,22:24)-auxa3*ibar(1:3,22:24)-auxb2*auxh01-auxc1*auxh02;
    	auxh04 = ibar(13:15,22:24)-auxa4*ibar(1:3,22:24)-auxb3*auxh01-auxc2*auxh02-auxd1*auxh03;
    	auxh05 = ibar(16:18,22:24)-auxa5*ibar(1:3,22:24)-auxb4*auxh01-auxc3*auxh02-auxd2*auxh03-auxe1*auxh04;
    	auxh06 = ibar(19:21,22:24)-auxa6*ibar(1:3,22:24)-auxb5*auxh01-auxc4*auxh02-auxd3*auxh03-auxe2*auxh04-auxf1*auxh05;
    	auxh = inv3by3(ibar(22:24,22:24)-auxa7*ibar(1:3,22:24)-auxb6*auxh01-auxc5*auxh02-auxd4*auxh03-auxe3*auxh04-auxf2*auxh05-auxg1*auxh06);
    	%
    	auxp1 = psibar(4:6)-auxa1*psibar(1:3);
    	auxp2 = psibar(7:9)-auxa2*psibar(1:3)-auxb1*auxp1;
    	auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
    	auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
    	auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
    	auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
    	auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    end
end

if fl==1
    % sub-system 1
    if jj<=nus1
        phihat(nls1:nus1) = ddq_tree_eff(q(nls1:nus1), ni(1), alp(nls1:nus1), a(nls1:nus1), b(nls1:nus1), th(nls1:nus1), [0;cumsum(ones(ni(1)-1,1))], r(nls1:nus1), dx(nls1:nus1), dy(nls1:nus1), dz(nls1:nus1), m(nls1:nus1), Icxx(nls1:nus1), Icyy(nls1:nus1), Iczz(nls1:nus1), Icxy(nls1:nus1), Icyz(nls1:nus1), Iczx(nls1:nus1), phil(nls1:nus1));
        % Psibar
        psibarAux2(1:3) = -dJ(1:3,nls1:nus1)*dq(nls1:nus1)-dJ(1:3,nls9p3:nus9)*dq(nls9p3:nus9);
        psibarAux1(1:3) = -J(1:3,nls1:nus1)*phihat(nls1:nus1)-J(1:3,nls9:nus9)*phihat(nls9:nus9);
    	psibar(1:3) = psibarAux1(1:3)+psibarAux2(1:3);
        %
        auxp1 = psibar(4:6)-auxa1*psibar(1:3);
        auxp2 = psibar(7:9)-auxa2*psibar(1:3)-auxb1*auxp1;
        auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
        auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 2
    elseif jj<=nus2
        phihat(nls2:nus2) = ddq_tree_eff(q(nls2:nus2), ni(2), alp(nls2:nus2), a(nls2:nus2), b(nls2:nus2), th(nls2:nus2), [0;cumsum(ones(ni(2)-1,1))], r(nls2:nus2), dx(nls2:nus2), dy(nls2:nus2), dz(nls2:nus2), m(nls2:nus2), Icxx(nls2:nus2), Icyy(nls2:nus2), Iczz(nls2:nus2), Icxy(nls2:nus2), Icyz(nls2:nus2), Iczx(nls2:nus2), phil(nls2:nus2));
        % Psibar
        psibarAux1(4:6) = -J(4:6,nls2:nus2)*phihat(nls2:nus2)-J(4:6,nls9:nus9)*phihat(nls9:nus9);
        psibarAux2(4:6) = -dJ(4:6,nls2:nus2)*dq(nls2:nus2)-dJ(4:6,nls9p3:nus9)*dq(nls9p3:nus9);
        psibar(4:6) = psibarAux1(4:6)+psibarAux2(4:6);
        %
        auxp1 = psibar(4:6)-auxa1*psibar(1:3);
        auxp2 = psibar(7:9)-auxa2*psibar(1:3)-auxb1*auxp1;
        auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
        auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 3
    elseif jj<=nus3
        phihat(nls3:nus3) = ddq_tree_eff(q(nls3:nus3), ni(3), alp(nls3:nus3), a(nls3:nus3), b(nls3:nus3), th(nls3:nus3), [0;cumsum(ones(ni(3)-1,1))], r(nls3:nus3), dx(nls3:nus3), dy(nls3:nus3), dz(nls3:nus3), m(nls3:nus3), Icxx(nls3:nus3), Icyy(nls3:nus3), Iczz(nls3:nus3), Icxy(nls3:nus3), Icyz(nls3:nus3), Iczx(nls3:nus3), phil(nls3:nus3));
        psibarAux1(7:9) = -J(7:9,nls3:nus3)*phihat(nls3:nus3)-J(7:9,nls9:nus9)*phihat(nls9:nus9);
        psibarAux2(7:9) = -dJ(7:9,nls3:nus3)*dq(nls3:nus3)-dJ(7:9,nls9p3:nus9)*dq(nls9p3:nus9);
        psibar(7:9) = psibarAux1(7:9)+psibarAux2(7:9);
        %
        auxp2 = psibar(7:9)-auxa2*psibar(1:3)-auxb1*auxp1;
        auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
        auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 4
    elseif jj<=nus4
        phihat(nls4:nus4) = ddq_tree_eff(q(nls4:nus4), ni(4), alp(nls4:nus4), a(nls4:nus4), b(nls4:nus4), th(nls4:nus4), [0;cumsum(ones(ni(4)-1,1))], r(nls4:nus4), dx(nls4:nus4), dy(nls4:nus4), dz(nls4:nus4), m(nls4:nus4), Icxx(nls4:nus4), Icyy(nls4:nus4), Iczz(nls4:nus4), Icxy(nls4:nus4), Icyz(nls4:nus4), Iczx(nls4:nus4), phil(nls4:nus4));
        psibarAux1(10:12) = -J(10:12,nls4:nus4)*phihat(nls4:nus4)-J(10:12,nls9:nus9)*phihat(nls9:nus9);
        psibarAux2(10:12) =-dJ(10:12,nls4:nus4)*dq(nls4:nus4)-dJ(10:12,nls9p3:nus9)*dq(nls9p3:nus9);
        psibar(10:12) = psibarAux1(10:12)+psibarAux2(10:12);
        %
    	auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
        auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 5
    elseif jj<=nus5
        phihat(nls5:nus5) = ddq_tree_eff(q(nls5:nus5), ni(5), alp(nls5:nus5), a(nls5:nus5), b(nls5:nus5), th(nls5:nus5), [0;cumsum(ones(ni(5)-1,1))], r(nls5:nus5), dx(nls5:nus5), dy(nls5:nus5), dz(nls5:nus5), m(nls5:nus5), Icxx(nls5:nus5), Icyy(nls5:nus5), Iczz(nls5:nus5), Icxy(nls5:nus5), Icyz(nls5:nus5), Iczx(nls5:nus5), phil(nls5:nus5));
        psibarAux1(13:15) = -J(13:15,nls5:nus5)*phihat(nls5:nus5)-J(13:15,nls9:nus9)*phihat(nls9:nus9);
        psibarAux2(13:15) =-dJ(13:15,nls5:nus5)*dq(nls5:nus5)-dJ(13:15,nls9p3:nus9)*dq(nls9p3:nus9);
        psibar(13:15) = psibarAux1(13:15)+psibarAux2(13:15);
    	%
    	auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 6
    elseif jj<=nus6
        phihat(nls6:nus6) = ddq_tree_eff(q(nls6:nus6), ni(6), alp(nls6:nus6), a(nls6:nus6), b(nls6:nus6), th(nls6:nus6), [0;cumsum(ones(ni(6)-1,1))], r(nls6:nus6), dx(nls6:nus6), dy(nls6:nus6), dz(nls6:nus6), m(nls6:nus6), Icxx(nls6:nus6), Icyy(nls6:nus6), Iczz(nls6:nus6), Icxy(nls6:nus6), Icyz(nls6:nus6), Iczx(nls6:nus6), phil(nls6:nus6));
        psibarAux1(16:18) = -J(16:18,nls6:nus6)*phihat(nls6:nus6)-J(16:18,nls9:nus9)*phihat(nls9:nus9);
        psibarAux2(16:18) =-dJ(16:18,nls6:nus6)*dq(nls6:nus6)-dJ(16:18,nls9p3:nus9)*dq(nls9p3:nus9);
        psibar(16:18) = psibarAux1(16:18)+psibarAux2(16:18);
        %
    	auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 7
    elseif jj<=nus7
        phihat(nls7:nus7) = ddq_tree_eff(q(nls7:nus7), ni(7), alp(nls7:nus7), a(nls7:nus7), b(nls7:nus7), th(nls7:nus7), [0;cumsum(ones(ni(7)-1,1))], r(nls7:nus7), dx(nls7:nus7), dy(nls7:nus7), dz(nls7:nus7), m(nls7:nus7), Icxx(nls7:nus7), Icyy(nls7:nus7), Iczz(nls7:nus7), Icxy(nls7:nus7), Icyz(nls7:nus7), Iczx(nls7:nus7), phil(nls7:nus7));
        psibarAux1(19:21) = -J(19:21,nls7:nus7)*phihat(nls7:nus7)-J(19:21,nls9:nus9)*phihat(nls9:nus9);
        psibarAux2(19:21) =-dJ(19:21,nls7:nus7)*dq(nls7:nus7)-dJ(19:21,nls9p3:nus9)*dq(nls9p3:nus9);
        psibar(19:21) = psibarAux1(19:21)+psibarAux2(19:21);
        %
    	auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    % sub-system 8
    elseif jj<=nus8
        phihat(nls8:nus8) = ddq_tree_eff(q(nls8:nus8), ni(8), alp(nls8:nus8), a(nls8:nus8), b(nls8:nus8), th(nls8:nus8), [0;cumsum(ones(ni(8)-1,1))], r(nls8:nus8), dx(nls8:nus8), dy(nls8:nus8), dz(nls8:nus8), m(nls8:nus8), Icxx(nls8:nus8), Icyy(nls8:nus8), Iczz(nls8:nus8), Icxy(nls8:nus8), Icyz(nls8:nus8), Iczx(nls8:nus8), phil(nls8:nus8));
        psibarAux1(22:24) = -J(22:24,nls8:nus8)*phihat(nls8:nus8)-J(22:24,nls9:nus9)*phihat(nls9:nus9);
        psibarAux2(22:24) =-dJ(22:24,nls8:nus8)*dq(nls8:nus8)-dJ(22:24,nls9p3:nus9)*dq(nls9p3:nus9);
        psibar(22:24) = psibarAux1(22:24)+psibarAux2(22:24);
        %
    	auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    %sub-system 9
    else
        phihat(nls9:nls9p2) = invmp.*phil(nls9:nls9p2);
    	phihat(nls9p3:nus9) = I92*phil(nls9p3:nus9);
    	psibarAux1(1:3) = -J(1:3,nls1:nus1)*phihat(nls1:nus1)-J(1:3,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(4:6) = -J(4:6,nls2:nus2)*phihat(nls2:nus2)-J(4:6,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(7:9) = -J(7:9,nls3:nus3)*phihat(nls3:nus3)-J(7:9,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(10:12) = -J(10:12,nls4:nus4)*phihat(nls4:nus4)-J(10:12,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(13:15) = -J(13:15,nls5:nus5)*phihat(nls5:nus5)-J(13:15,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(16:18) = -J(16:18,nls6:nus6)*phihat(nls6:nus6)-J(16:18,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(19:21) = -J(19:21,nls7:nus7)*phihat(nls7:nus7)-J(19:21,nls9:nus9)*phihat(nls9:nus9);
        psibarAux1(22:24) = -J(22:24,nls8:nus8)*phihat(nls8:nus8)-J(22:24,nls9:nus9)*phihat(nls9:nus9);
        %
    	psibarAux2(1:3) = -dJ(1:3,nls1:nus1)*dq(nls1:nus1)-dJ(1:3,nls9p3:nus9)*dq(nls9p3:nus9);
    	psibarAux2(4:6) = -dJ(4:6,nls2:nus2)*dq(nls2:nus2)-dJ(4:6,nls9p3:nus9)*dq(nls9p3:nus9);
        psibarAux2(7:9) = -dJ(7:9,nls3:nus3)*dq(nls3:nus3)-dJ(7:9,nls9p3:nus9)*dq(nls9p3:nus9);
        psibarAux2(10:12) =-dJ(10:12,nls4:nus4)*dq(nls4:nus4)-dJ(10:12,nls9p3:nus9)*dq(nls9p3:nus9);
        psibarAux2(13:15) =-dJ(13:15,nls5:nus5)*dq(nls5:nus5)-dJ(13:15,nls9p3:nus9)*dq(nls9p3:nus9);
        psibarAux2(16:18) =-dJ(16:18,nls6:nus6)*dq(nls6:nus6)-dJ(16:18,nls9p3:nus9)*dq(nls9p3:nus9);
        psibarAux2(19:21) =-dJ(19:21,nls7:nus7)*dq(nls7:nus7)-dJ(19:21,nls9p3:nus9)*dq(nls9p3:nus9);
        psibarAux2(22:24) =-dJ(22:24,nls8:nus8)*dq(nls8:nus8)-dJ(22:24,nls9p3:nus9)*dq(nls9p3:nus9);
        psibar = psibarAux1+psibarAux2;
        %
        auxp1 = psibar(4:6)-auxa1*psibar(1:3);
        auxp2 = psibar(7:9)-auxa2*psibar(1:3)-auxb1*auxp1;
        auxp3 = psibar(10:12) -auxa3*psibar(1:3) - auxb2*auxp1 - auxc1*auxp2;
        auxp4 = psibar(13:15) - auxa4*psibar(1:3) - auxb3*auxp1 - auxc2*auxp2 - auxd1*auxp3;
        auxp5 = psibar(16:18) - auxa5*psibar(1:3) - auxb4*auxp1 - auxc3*auxp2 - auxd2*auxp3 -auxe1*auxp4;
        auxp6 = psibar(19:21) - auxa6*psibar(1:3) - auxb5*auxp1 - auxc4*auxp2 - auxd3*auxp3 -auxe2*auxp4 -auxf1*auxp5;
        auxp7 = psibar(22:24) - auxa7*psibar(1:3) - auxb6*auxp1 - auxc5*auxp2 - auxd4*auxp3 -auxe3*auxp4 -auxf2*auxp5 -auxg1*auxp6;
    end
end

%% Optimised the code below by performing block Gaussian elimination.
% lam = ibar\psibar;

% Back substitution

lam(22:24) = auxh*auxp7;
lam(19:21) = auxg*(auxp6-auxh06*lam(22:24));
lam(16:18) = auxf*(auxp5-auxg05*lam(19:21)-auxh05*lam(22:24));
lam(13:15) = auxe*(auxp4-auxf04*lam(16:18)-auxg04*lam(19:21)-auxh04*lam(22:24));
lam(10:12) = auxd*(auxp3-auxe03*lam(13:15)-auxf03*lam(16:18)-auxg03*lam(19:21)-auxh03*lam(22:24));
lam(7:9)   = auxc*(auxp2-auxd02*lam(10:12)-auxe02*lam(13:15)-auxf02*lam(16:18)-auxg02*lam(19:21)-auxh02*lam(22:24));
lam(4:6)   = auxb*(auxp1-auxc0*lam(7:9)-auxd01*lam(10:12)-auxe01*lam(13:15)-auxf01*lam(16:18)-auxg01*lam(19:21)-auxh01*lam(22:24));
lam(1:3)   = auxa*(psibar(1:3)-ibar(1:3,4:6)*lam(4:6)-ibar(1:3,7:9)*lam(7:9)-ibar(1:3,10:12)*lam(10:12)-ibar(1:3,13:15)*lam(13:15)-ibar(1:3,16:18)*lam(16:18)-ibar(1:3,19:21)*lam(19:21)-ibar(1:3,22:24)*lam(22:24));

%% Joint accelerations

% Optimised the code below using the sparsity of ihatTr
qdd(nls1:nus1) = phihat(nls1:nus1)+ihatTr(nls1:nus1,1:3)*lam(1:3);
qdd(nls2:nus2) = phihat(nls2:nus2)+ihatTr(nls2:nus2,4:6)*lam(4:6);
qdd(nls3:nus3) = phihat(nls3:nus3)+ihatTr(nls3:nus3,7:9)*lam(7:9);
qdd(nls4:nus4) = phihat(nls4:nus4)+ihatTr(nls4:nus4,10:12)*lam(10:12);
qdd(nls5:nus5) = phihat(nls5:nus5)+ihatTr(nls5:nus5,13:15)*lam(13:15);
qdd(nls6:nus6) = phihat(nls6:nus6)+ihatTr(nls6:nus6,16:18)*lam(16:18);
qdd(nls7:nus7) = phihat(nls7:nus7)+ihatTr(nls7:nus7,19:21)*lam(19:21);
qdd(nls8:nus8) = phihat(nls8:nus8)+ihatTr(nls8:nus8,22:24)*lam(22:24);
qdd(nls9:nus9) = phihat(nls9:nus9)+ihatTr(nls9:nus9,:)*lam;

end
