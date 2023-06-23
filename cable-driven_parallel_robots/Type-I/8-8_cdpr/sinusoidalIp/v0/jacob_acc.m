% RSSLM-CDPR-Type-I jacob_acc module. This module computes the state Jacobian matrix, i.e., the changes in joint accelerations w.r.t changes in the joint positions and velocities.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to invdyn_tree_eff, torque, acceleration

% System: 8-8 CDPR with (sinusoidal or linear) cable feed

function dfdy = jacob_acc(t,y)

% Global variables -- required
global n b th alp a r dx dy dz m Icxx Icyy Iczz Icxy Icyz Iczx;
global ddxdt ddydt ddzdt dmdt dIcxxdt dIcyydt dIczzdt dIcxydt dIcyzdt dIczxdt;

q=y(1:n);
dq=y(n+1:2*n);

% Increment/Pertubation
h    =0.5e-6;
invh =1e6;
twoh =h+h;

% Initialisation
Ja1   = [zeros(n,n), eye(n,n)];
Ja2q  = zeros(n,n);
Ja2dq = zeros(n,n);

% The tolerance (1e-5) at which the constraint Jacobian matrix and it's 
% derivative would be updated is greater than the increment (0.5 e-6)
% chosen here. Hence, J and dJ are assumed to be constant in computing dfdy.

% Optimal computation of the spring-damper forces/torques
tau_d = torque(q,dq);

% Finding C dq+tug using inverse dynamics algorithm
[tu] = invdyn_tree_eff(q, dq, zeros(n,1), b, th);

% Finding dMdt dq 
[tm] = dmdtqdot(q,dq,b,th);

% Optimal computation of the joint accelerations
[phihat, ihatTr, ibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxa6, auxa7, auxb, auxb1, auxb2, auxb3, auxb4, auxb5, auxb6, auxc0, auxc, auxc1, auxc2, auxc3, auxc4, auxc5, auxd01, auxd02, auxd, auxd1, auxd2, auxd3, auxd4, auxe01, auxe02, auxe03, auxe, auxe1, auxe2, auxe3, auxf01, auxf02, auxf03, auxf04, auxf, auxf1, auxf2, auxg01, auxg02, auxg03, auxg04, auxg05, auxg, auxg1, auxh01, auxh02, auxh03, auxh04, auxh05, auxh06, auxh, auxp1, auxp2, auxp3, auxp4, auxp5, auxp6, auxp7] = acceleration_init(q, dq, tu, tm, tau_d);
psibar = psibarAux1+psibarAux2;

for ii = 1:n 
    q(ii) = q(ii)+h;
        
    % Finding C dq+tug using inverse dynamics algorithm
    [tu1] = invdyn_tree_eff_sub_system(ii, tu, q, dq, b, th, alp, a, r, dx, dy, dz, m, Icxx, Icyy, Iczz, Icxy, Icyz, Iczx);
    % Finding dMdt dq
    [tm1] = dmdtqdot_sub_system(ii, tm, q, dq, b, th, alp, a, r, dx, dy, dz, m, ddxdt, ddydt, ddzdt, dmdt, dIcxxdt, dIcyydt, dIczzdt, dIcxydt, dIcyzdt, dIczxdt);
    % Finding spring-damper forces
    tau_d1 = tau_d;
    tau_d1(ii) = sprDamp_joint(ii,q(ii),dq(ii));
    % Variation of the joint acceleration w.r.t joint positions and velocities
    [ddq1] = acceleration_var(ii, 0, q, dq, tu1, tm1, tau_d1, phihat, ihatTr, ibar, psibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxa6, auxa7, auxb, auxb1, auxb2, auxb3, auxb4, auxb5, auxb6, auxc0, auxc, auxc1, auxc2, auxc3, auxc4, auxc5, auxd01, auxd02, auxd, auxd1, auxd2, auxd3, auxd4, auxe01, auxe02, auxe03, auxe, auxe1, auxe2, auxe3, auxf01, auxf02, auxf03, auxf04, auxf, auxf1, auxf2, auxg01, auxg02, auxg03, auxg04, auxg05, auxg, auxg1, auxh01, auxh02, auxh03, auxh04, auxh05, auxh06, auxh, auxp1, auxp2, auxp3, auxp4, auxp5, auxp6, auxp7);
    
    q(ii) = q(ii)-twoh;
    % Finding C dq+tug using inverse dynamics algorithm
    [tu2] = invdyn_tree_eff_sub_system(ii, tu, q, dq, b, th, alp, a, r, dx, dy, dz, m, Icxx, Icyy, Iczz, Icxy, Icyz, Iczx);
    % Finding dMdt dq
    [tm2] = dmdtqdot_sub_system(ii, tm, q, dq, b, th, alp, a, r, dx, dy, dz, m, ddxdt, ddydt, ddzdt, dmdt, dIcxxdt, dIcyydt, dIczzdt, dIcxydt, dIcyzdt, dIczxdt);
    % Finding input joint torque
    tau_d2 = tau_d;
    tau_d2(ii) = sprDamp_joint(ii,q(ii),dq(ii));
    
    % Variation of the joint acceleration w.r.t joint positions and velocities
    [ddq2] = acceleration_var(ii, 0, q, dq, tu2, tm2, tau_d2, phihat, ihatTr, ibar, psibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxa6, auxa7, auxb, auxb1, auxb2, auxb3, auxb4, auxb5, auxb6, auxc0, auxc, auxc1, auxc2, auxc3, auxc4, auxc5, auxd01, auxd02, auxd, auxd1, auxd2, auxd3, auxd4, auxe01, auxe02, auxe03, auxe, auxe1, auxe2, auxe3, auxf01, auxf02, auxf03, auxf04, auxf, auxf1, auxf2, auxg01, auxg02, auxg03, auxg04, auxg05, auxg, auxg1, auxh01, auxh02, auxh03, auxh04, auxh05, auxh06, auxh, auxp1, auxp2, auxp3, auxp4, auxp5, auxp6, auxp7);
    
    Ja2q(:,ii) = invh*(ddq1-ddq2);
    q(ii) = q(ii)+h;
    
    dq(ii) = dq(ii)+h;

    % Finding C dq+tug using inverse dynamics algorithm
    [tu3] = invdyn_tree_eff_sub_system(ii, tu, q, dq, b, th, alp, a, r, dx, dy, dz, m, Icxx, Icyy, Iczz, Icxy, Icyz, Iczx);
    % Finding dMdt dq
    [tm3] = dmdtqdot_sub_system(ii, tm, q, dq, b, th, alp, a, r, dx, dy, dz, m, ddxdt, ddydt, ddzdt, dmdt, dIcxxdt, dIcyydt, dIczzdt, dIcxydt, dIcyzdt, dIczxdt);
    % Finding input joint torque
    tau_d3 = tau_d;
    tau_d3(ii) = sprDamp_joint(ii,q(ii),dq(ii));
    
    % Variation of the joint acceleration w.r.t joint positions and velocities
    [ddq3] = acceleration_var(ii, 1, q, dq, tu3, tm3, tau_d3, phihat, ihatTr, ibar, psibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxa6, auxa7, auxb, auxb1, auxb2, auxb3, auxb4, auxb5, auxb6, auxc0, auxc, auxc1, auxc2, auxc3, auxc4, auxc5, auxd01, auxd02, auxd, auxd1, auxd2, auxd3, auxd4, auxe01, auxe02, auxe03, auxe, auxe1, auxe2, auxe3, auxf01, auxf02, auxf03, auxf04, auxf, auxf1, auxf2, auxg01, auxg02, auxg03, auxg04, auxg05, auxg, auxg1, auxh01, auxh02, auxh03, auxh04, auxh05, auxh06, auxh, auxp1, auxp2, auxp3, auxp4, auxp5, auxp6, auxp7);
    
    dq(ii) = dq(ii)-twoh;
    % Finding C dq+tug using inverse dynamics algorithm
    [tu4] = invdyn_tree_eff_sub_system(ii, tu, q, dq, b, th, alp, a, r, dx, dy, dz, m, Icxx, Icyy, Iczz, Icxy, Icyz, Iczx);
     % Finding dMdt dq
    [tm4] = dmdtqdot_sub_system(ii, tm, q, dq, b, th, alp, a, r, dx, dy, dz, m, ddxdt, ddydt, ddzdt, dmdt, dIcxxdt, dIcyydt, dIczzdt, dIcxydt, dIcyzdt, dIczxdt);
    % Finding input joint torque
    tau_d4 = tau_d;
    tau_d4(ii) = sprDamp_joint(ii,q(ii),dq(ii));
    
    % Variation of the joint acceleration w.r.t joint positions and velocities
    [ddq4] = acceleration_var(ii, 1, q, dq, tu4, tm4, tau_d4, phihat, ihatTr, ibar, psibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxa6, auxa7, auxb, auxb1, auxb2, auxb3, auxb4, auxb5, auxb6, auxc0, auxc, auxc1, auxc2, auxc3, auxc4, auxc5, auxd01, auxd02, auxd, auxd1, auxd2, auxd3, auxd4, auxe01, auxe02, auxe03, auxe, auxe1, auxe2, auxe3, auxf01, auxf02, auxf03, auxf04, auxf, auxf1, auxf2, auxg01, auxg02, auxg03, auxg04, auxg05, auxg, auxg1, auxh01, auxh02, auxh03, auxh04, auxh05, auxh06, auxh, auxp1, auxp2, auxp3, auxp4, auxp5, auxp6, auxp7);
    
    Ja2dq(:,ii) = invh*(ddq3-ddq4);
    dq(ii) = dq(ii)+h;
end

dfdy = [Ja1; [Ja2q, Ja2dq]];

end
