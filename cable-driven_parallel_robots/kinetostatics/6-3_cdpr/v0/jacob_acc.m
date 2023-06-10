% RSSLM-CDPR-Kinetostatics jacob_acc module. This module computes the state Jacobian matrix, i.e., the changes in joint accelerations w.r.t changes in the joint positions and velocities.

% Contibutors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to invdyn_tree_eff, torque, acceleration

% System: 6-3 CDPR

function dfdy = jacob_acc(t,y)

% Global variables -- required
global n q dq b th alp a r dx dy dz m Icxx Icyy Iczz Icxy Icyz Iczx;

% Increment/Pertubation
h    =0.5e-6;
invh =1e6;
twoh =h+h;

%% Initialisation
Ja1   = [zeros(n,n), eye(n,n)];
Ja2q  = zeros(n,n);
Ja2dq = zeros(n,n);

% The tolerance (1e-3) at which the constraint Jacobian matrix and it's 
% derivative would be updated is greater than the increment (0.5 e-6)
% chosen here. Hence, J and dJ are assumed to be constant while computing the state Jacobian matrix.

% Optimal computation of the spring-damper forces/torques
tau_d = torque(q,dq);

% Optimal computation of the joint accelerations
[tu] = invdyn_tree_eff(q, dq, zeros(n,1), b, th);
[phihat, ihatTr, ibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxb, auxb1, auxb2, auxb3, auxb4, auxc0, auxc, auxc1, auxc2, auxc3, auxd01, auxd02, auxd, auxd1, auxd2, auxe01, auxe02, auxe03, auxe, auxe1, auxf01, auxf02, auxf03, auxf04, auxp1, auxp2, auxp3, auxp4, auxp5] = acceleration_init(q, dq, tu, tau_d);
psibar = psibarAux1+psibarAux2;

for ii = 1:n 
    q(ii) = q(ii)+h;
        
    % Finding C dq+tug using inverse dynamics algorithm
    [tu1] = invdyn_tree_eff_sub_system(ii, tu, q, dq, b, th, alp, a, r, dx, dy, dz, m, Icxx, Icyy, Iczz, Icxy, Icyz, Iczx);
    % Finding spring-damper forces
    tau_d1 = tau_d;
    tau_d1(ii) = sprDamp_joint(ii,q(ii),dq(ii));
    % Variation of the joint accelerations w.r.t joint positions
    [ddq1] = acceleration_var(ii, 0, q, dq, tu1, tau_d1, phihat, ihatTr, ibar, psibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxb, auxb1, auxb2, auxb3, auxb4, auxc0, auxc, auxc1, auxc2, auxc3, auxd01, auxd02, auxd, auxd1, auxd2, auxe01, auxe02, auxe03, auxe, auxe1, auxf01, auxf02, auxf03, auxf04, auxp1, auxp2, auxp3, auxp4, auxp5);
    
    q(ii) = q(ii)-twoh;
    % Finding C dq+tug using inverse dynamics algorithm
    [tu2] = invdyn_tree_eff_sub_system(ii, tu, q, dq, b, th, alp, a, r, dx, dy, dz, m, Icxx, Icyy, Iczz, Icxy, Icyz, Iczx);
    % Finding input joint torque
    tau_d2 = tau_d;
    tau_d2(ii) = sprDamp_joint(ii,q(ii),dq(ii));
    
    % Variation of the joint accelerations w.r.t joint positions
    [ddq2] = acceleration_var(ii, 0, q, dq, tu2, tau_d2, phihat, ihatTr, ibar, psibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxb, auxb1, auxb2, auxb3, auxb4, auxc0, auxc, auxc1, auxc2, auxc3, auxd01, auxd02, auxd, auxd1, auxd2, auxe01, auxe02, auxe03, auxe, auxe1, auxf01, auxf02, auxf03, auxf04, auxp1, auxp2, auxp3, auxp4, auxp5);
    
    Ja2q(:,ii) = invh*(ddq1-ddq2);
    q(ii) = q(ii)+h;
    
    dq(ii) = dq(ii)+h;

    % Finding C dq+tug using inverse dynamics algorithm
    [tu3] = invdyn_tree_eff_sub_system(ii, tu, q, dq, b, th, alp, a, r, dx, dy, dz, m, Icxx, Icyy, Iczz, Icxy, Icyz, Iczx);
    % Finding input joint torque
    tau_d3 = tau_d;
    tau_d3(ii) = sprDamp_joint(ii,q(ii),dq(ii));
    
    % Variation of the joint accelerations w.r.t joint velocities
    [ddq3] = acceleration_var(ii, 1, q, dq, tu3, tau_d3, phihat, ihatTr, ibar, psibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxb, auxb1, auxb2, auxb3, auxb4, auxc0, auxc, auxc1, auxc2, auxc3, auxd01, auxd02, auxd, auxd1, auxd2, auxe01, auxe02, auxe03, auxe, auxe1, auxf01, auxf02, auxf03, auxf04, auxp1, auxp2, auxp3, auxp4, auxp5);
    
    dq(ii) = dq(ii)-twoh;
    % Finding C dq+tug using inverse dynamics algorithm
    [tu4] = invdyn_tree_eff_sub_system(ii, tu, q, dq, b, th, alp, a, r, dx, dy, dz, m, Icxx, Icyy, Iczz, Icxy, Icyz, Iczx);
    % Finding input joint torque
    tau_d4 = tau_d;
    tau_d4(ii) = sprDamp_joint(ii,q(ii),dq(ii));
    
    % Variation of the joint accelerations w.r.t joint velocities
    [ddq4] = acceleration_var(ii, 1, q, dq, tu4, tau_d4, phihat, ihatTr, ibar, psibar, psibarAux1, psibarAux2, auxa, auxa1, auxa2, auxa3, auxa4, auxa5, auxb, auxb1, auxb2, auxb3, auxb4, auxc0, auxc, auxc1, auxc2, auxc3, auxd01, auxd02, auxd, auxd1, auxd2, auxe01, auxe02, auxe03, auxe, auxe1, auxf01, auxf02, auxf03, auxf04, auxp1, auxp2, auxp3, auxp4, auxp5);
    
    Ja2dq(:,ii) = invh*(ddq3-ddq4);
    dq(ii) = dq(ii)+h;
end

dfdy = [Ja1; [Ja2q, Ja2dq]];

end