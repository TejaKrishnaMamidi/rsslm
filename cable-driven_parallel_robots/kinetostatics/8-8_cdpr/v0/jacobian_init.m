% RSSLM-CDPR-Kinetostatics jacobian_init module. The entries of constraint Jacobian matrix are computed here.

% Contibutors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to constraints.m and constraints_init.m

function []=jacobian_init()

% System: 8-8 CDPR

% Global variables -- required
global b th n q dq tt tb so st vt Qf p;
global nls1 nls2 nls3 nls4 nls5 nls6 nls7 nls8 nls9 nus1 nus2 nus3 nus4 nus5 nus6 nus7 nus8 nus9;

% Global variables -- defined
global qi dqi J dJ;

% Current position and velocity
qi=q;
dqi=dq;

% Initialisation
J = zeros(24,n);
dJ = zeros(24,n);

% Increment
h=0.5e-6;
invh=1e6;
twoh=h+h; 

[thc,bc]=constraints_init(th, b);

%% Computing elements of J and dJ

% sub-system 1
for ii = nls1:nus1
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus1, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus1, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(:,ii)=invh*(eta1-eta2);
    dJ(:,ii)=invh*(deta1-deta2);
end
% sub-system 2
for ii = nls2:nus2
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus2, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus2, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(:,ii)=invh*(eta1-eta2);
    dJ(:,ii)=invh*(deta1-deta2);
end
% sub-system 3
for ii = nls3:nus3
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus3, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus3, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(:,ii)=invh*(eta1-eta2);
    dJ(:,ii)=invh*(deta1-deta2);
end
% sub-system 4
for ii = nls4:nus4
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus4, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus4, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(:,ii)=invh*(eta1-eta2);
    dJ(:,ii)=invh*(deta1-deta2);
end
% sub-system 5
for ii = nls5:nus5
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus5, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus5, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(:,ii)=invh*(eta1-eta2);
    dJ(:,ii)=invh*(deta1-deta2);
end
% sub-system 6
for ii = nls6:nus6
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus6, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus6, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(:,ii)=invh*(eta1-eta2);
    dJ(:,ii)=invh*(deta1-deta2);
end
%
% sub-system 7
for ii = nls7:nus7
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus7, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus7, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(:,ii)=invh*(eta1-eta2);
    dJ(:,ii)=invh*(deta1-deta2);
end
%
% sub-system 8
for ii = nls8:nus8
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus8, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus8, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(:,ii)=invh*(eta1-eta2);
    dJ(:,ii)=invh*(deta1-deta2);
end
%
% sub-system 9
J(:,nls9+1) = [-1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0];
J(:,nls9+2) = [0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0];
J(:,nls9) = [0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1];
for ii = (nls9+3):nus9
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, -1, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, -1, qi, dqi, tt, tb, so, st, vt, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(:,ii)=invh*(eta1-eta2);
    dJ(:,ii)=invh*(deta1-deta2);
end
%

end
