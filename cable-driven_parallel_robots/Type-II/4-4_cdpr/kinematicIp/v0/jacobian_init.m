% RSSLM-CDPR-Type-II-KinematicIp jacobian_init module. The entries of constraint Jacobian matrix are computed here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to constraints.m and constraints_init.m

function []=jacobian_init()

% System: 4-4 CDPR with movements of cables' exit points

% Global variables -- required
global b th n nc q dq tt tb so Qf p nls1 nls2 nls3 nls4 nls5 nus1 nus2 nus3 nus4 nus5;

% Global variables -- defined/modified
global J dJ qi dqi;

% Current position and velocity
qi=q;
dqi=dq;

% Initialisation
J = zeros(nc,n);
dJ = zeros(nc,n);

% Increment
h=0.5e-6;
invh=1e6;
twoh=h+h; 

[thc,bc]=constraints_init(th, b);

%% Computing elements of J and dJ

% sub-system 1
for ii = nls1:nus1
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus1, qi, dqi, tt, tb, so, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus1, qi, dqi, tt, tb, so, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(1:12,ii)=invh*(eta1-eta2);
    dJ(1:12,ii)=invh*(deta1-deta2);
end
J(13:15,nls1+1) = [1; 0; 0];
J(13:15,nls1+2) = [0; 1; 0];
J(13:15,nls1)=[0;0;1];

% sub-system 2
for ii = nls2:nus2
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus2, qi, dqi, tt, tb, so, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus2, qi, dqi, tt, tb, so, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(1:12,ii)=invh*(eta1-eta2);
    dJ(1:12,ii)=invh*(deta1-deta2);
end
J(16:18,nls2+1) = [1; 0; 0];
J(16:18,nls2+2) = [0; 1; 0];
J(16:18,nls2) = [0; 0; 1];

% sub-system 3
for ii = nls3:nus3
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus3, qi, dqi, tt, tb, so, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus3, qi, dqi, tt, tb, so, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(1:12,ii)=invh*(eta1-eta2);
    dJ(1:12,ii)=invh*(deta1-deta2);
end
J(19:21,nls3+1) = [1; 0; 0];
J(19:21,nls3+2) = [0; 1; 0];
J(19:21,nls3) = [0; 0; 1];

% sub-system 4
for ii = nls4:nus4
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, nus4, qi, dqi, tt, tb, so, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, nus4, qi, dqi, tt, tb, so, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(1:12,ii)=invh*(eta1-eta2);
    dJ(1:12,ii)=invh*(deta1-deta2);
end
J(22:24,nls4+1) = [1; 0; 0];
J(22:24,nls4+2) = [0; 1; 0];
J(22:24,nls4) = [0; 0; 1];

% sub-system 5
J(1:12,nls5+1) = [-1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0];
J(1:12,nls5+2) = [0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1; 0];
J(1:12,nls5) = [0; 0; -1; 0; 0; -1; 0; 0; -1; 0; 0; -1];
for ii = (nls5+3):nus5
    qi(ii)=qi(ii)+h;
    [eta1,deta1]=constraints(ii, -1, qi, dqi, tt, tb, so, Qf, p, bc, thc);
    qi(ii)=qi(ii)-twoh;
    [eta2,deta2]=constraints(ii, -1, qi, dqi, tt, tb, so, Qf, p, bc, thc);
    qi(ii)=qi(ii)+h;
    J(1:12,ii)=invh*(eta1-eta2);
    dJ(1:12,ii)=invh*(deta1-deta2);
end

end
