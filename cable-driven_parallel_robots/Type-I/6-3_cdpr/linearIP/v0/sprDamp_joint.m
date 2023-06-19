% RSSLM-CDPR-Type-I sprDamp_joint module. The module helps finding the changes in the generalised forces due to spring and damper elements with perturbations of joint positions and velocities.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

% System: 6-3 CDPR with cable feed

function [tauii] = sprDamp_joint(ii, q, dq)

% Global variables -- required
global ali ka ktr ca ctr fltr cdmp r dvi;
global nus1 nus2 nus3 nus4 nus5 nus6;

% Choice of the sub-system
if ii<= nus1
    jj = 1;
elseif ii<=nus2
    jj = 2;
elseif ii<=nus3
    jj = 3;
elseif ii<=nus4
    jj = 4;
elseif ii<=nus5
    jj = 5;
elseif ii<=nus6
    jj = 6;
else
    jj = 0;
end

% Associated change in the generalised forces
if jj==0
    tauii=-cdmp*dq^2*sign(dq);
elseif r(ii)==0
    tauii=ka(jj)*(2*ali(jj) - q) - ca(jj)*(dq-2*dvi(jj));
elseif r(ii)~=0
    tauii=ktr(jj)*(fltr(ii) - q) - ctr(jj)*dq;
end

end
