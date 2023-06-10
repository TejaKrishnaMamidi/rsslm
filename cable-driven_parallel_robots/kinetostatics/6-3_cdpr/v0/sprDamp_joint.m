% RSSLM-CDPR-Kinetostatics sprDamp_joint module. The module helps finding the changes in the generalised forces due to spring and damper elements with perturbations of joint positions and velocities.

% Contibutors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

% System: 6-3 CDPR

function [tauii] = sprDamp_joint(ii, q, dq)

% Global variables -- required
global ali ka ktr ca ctr fltr cdmp r;
global nus1 nus2 nus3 nus4 nus5 nus6;

% Choice of the sub-system
if ii<= nus1
    fla = ali(1);
elseif ii<=nus2
    fla = ali(2);
elseif ii<=nus3
    fla = ali(3);
elseif ii<=nus4
    fla = ali(4);
elseif ii<=nus5
    fla = ali(5);
elseif ii<=nus6
    fla = ali(6);
else
    fla = 0;
end

% Associated change in the generalised forces
if fla==0
    tauii=-cdmp*dq^2*sign(dq);
elseif r(ii)==0
    tauii=ka*(2*fla - q) - ca*dq;
elseif r(ii)~=0
    tauii=ktr*(fltr(ii) - q) - ctr*dq;
end

end
