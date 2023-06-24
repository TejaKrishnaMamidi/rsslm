% RSSLM-CDPR-Type-II-DynamicIp force_trj module. This module contains the input variations of the thrust forces of the quadcopters.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to ddq_tree_eff

% System: 4-4 CDPR with cables attached to quadcopters

function [fip] = force_trj(qcn, rot, t)

%% Initialisation

f1=2.5; f2=2.5; f3=2.5; f4=2.5;
M1=f1*1056e-4; M2=f2*1056e-4; 
M3=-M1; M4=-M2;
a=169e-3;

%% Assignment

tsw1=5e-1;
tsw2=3;

mag=[rot.*(f1+f2+f3+f4); (f2-f4)*a; (f3-f1)*a; M1+M2+M3+M4];

% 
if t<=tsw1
    fip = sin((pi/2)*t/tsw1).*mag;
elseif qcn==1 && t>tsw2 % Failure of the quadcopter
    fip = zeros(6,1);
else
    fip = mag;
end

end
