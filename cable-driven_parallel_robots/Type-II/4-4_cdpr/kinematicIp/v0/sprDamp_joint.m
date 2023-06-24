% RSSLM-CDPR-Type-II-KinematicIp sprDamp_joint module. The module helps finding the changes in the generalised forces due to spring and damper elements with perturbations of joint positions and velocities.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

% System: 4-4 CDPR with movements of cables' exit points

function [taud] = sprDamp_joint(ii, q, dq, taudConst)

% Global variables -- required
global ali ka ktr kgr ca ctr cgr fltr cdmp r nus1 nus2 nus3 nus4 nls1 nls2 nls3 nls4 nls5;

% Initialisation
taud = taudConst;

% Choice of the sub-system and exit points
if ii<=nls1+2
    jj=0;
elseif ii<= nus1
    jj=1;
elseif ii<=nls2+2
    jj=0;  
elseif ii<=nus2 
    jj=2;
elseif ii<=nls3+2    
    jj=0;      
elseif ii<=nus3 
    jj=3;
elseif ii<=nls4+2    
    jj=0;
elseif ii<=nus4
    jj=4;
else
    jj=0;
    % Moving platform -- damping characteristics
    taud(ii)=-cdmp*dq(ii)^2*sign(dq(ii));
    % Moving platform -- ground reactions
    if ii==nls5 && q(ii) <= fltr(nls5)
        taud(ii)=kgr*(fltr(nls5)-q(ii))-cgr*dq(ii);
    end
end

% Associated changes in the generalised forces
if r(ii)==0 && jj~=0
    taud(ii)=ka(jj)*(2*ali(jj) - q(ii)) - ca(jj)*dq(ii);
elseif r(ii)~=0 && jj~=0
    taud(ii)=ktr(jj)*(fltr(ii) - q(ii)) - ctr(jj)*dq(ii);     
end

end
