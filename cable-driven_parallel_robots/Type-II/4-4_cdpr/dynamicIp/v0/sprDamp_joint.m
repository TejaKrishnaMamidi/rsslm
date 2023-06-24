% RSSLM-CDPR-Type-II-DynamicIp sprDamp_joint module. The module helps finding the changes in the generalised forces due to spring and damper elements with perturbations of joint positions and velocities.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

% System: 4-4 CDPR with cables attached to quadcopters

function [taud] = sprDamp_joint(ii, q, dq, taudConst, t)

% Global variables -- required
global ali ka ktr kgr ca ctr cgr fltr cdmp r nus1 nus2 nus3 nus4 nus5 nus6 nus7 nus8 nus9 nls1 nls2 nls3 nls4 nls5 nls6 nls7 nls8 nls9;

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
elseif ii<=nus5
    jj=0;
    % Moving platform -- damping characteristics
    taud(ii)=-cdmp*dq(ii)^2*sign(dq(ii));
    % Moving platform -- ground reactions
    if ii==nls5 && q(ii) <= fltr(nls5)
        taud(ii)=kgr*(fltr(nls5)-q(ii))-cgr*dq(ii);
    end
elseif ii<=nus6
    jj=0;
    %c63 = cos(q(nus6)); s63=sin(q(nus6));
    c62 = cos(q(nus6-1)-pi/2); s62=sin(q(nus6-1)-pi/2);
    c61 = cos(q(nus6-2)); s61=sin(q(nus6-2)); 
    %rotq1 = [c62*c63,               -c62*s63,               s62
    %         c63*s61*s62 + c61*s63,  c61*c63 - s61*s62*s63, -s61*c62
    %        -(c61*c63*s62)+s61*s63, c63*s61+c61*s62*s63,   c61*c62];
    % Quadcopter-1 -- thrust forces
    taud(nls6:nus6)=force_trj(1, [c61*c62; s62; -s61*c62], t);
    % Quadcopter-1 -- ground reactions
%    if ii==nls6 && q(ii) <= fltr(nls6)
%        taud(ii)=taud(ii)+kgr*(fltr(nls6)-q(ii))-cgr*dq(ii);
%    end    
elseif ii<=nus7
    jj=0;
    %c73 = cos(q(nus7)); s73=sin(q(nus7));
    c72 = cos(q(nus7-1)-pi/2); s72=sin(q(nus7-1)-pi/2);
    c71 = cos(q(nus7-2)); s71=sin(q(nus7-2)); 
    %rotq2 = [c72*c73,               -c72*s73,               s72
    %         c73*s71*s72 + c71*s73,  c71*c73 - s71*s72*s73, -s71*c72
    %        -(c71*c73*s72)+s71*s73, c73*s71+c71*s72*s73,   c71*c72];
    % Quadcopter-2 -- thrust forces
    taud(nls7:nus7)=force_trj(2, [c71*c72; s72; -s71*c72], t);
    % Quadcopter-2 -- ground reactions
%    if ii==nls7 && q(ii) <= fltr(nls7)
%        taud(ii)=taud(ii)+kgr*(fltr(nls7)-q(ii))-cgr*dq(ii);
%    end    
elseif ii<=nus8
    jj=0;
    %c83 = cos(q(nus8)); s83=sin(q(nus8));
    c82 = cos(q(nus8-1)-pi/2); s82=sin(q(nus8-1)-pi/2);
    c81 = cos(q(nus8-2)); s81=sin(q(nus8-2)); 
    %rotq3 = [c82*c83,               -c82*s83,               s82
    %         c83*s81*s82 + c81*s83,  c81*c83 - s81*s82*s83, -s81*c82
    %        -(c81*c83*s82)+s81*s83, c83*s81+c81*s82*s83,   c81*c82];
    % Quadcopter-3 -- thrust forces
    taud(nls8:nus8)=force_trj(3, [c81*c82; s82; -s81*c82], t);
    % Quadcopter-3 -- ground reactions
%    if ii==nls8 && q(ii) <= fltr(nls8)
%        taud(ii)=taud(ii)+kgr*(fltr(nls8)-q(ii))-cgr*dq(ii);
%    end    
elseif ii<=nus9
    jj=0;
    %c93 = cos(q(nus9)); s93=sin(q(nus9));
    c92 = cos(q(nus9-1)-pi/2); s92=sin(q(nus9-1)-pi/2);
    c91 = cos(q(nus9-2)); s91=sin(q(nus9-2)); 
    %rotq4 = [c92*c93,               -c92*s93,               s92
    %         c93*s91*s92 + c91*s93,  c91*c93 - s91*s92*s93, -s91*c92
    %        -(c91*c93*s92)+s91*s93, c93*s91+c91*s92*s93,   c91*c92];
    % Quadcopter-4 -- thrust forces
    taud(nls9:nus9)=force_trj(4, [c91*c92; s92; -s91*c92], t);
    % Quadcopter-4 -- ground reactions
%    if ii==nls9 && q(ii) <= fltr(nls9)
%        taud(ii)=taud(ii)+kgr*(fltr(nls9)-q(ii))-cgr*dq(ii);
%    end
else    
    jj=0;
end

% Associated changes in the generalised forces of spring-damper elements

if r(ii)==0 && jj~=0
    taud(ii)=ka(jj)*(2*ali(jj) - q(ii)) - ca(jj)*dq(ii);
elseif r(ii)~=0 && jj~=0
    taud(ii)=ktr(jj)*(fltr(ii) - q(ii)) - ctr(jj)*dq(ii);     
end

end
