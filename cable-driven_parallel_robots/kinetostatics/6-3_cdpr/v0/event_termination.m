% RSSLM-CDPR-Kinetostatics event_termination module. This module terminates the simulation when the specified event occurs during the dynamic evolution of the configuration of CDPR.

% Contibutors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

% System: 6-3 CDPR

function [cumValue, isterminal, direction] = event_termination(t,y)

% Global variables -- required
global n ddq posIndex pos jtor cflam fid1 fid2 fid3 nse;

% Uncomment for computing the residue of the equation of motion
%global lam tau_d I J tu;

dq(1:n,1)= y(n+1:2*n);

cumValue=1;
value=ones(2,1);

%  Maximum of the absolute of joint velocities and accelerations
value(1) = max(abs(dq));
value(2) = max(abs(ddq));

% For display
disp('Max. velocity or acceleration:');
fprintf('%.16e\n', max(abs(value(1))));
fprintf('%.16e\n', max(abs(value(2))));
%fprintf('%.16e\n', max(abs(-I*ddq+tau_d-tu+J'*lam)));

%% To ensure a static equilibrium configuration is attained.

% Check for the joint velocities
if value(1)<1e-2
   value(1)=0; 
end
% Check for the joint accelerations
if value(2)<1e-2
    value(2)=0;
end
% Check for the joint velocities and accelerations to be less than a numerical value for five consecutive time instances.
if value(1)==0 && value(2)==0
    nse = nse -1;
    if nse == 0
        cumValue=0;
        if posIndex~=1
            for ii=1:posIndex
                fprintf(fid1,"%.16e ",pos(:,ii));
                fprintf(fid1,"\n");
                fprintf(fid2,"%.16e ",jtor(:,ii));
                fprintf(fid2,"\n");
                fprintf(fid3,"%.16e ",cflam(:,ii));
                fprintf(fid3,"\n");
            end
            posIndex=1;
        end
    end
else
    nse = 5;
end
%
fprintf('%.16e\n', nse);
%
isterminal=1;
direction=0;
%

end
