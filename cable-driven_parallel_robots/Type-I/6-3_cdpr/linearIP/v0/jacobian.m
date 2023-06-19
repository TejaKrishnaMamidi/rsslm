% RSSLM-CDPR-Type-I (or Kinetostatics) jacobian module. This code decides when to update the constraint Jacobain matrix and it's derivative.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to jacobian_init.m

% System: 6-3 CDPR with feed

function []=jacobian()

% Global variables -- required
global qi dqi q dq;

dpos = max(abs(qi-q));
dvel = max(abs(dqi-dq));

%disp([dpos,dvel]);

% If changes in any of the joint positions or velocities are more than the below tolerance, then the constraint Jacobian matirx is updated.  
tol=1e-5;

if dpos >=tol || dvel >=tol
    jacobian_init()
end

end
