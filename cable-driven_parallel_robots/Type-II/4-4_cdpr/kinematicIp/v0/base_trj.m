% RSSLM-CDPR-Type-II-KinematicIp base_trj module. This module contains the definitions of the input trajectories of the cables.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to ddq_tree_eff

% System: 4-4 CDPR with movements of the cables' exit points

function [bitrj] = base_trj(t)

% Global variable -- required
global tf;

ti1 = 0;
tf1 = 2;

ti2 = tf1;
tf2 = 3;

ti3 = tf2;
tf3 = 4;

ti4 = tf3;
tf4 = tf;

% Second derivative of the trajectories with respect to time

if t < tf1 
	utrj = 6*(ti1+tf1-2*t)/(tf1-ti1)^3;
	btrj = utrj*[0;0;1];
elseif t < tf2
	utrj = 6*(ti2+tf2-2*t)/(tf2-ti2)^3;
	btrj = utrj*[1/(2*sqrt(2)); 1/(2*sqrt(2)); sqrt(3)/2];
elseif t < tf3
	utrj = 6*(ti3+tf3-2*t)/(tf3-ti3)^3;
	btrj = utrj*[1/(2*sqrt(2)); 1/(2*sqrt(2)); -sqrt(3)/2];
else
	utrj = 6*(ti4+tf4-2*t)/(tf4-ti4)^3;
	btrj = utrj*[0;0;-1];
end

bitrj = -[btrj; btrj; btrj; btrj];

end