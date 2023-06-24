% RSSLM-CDPR-Type-II-KinematicIp data module. This module contains the information of geometric, elastic, and inertia parameters of the model of the mechanical system under study.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

function []=data() 

% Global variables -- defined
global nen nei ni n nc ali;

%System: 4-4 CDPR with each cable modelled by multiple modified rigid finite elements (MRFEs) with movements of the exit points of cables.

% Number of MRFEs associated with each sub-system (cable) of the 4-4 CDPR. 
nen = 5;
nei = nen*ones(4,1);

% Number of links for each of the sub-systems: 4 cables, 1 moving platform
ni=[3.*nei+3;6];

% Total number of links
n=sum(ni);

% Number of constraints
nc=24;

% Cable lengths (Half of the actual lengths)
li=375e-3*ones(4,1);

% Link lengths
ali=li./nei;

end
