% RSSLM-CDPR-Type-II-DynamicIp data module. TThis module contains the information of the number of elements used for modelling a cable, their lengths, and the total number of links of the analysed multi-body system.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

function []=data() 

% Global variables -- defined
global nen nei ni n ali;

%System: 4-4 CDPR with each cable modelled by multiple modified rigid finite elements (MRFEs) and its cables attached to quadcopters

% Number of MRFEs associated with each sub-system (cable) of the 4-4 CDPR. 
nen = 5;
nei = nen*ones(4,1);

% Number of links for each of the sub-systems: 4 cables, 1 moving platform, 4 quadrotors
ni=[3.*nei+3;6.*ones(5,1)];

% Total number of links
n=sum(ni);

% Cable lengths (Half of the actual lengths)
li=375e-3*ones(4,1);

% Link lengths
ali=li./nei;

end
