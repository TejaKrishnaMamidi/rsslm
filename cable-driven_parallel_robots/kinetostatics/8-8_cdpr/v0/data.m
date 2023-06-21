% RSSLM-CDPR-Kinetostatics data module. This module contains the information of the number of elements used for modelling a cable, their lengths, and the total number of links of the analysed multi-body system.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

function []=data() 

% Global variables
global nen nei ni n ali;

%System: 8-8 CDPR with each cable modelled by multiple modified rigid finite elements (MRFEs)

% Number of MRFEs associated with each sub-system (cable) of the 8-8 CDPR. 
nen = 20;
nei = nen*[1;1;1;1;1;1;1;1];

% Number of links for each of the sub-systems
ni=[3*nei(1:8);6];

% Total number of links
n=sum(ni);

% Cable lengths (Half of the actual lengths)
li=[5.24107; 4.91948; 5.08018; 4.48413; 5.155; 4.21081; 4.33162; 4.32778];

% Link lengths
ali=li./nei;

end
