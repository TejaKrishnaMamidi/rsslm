% RSSLM-CDPR-Kinetostatics data module. This module contains the information of the number of elements used for modelling a cable, their lengths, and the total number of links of the analysed multi-body system.

% Contibutors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

function []=data() 

% Global variables
global nen nei ni n ali;

% System: 6-3 CDPR with each cable modelled by multiple modified rigid finite elements (MRFEs)

% Number of MRFEs associated with each sub-system (cable) of the 6-3 CDPR. 
nen = 20;
nei = nen*[1;1;1;1;1;1];

% Number of links for each of the sub-systems
ni=[3*nei(1:6);6];

% Total number of links
n=sum(ni);

% Cable lengths (Half of the actual lengths)
li=[165.548038578429090884062913499; 163.725847309351392376191313313; 166.768033917156933117799048395; 164.065734366576274820655872905; 168.552581236984964950104483129; 166.526851825357787870363936235];

% Link lengths
ali=li./nei;

end
