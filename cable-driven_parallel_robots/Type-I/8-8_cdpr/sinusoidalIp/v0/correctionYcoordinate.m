% RSSLM-CDPR-Type-I (or Kinetostatics) correctionYcoordinate module. The module corrects the y-coordinates of the exit points of the cables.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

function [biy]=correctionYcoordinate(j)

% No function calls

% System: 8-8 CDPR

% Global variables -- required
global nus1 nus2 nus3 nus4 nus5 nus6 nus7 nus8;

if j < nus1
    biy=-5.24398;
elseif j < nus2
    biy=-5.10296;
elseif j < nus3
    biy=5.23598;
elseif j < nus4
    biy=5.3476;
elseif j < nus5
    biy=5.37281;
elseif j < nus6
    biy=5.20584;
elseif j < nus7
    biy=-5.13255;
elseif j < nus8
    biy=-5.26946;
else
    biy=0;
end

end

