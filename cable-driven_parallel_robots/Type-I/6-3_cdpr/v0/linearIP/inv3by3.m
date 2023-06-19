% RSSLM-CDPR-Type-I inv3by3 module. A pre-computed symbolic expressions of the inverse of a non-singular matrix of dimension 3 by 3 is used for determining the matrix inverse.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

function [inva] = inv3by3(matA)

detA = 2*matA(1,2)*matA(2,3)*matA(1,3) - (matA(1,3)^2)*matA(2,2) - (matA(2,3)^2)*matA(1,1) - (matA(1,2)^2)*matA(3,3) + matA(1,1)*matA(2,2)*matA(3,3);

%if min(abs(svds(matA)))>1e-15
if abs(detA)>1e-15
    inva12 = matA(1,3)* matA(2,3)-matA(1,2)*matA(3,3);
    inva13 = matA(1,2)*matA(2,3)-matA(1,3)* matA(2,2);
    inva23 = matA(1,3)*matA(1,2)-matA(1,1)*matA(2,3);
    inva = (1/detA)*[matA(2,2)*matA(3,3)-matA(2,3)^2, inva12 , inva13; inva12, matA(1,1)*matA(3,3)-matA(1,3)^2, inva23; inva13, inva23, matA(1,1)*matA(2,2)-matA(1,2)^2];
else
    inva = pinv(matA);
end

end
