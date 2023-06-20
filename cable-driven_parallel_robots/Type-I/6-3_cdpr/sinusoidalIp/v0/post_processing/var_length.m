% RSSLM-CDPR-Type-I var_length module. The variables dependent on the length of cables are computed here.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Prof. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% No function calls

%System: 6-3 CDPR with each cable modelled by multiple modified rigid finite elements and fed or retreived.

function [dy, dz]=var_length(ali) 

% Global variables -- required
global nls1 nls2 nls3 nls4 nls5 nls6 nus1 nus2 nus3 nus4 nus5 nus6;

% Global variables -- defined
global dy dz;

% Updating the parameters affected by changes in cables lengths


% sub-system 1

for ii=nls1:3:nus1
	iip1=ii+1;
	iip2=iip1+1;
	dy(iip1)=ali(1)/2;
	dz(iip2)=-ali(1)/2;
end

% sub-system 2


for jj=nls2:3:nus2
	jjp1=jj+1;
	jjp2=jjp1+1;
	dy(jjp1)=ali(2)/2;
	dz(jjp2)=-ali(2)/2;
end

% sub-system 3

for kk=nls3:3:nus3
	kkp1=kk+1;
	kkp2=kkp1+1;
	dy(kkp1)=ali(3)/2;
	dz(kkp2)=-ali(3)/2;
end

% sub-system 4

for ll=nls4:3:nus4
	llp1=ll+1;
	llp2=llp1+1;
	dy(llp1)=ali(4)/2;
	dz(llp2)=-ali(4)/2;
end

% sub-system 5

for mm=nls5:3:nus5
	mmp1=mm+1;
	mmp2=mmp1+1;
	dy(mmp1)=ali(5)/2;
	dz(mmp2)=-ali(5)/2;
end

% sub-system 6

for nn=nls6:3:nus6
	nnp1=nn+1;
	nnp2=nnp1+1;
	dy(nnp1)=ali(6)/2;
	dz(nnp2)=-ali(6)/2;
end

end


