% RSSLM-CDPR-Type-I animate module. This module depicts the motion of the system under study.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to data, var_length, and for_kine.

% System: 8-8 CDPR with cable feed

function [] = animate()
disp('------------------------------------------------------------------');
disp('Animating the simulation data');

% Dependency
data;

% global variables -- required
global n ni alp a b th bt r dx;
global nls1 nls2 nls3 nls4 nls5 nls6 nls7 nls8 nus1 nus2 nus3 nus4 nus5 nus6 nus7 nus8;

% Using the recorded data from the simulation
load posvelacc.dat posvelacc;
load elLen.dat elLen;

% Removes any zero entries from the recorded data.
zeroIndex = [];
lenposvelacc = size(posvelacc,1);

for i=1:lenposvelacc
    if posvelacc(i,2:(3*n+1))==0
        zeroIndex=[zeroIndex; i];
        disp(i);
    end
end

count=0;
for j=1:length(zeroIndex)
    posvelacc(zeroIndex(j)+count,:)=[];
    elLen(zeroIndex(j)+count,:)=[];
    count=count-1;
end

Y=posvelacc(:,2:(3*n+1));
T=posvelacc(:,1);

%load statevar.dat statevar
%load timevar.dat timevar

%Y=statevar;
%T=timevar;

% Plot limits
xmin=-10;
xmax=10;
zmin=0;
zmax=10;
ymin=-10;
ymax=10;

% Initialisation

BS1X = zeros(1,2*ni(1));
BS1Y = zeros(1,2*ni(1));
BS1Z = zeros(1,2*ni(1));

BS2X = zeros(1,2*ni(2));
BS2Y = zeros(1,2*ni(2));
BS2Z = zeros(1,2*ni(2));

BS3X = zeros(1,2*ni(3));
BS3Y = zeros(1,2*ni(3));
BS3Z = zeros(1,2*ni(3));

BS4X = zeros(1,2*ni(4));
BS4Y = zeros(1,2*ni(4));
BS4Z = zeros(1,2*ni(4));

BS5X = zeros(1,2*ni(5));
BS5Y = zeros(1,2*ni(5));
BS5Z = zeros(1,2*ni(5));

BS6X = zeros(1,2*ni(6));
BS6Y = zeros(1,2*ni(6));
BS6Z = zeros(1,2*ni(6));

BS7X = zeros(1,2*ni(7));
BS7Y = zeros(1,2*ni(7));
BS7Z = zeros(1,2*ni(7));

BS8X = zeros(1,2*ni(8));
BS8Y = zeros(1,2*ni(8));
BS8Z = zeros(1,2*ni(8));

% Uncomment to record the evolution of the configurations of cables
% fidcs1 = fopen('cable1datan20.dat','w');
% fidcs2 = fopen('cable2datan20.dat','w');
% fidcs3 = fopen('cable3datan20.dat','w');
% fidcs4 = fopen('cable4datan20.dat','w');
% fidcs5 = fopen('cable5datan20.dat','w');
% fidcs6 = fopen('cable6datan20.dat','w');
% fidcs7 = fopen('cable7datan20.dat','w');
% fidcs8 = fopen('cable8datan20.dat','w');
% fidOrient = fopen('orientdatan20.dat','w');
% fidPose = fopen('posedatan20.dat','w');
% fidt = fopen('timedatan20.dat','w');

figure('Name','Animation Window','NumberTitle','off');
for i=1:length(T)
%for i=1
    q=Y(i,1:n);
    dq=Y(i,n+1:2*n)';
    [dy, dz]= var_length(elLen(i,:));
    [so, ~, ~, ~, st]=for_kine(q, dq, n, alp, a, b, th, bt, r, dx, dy, dz);
    
    ind = 1;
    % sub-system 1
    for ii=nls1:3:nus1
	BS1X(ind)=so(1,ii);
	BS1X(ind+1)=st(1,ii);
	BS1X(ind+2)=so(1,ii+1);
	BS1X(ind+3)=st(1,ii+1);
	BS1X(ind+4)=st(1,ii+1);
	BS1X(ind+5)=so(1,ii+2);
	BS1Y(ind)=so(2,ii);
	BS1Y(ind+1)=st(2,ii);
	BS1Y(ind+2)=so(2,ii+1);
	BS1Y(ind+3)=st(2,ii+1);
	BS1Y(ind+4)=st(2,ii+1);
	BS1Y(ind+5)=so(2,ii+2);
	BS1Z(ind)=so(3,ii);
	BS1Z(ind+1)=st(3,ii);
	BS1Z(ind+2)=so(3,ii+1);
	BS1Z(ind+3)=st(3,ii+1);
	BS1Z(ind+4)=st(3,ii+1);
	BS1Z(ind+5)=so(3,ii+2);
	ind=ind+6;
    end
     
    ind = 1;		
    % sub-system 2
    for jj=nls2:3:nus2
	BS2X(ind)=so(1,jj);
	BS2X(ind+1)=st(1,jj);
	BS2X(ind+2)=so(1,jj+1);
	BS2X(ind+3)=st(1,jj+1);
	BS2X(ind+4)=st(1,jj+1);
	BS2X(ind+5)=so(1,jj+2);
	BS2Y(ind)=so(2,jj);
	BS2Y(ind+1)=st(2,jj);
	BS2Y(ind+2)=so(2,jj+1);
	BS2Y(ind+3)=st(2,jj+1);
	BS2Y(ind+4)=st(2,jj+1);
	BS2Y(ind+5)=so(2,jj+2);
	BS2Z(ind)=so(3,jj);
	BS2Z(ind+1)=st(3,jj);
	BS2Z(ind+2)=so(3,jj+1);
	BS2Z(ind+3)=st(3,jj+1);
	BS2Z(ind+4)=st(3,jj+1);
	BS2Z(ind+5)=so(3,jj+2);
	ind=ind+6;
    end	

    ind = 1;
    % sub-system 3
    for kk=nls3:3:nus3
	BS3X(ind)=so(1,kk);
	BS3X(ind+1)=st(1,kk);
	BS3X(ind+2)=so(1,kk+1);
	BS3X(ind+3)=st(1,kk+1);
	BS3X(ind+4)=st(1,kk+1);
	BS3X(ind+5)=so(1,kk+2);
	BS3Y(ind)=so(2,kk);
	BS3Y(ind+1)=st(2,kk);
	BS3Y(ind+2)=so(2,kk+1);
	BS3Y(ind+3)=st(2,kk+1);
	BS3Y(ind+4)=st(2,kk+1);
	BS3Y(ind+5)=so(2,kk+2);
	BS3Z(ind)=so(3,kk);
	BS3Z(ind+1)=st(3,kk);
	BS3Z(ind+2)=so(3,kk+1);
	BS3Z(ind+3)=st(3,kk+1);
	BS3Z(ind+4)=st(3,kk+1);
	BS3Z(ind+5)=so(3,kk+2);
	ind=ind+6;
    end	
    
    ind = 1;
    % sub-system 4
    for ll=nls4:3:nus4
	BS4X(ind)=so(1,ll);
	BS4X(ind+1)=st(1,ll);
	BS4X(ind+2)=so(1,ll+1);
	BS4X(ind+3)=st(1,ll+1);
	BS4X(ind+4)=st(1,ll+1);
	BS4X(ind+5)=so(1,ll+2);
	BS4Y(ind)=so(2,ll);
	BS4Y(ind+1)=st(2,ll);
	BS4Y(ind+2)=so(2,ll+1);
	BS4Y(ind+3)=st(2,ll+1);
	BS4Y(ind+4)=st(2,ll+1);
	BS4Y(ind+5)=so(2,ll+2);
	BS4Z(ind)=so(3,ll);
	BS4Z(ind+1)=st(3,ll);
	BS4Z(ind+2)=so(3,ll+1);
	BS4Z(ind+3)=st(3,ll+1);
	BS4Z(ind+4)=st(3,ll+1);
	BS4Z(ind+5)=so(3,ll+2);
	ind=ind+6;
    end

    ind = 1;	
    % sub-system 5
    for mm=nls5:3:nus5
	BS5X(ind)=so(1,mm);
	BS5X(ind+1)=st(1,mm);
	BS5X(ind+2)=so(1,mm+1);
	BS5X(ind+3)=st(1,mm+1);
	BS5X(ind+4)=st(1,mm+1);
	BS5X(ind+5)=so(1,mm+2);
	BS5Y(ind)=so(2,mm);
	BS5Y(ind+1)=st(2,mm);
	BS5Y(ind+2)=so(2,mm+1);
	BS5Y(ind+3)=st(2,mm+1);
	BS5Y(ind+4)=st(2,mm+1);
	BS5Y(ind+5)=so(2,mm+2);
	BS5Z(ind)=so(3,mm);
	BS5Z(ind+1)=st(3,mm);
	BS5Z(ind+2)=so(3,mm+1);
	BS5Z(ind+3)=st(3,mm+1);
	BS5Z(ind+4)=st(3,mm+1);
	BS5Z(ind+5)=so(3,mm+2);
	ind=ind+6;
    end
	
    ind = 1;	
    % sub-system 6
    for nn=nls6:3:nus6
	BS6X(ind)=so(1,nn);
	BS6X(ind+1)=st(1,nn);
	BS6X(ind+2)=so(1,nn+1);
	BS6X(ind+3)=st(1,nn+1);
	BS6X(ind+4)=st(1,nn+1);
	BS6X(ind+5)=so(1,nn+2);
	BS6Y(ind)=so(2,nn);
	BS6Y(ind+1)=st(2,nn);
	BS6Y(ind+2)=so(2,nn+1);
	BS6Y(ind+3)=st(2,nn+1);
	BS6Y(ind+4)=st(2,nn+1);
	BS6Y(ind+5)=so(2,nn+2);
	BS6Z(ind)=so(3,nn);
	BS6Z(ind+1)=st(3,nn);
	BS6Z(ind+2)=so(3,nn+1);
	BS6Z(ind+3)=st(3,nn+1);
	BS6Z(ind+4)=st(3,nn+1);
	BS6Z(ind+5)=so(3,nn+2);
    ind=ind+6;
    end

    ind = 1;	
    % sub-system 7
    for oo=nls7:3:nus7
	BS7X(ind)=so(1,oo);
	BS7X(ind+1)=st(1,oo);
	BS7X(ind+2)=so(1,oo+1);
	BS7X(ind+3)=st(1,oo+1);
	BS7X(ind+4)=st(1,oo+1);
	BS7X(ind+5)=so(1,oo+2);
	BS7Y(ind)=so(2,oo);
	BS7Y(ind+1)=st(2,oo);
	BS7Y(ind+2)=so(2,oo+1);
	BS7Y(ind+3)=st(2,oo+1);
	BS7Y(ind+4)=st(2,oo+1);
	BS7Y(ind+5)=so(2,oo+2);
	BS7Z(ind)=so(3,oo);
	BS7Z(ind+1)=st(3,oo);
	BS7Z(ind+2)=so(3,oo+1);
	BS7Z(ind+3)=st(3,oo+1);
	BS7Z(ind+4)=st(3,oo+1);
	BS7Z(ind+5)=so(3,oo+2);
    ind=ind+6;
    end

    ind = 1;	
    % sub-system 8
    for pp=nls8:3:nus8
	BS8X(ind)=so(1,pp);
	BS8X(ind+1)=st(1,pp);
	BS8X(ind+2)=so(1,pp+1);
	BS8X(ind+3)=st(1,pp+1);
	BS8X(ind+4)=st(1,pp+1);
	BS8X(ind+5)=so(1,pp+2);
	BS8Y(ind)=so(2,pp);
	BS8Y(ind+1)=st(2,pp);
	BS8Y(ind+2)=so(2,pp+1);
	BS8Y(ind+3)=st(2,pp+1);
	BS8Y(ind+4)=st(2,pp+1);
	BS8Y(ind+5)=so(2,pp+2);
	BS8Z(ind)=so(3,pp);
	BS8Z(ind+1)=st(3,pp);
	BS8Z(ind+2)=so(3,pp+1);
	BS8Z(ind+3)=st(3,pp+1);
	BS8Z(ind+4)=st(3,pp+1);
	BS8Z(ind+5)=so(3,pp+2);
    ind=ind+6;
    end
        % Positions of the trailing ends of the cables
%     fprintf('%.16e ',so(1:3,nus1));
%     fprintf('\n');
%     fprintf('%.16e ',so(1:3,nus2));
%     fprintf('\n');
%     fprintf('%.16e ',so(1:3,nus3));
%     fprintf('\n');
%     fprintf('%.16e ',so(1:3,nus4));
%     fprintf('\n');
%     fprintf('%.16e ',so(1:3,nus5));
%     fprintf('\n');
%     fprintf('%.16e ',so(1:3,nus6));
%     fprintf('\n');
%     fprintf('%.16e ',so(1:3,nus7));
%     fprintf('\n');
%     fprintf('%.16e ',so(1:3,nus8));
%     fprintf('\n');
    %
    % Extract cable data -- uncomment the below lines if the evolution of the shapes of cables are required.
%    for ii=1:1:length(BS1X)
%        fprintf(fidcs1,'%.16e ',BS1X(ii),BS1Y(ii),BS1Z(ii));
%        fprintf(fidcs1,'\n');
%    end
    %
%    for ii=1:1:length(BS2X)
%        fprintf(fidcs2,'%.16e ',BS2X(ii),BS2Y(ii),BS2Z(ii));
%        fprintf(fidcs2,'\n');
%    end
    %
%    for ii=1:1:length(BS3X)
%        fprintf(fidcs3,'%.16e ',BS3X(ii),BS3Y(ii),BS3Z(ii));
%        fprintf(fidcs3,'\n');
%    end
    %
%    for ii=1:1:length(BS4X)
%        fprintf(fidcs4,'%.16e ',BS4X(ii),BS4Y(ii),BS4Z(ii));
%        fprintf(fidcs4,'\n');
%    end
    %
%    for ii=1:1:length(BS5X)
%        fprintf(fidcs5,'%.16e ',BS5X(ii),BS5Y(ii),BS5Z(ii));
%        fprintf(fidcs5,'\n');
%    end
    %
%    for ii=1:1:length(BS6X)
%        fprintf(fidcs6,'%.16e ',BS6X(ii),BS6Y(ii),BS6Z(ii));
%        fprintf(fidcs6,'\n');
%    end
    %
%    for ii=1:1:length(BS7X)
%        fprintf(fidcs7,'%.16e ',BS7X(ii),BS7Y(ii),BS7Z(ii));
%        fprintf(fidcs7,'\n');
%    end
    %
%    for ii=1:1:length(BS8X)
%        fprintf(fidcs8,'%.16e ',BS8X(ii),BS8Y(ii),BS8Z(ii));
%        fprintf(fidcs8,'\n');
%    end
    %
%    fprintf(fidOrient,'%.16e ',q(n-2),q(n-1),q(n));
%    fprintf(fidOrient,'\n');
    %
%    fprintf(fidPose,'%.16e ',q(n-4),q(n-3),q(n-5));
%    fprintf(fidPose,'\n');
    %
%    fprintf(fidt,'%.16e',T(i));
%    fprintf(fidt,'\n');
    %
    t=round(T(i),2);
    t=num2str(t);
    plot3(BS1X,BS1Y,BS1Z,'r', BS2X,BS2Y,BS2Z,'b', BS3X,BS3Y,BS3Z,'r', BS4X,BS4Y,BS4Z,'b', BS5X,BS5Y,BS5Z,'r', BS6X,BS6Y,BS6Z,'b', BS7X,BS7Y,BS7Z,'r', BS8X,BS8Y,BS8Z,'b', 'lineWidth', 1);
    axis([ xmin xmax  ymin ymax zmin zmax]);
    set (gca,'fontsize',15,'fontweight','normal','fontname','times new romans','linewidth',0.5,'Box', 'off','TickDir','out' );
    xlabel('X (m)','fontweight','normal','fontsize',15);
    ylabel('Y (m)','fontweight','normal','fontsize',15);
    zlabel('Z (m)','fontweight','normal','fontsize',15);
    title(['Current time t=',t,' s'],'fontweight','normal','fontsize',15);
    grid on;
    view(3);
    daspect([1, 1, 1]);
    set(gcf, 'Position', [0 0 800 450]);
    drawnow;
    %saveas(gcf,"88cdprfast_initial.png");
end
fclose all;
end
