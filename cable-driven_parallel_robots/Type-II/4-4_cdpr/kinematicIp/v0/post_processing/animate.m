% RSSLM-CDPR-Type-II-KinematicIp animate module. This module depicts the motion of the system under study.

% Contributors: Dr. Teja Krishna Mamidi, Prof. Sandipan Bandyopadhyay @IIT Madras, 
% Acknowledgments: Dr. Suril V. Shah and Prof. S. K. Saha @IIT Delhi

% Function calls to data, inputs, and for_kine.

% System: 4-4 CDPR with movements of cables' exit points

function [] = animate()
disp('------------------------------------------------------------------');
disp('Animating the simulation data');

% Dependency
data;
inputs;

% global variables -- required
global n ni alp a b th bt r dx dy dz;
global nls1 nls2 nls3 nls4 nus1 nus2 nus3 nus4 nus5;

% Uncomment the below lines if the data is to be read from posvelacc.dat.
% load posvelacc.dat posvelacc;
% 
% zeroIndex = [];
% lenposvelacc = size(posvelacc,1);
% 
% for i=1:lenposvelacc
%     if posvelacc(i,2:(3*n+1))==0
%         zeroIndex=[zeroIndex; i];
%         disp(i);
%     end
% end
% 
% count=0;
% for j=1:length(zeroIndex)
%     posvelacc(zeroIndex(j)+count,:)=[];
%     count=count-1;
% end
% 
% Y=posvelacc(:,2:(3*n+1));
% T=posvelacc(:,1);

% Using the recorded data from the simulation -- comment the below five lines if posvelacc.dat was to be used.
load statevar.dat statevar;
load timevar.dat timevar;

Y=statevar;
T=timevar;

xmin=-1;
xmax=3;
zmin=-1;
zmax=3;
ymin=-1;
ymax=3;

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

% Uncomment to record the evolution of the configurations of cables
% fidcs1 = fopen('cable1datan5.dat','w');
% fidcs2 = fopen('cable2datan5.dat','w');
% fidcs3 = fopen('cable3datan5.dat','w');
% fidcs4 = fopen('cable4datan5.dat','w');
% fidmp = fopen('mpdatan5.dat','w');
% fidt = fopen('timedatan5.dat','w');

figure('Name','Animation Window','NumberTitle','off');
for i=1:length(T)
%for i=1
    q=Y(i,1:n);
    dq=Y(i,n+1:2*n)';
    [so, ~, ~, ~, st]=for_kine(q, dq, nus4, alp, a, b, th, bt, r, dx, dy, dz);
    
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
    
    %sub-system 5
    c73 = cos(q(nus5)); s73=sin(q(nus5));
    c72 = cos(q(nus5-1)-pi/2); s72=sin(q(nus5-1)-pi/2);
    c71 = cos(q(nus5-2)); s71=sin(q(nus5-2)); 
    rotmp = [c72*c73,              -c72*s73,                   s72
             c73*s71*s72 + c71*s73,  c71*c73 - s71*s72*s73,    -s71*c72
            -(c71*c73*s72)+s71*s73, c73*s71+c71*s72*s73,       c71*c72];
    mpcom = [q(nus5-4);q(nus5-3);q(nus5-5)];
    posa1 = mpcom + rotmp*[-3/10; -2/5; 1/10];
    posa2 = mpcom + rotmp*[3/10; -2/5; 1/10];
    posa3 = mpcom + rotmp*[-3/10; 2/5; 1/10]; 
    posa4 = mpcom + rotmp*[3/10; 2/5; 1/10];
    posa5 = mpcom + rotmp*[-3/10; -2/5; -1/10];
    posa6 = mpcom + rotmp*[3/10; -2/5; -1/10];
    posa7 = mpcom + rotmp*[-3/10; 2/5; -1/10]; 
    posa8 = mpcom + rotmp*[3/10; 2/5; -1/10];
    BS51X=[posa1(1);posa2(1);posa4(1);posa3(1);posa1(1);posa5(1);posa6(1);posa8(1);posa7(1);posa5(1)];
    BS51Y=[posa1(2);posa2(2);posa4(2);posa3(2);posa1(2);posa5(2);posa6(2);posa8(2);posa7(2);posa5(2)];
    BS51Z=[posa1(3);posa2(3);posa4(3);posa3(3);posa1(3);posa5(3);posa6(3);posa8(3);posa7(3);posa5(3)];
    BS52X=[posa2(1);posa6(1);posa8(1);posa4(1);posa3(1);posa7(1)];
    BS52Y=[posa2(2);posa6(2);posa8(2);posa4(2);posa3(2);posa7(2)];
    BS52Z=[posa2(3);posa6(3);posa8(3);posa4(3);posa3(3);posa7(3)];
    %
%     Extract cable data -- uncomment the below lines if the evolution of the shapes of cables are required.
%     for ii=1:1:length(BS1X)
%         fprintf(fidcs1,'%.16e ',BS1X(ii),BS1Y(ii),BS1Z(ii));
%         fprintf(fidcs1,'\n');
%     end
% 
%     for ii=1:1:length(BS2X)
%         fprintf(fidcs2,'%.16e ',BS2X(ii),BS2Y(ii),BS2Z(ii));
%         fprintf(fidcs2,'\n');
%     end
%  
%     for ii=1:1:length(BS3X)
%         fprintf(fidcs3,'%.16e ',BS3X(ii),BS3Y(ii),BS3Z(ii));
%         fprintf(fidcs3,'\n');
%     end
% 
%     for ii=1:1:length(BS4X)
%         fprintf(fidcs4,'%.16e ',BS4X(ii),BS4Y(ii),BS4Z(ii));
%         fprintf(fidcs4,'\n');
%     end
%     %
%     fprintf(fidmp,'%.16e ',posa1);
%     fprintf(fidmp,'%.16e ',posa2);
%     fprintf(fidmp,'%.16e ',posa3);
%     fprintf(fidmp,'%.16e ',posa4);
%     fprintf(fidmp,'%.16e ',posa5);
%     fprintf(fidmp,'%.16e ',posa6);
%     fprintf(fidmp,'%.16e ',posa7);
%     fprintf(fidmp,'%.16e ',posa8);
%     fprintf(fidmp,'\n');
%     %
%     fprintf(fidt,'%.16e',T(i));
%     fprintf(fidt,'\n');
    %
    t=round(T(i),2);
    t=num2str(t);
    plot3(BS1X(6:end),BS1Y(6:end),BS1Z(6:end),'r', BS2X(6:end),BS2Y(6:end),BS2Z(6:end),'b', BS3X(6:end),BS3Y(6:end),BS3Z(6:end),'g', BS4X(6:end),BS4Y(6:end),BS4Z(6:end),'m', BS51X,BS51Y,BS51Z,'k', BS52X,BS52Y,BS52Z,'k','lineWidth', 1);
    axis([ xmin xmax  ymin ymax zmin zmax]);
    set (gca,'fontsize',15,'fontweight','normal','fontname','times new romans','linewidth',0.5,'Box', 'off','TickDir','out' );
    xlabel('X (m)','fontweight','normal','fontsize',15);
    ylabel('Y (m)','fontweight','normal','fontsize',15);
    zlabel('Z (m)','fontweight','normal','fontsize',15);
    title(['Current time t=',t,' s'],'fontweight','normal','fontsize',15);
    grid on;
    view(90,0);
    %view(2);
    daspect([1, 1, 1]);
    %set(gcf, 'Position', [0 0 800 450]);
    drawnow;
    %saveas(gcf,"63cdprfast_initial.png");
end
fclose all;
end
