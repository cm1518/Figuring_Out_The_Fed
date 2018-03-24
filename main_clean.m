%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File to reproduce main figures of "Figuring out the Fed" by Christian
%Matthes
%output can be modified to generate other figures in the paper
%contact: christian.matthes@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

%load data

load data_final


%load posterior mode estimate

load posterior_mode

%setting some parameters needed for later
T=length(y);
z_initial=0;
g_initial=0;
p2_initial=zeros(2,1);


%call function to solve model 

[probspce, lcpce, ldpce, p2pce, x2pce, x1pce, E_dpce, E_cpce, M,C,Md,Cd,Fd]=filter_final(paramf,pipce,y,i,T, z_initial, g_initial,p2_initial);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reproducing the main figures from the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(dates,probspce(1,2:end),'LineWidth',2)
title('probability of commitment model - max posterior estimate','Fontsize',12)
axis tight

figure;
subplot(2,1,1)
plot(dates(2:end),400*pipce(2:end),dates(2:end),400*probspce(1,2:end-1).*E_cpce(1,1:end-1)+400*probspce(2,2:end-1).*E_dpce(1,1:end-1),'r','LineWidth',2);
legend('inflation rate','private sector expectations')
title('inflation expectations - Posterior Mode')
axis tight
subplot(2,1,2)
plot(dates(2:end),y(2:end),dates(2:end),probspce(1,2:end-1).*E_cpce(2,1:end-1)+probspce(2,2:end-1).*E_dpce(2,1:end-1),'r','LineWidth',2);
legend('output','private sector expectations')
title('output expectations - Posterior Mode')
axis tight


