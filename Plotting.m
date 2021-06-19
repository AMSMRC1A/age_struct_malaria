% close all
clear all
clc
format long
global P
global colour_r1 colour_r2

load('Output/config.mat','P')
tfinal = P.t(end);
t = P.t;
a = P.a;
dt = P.dt;

%% human population size in time by stages
SH = FileRead('Output/','SH','Rave');
EH = FileRead('Output/','EH','Rave');
AH = FileRead('Output/','AH','Rave');
DH = FileRead('Output/','DH','Rave');
RH = FileRead('Output/','RH','Rave');
NH = SH+EH+AH+DH+RH;
figure_setups;
plot(t,SH,'b-'); hold on
plot(t,EH,'-','Color',colour_r1); 
plot(t,AH,'-','Color',colour_r2); 
plot(t,DH,'r-'); 
plot(t,RH,'g-'); 
plot(t,NH,'-'); 
legend('SH-age','EH-age','AH-age', 'DH-age','RH-age','$N_H$');
title('human population size by stages')
axis_years(gca,tfinal); % change to x-axis to years if needed
grid on
axis([0 tfinal 0 1.1])

%% Mosquito infection in time by stages
SM = FileRead('Output/','SM','Read');
EM = FileRead('Output/','EM','Read');
IM = FileRead('Output/','IM','Read');
NM = SM+EM+IM;
figure_setups;
plot(t/365,SM,'b-'); hold on;
plot(t/365,EM,'-','Color',colour_r1);
plot(t/365,IM,'r-');
plot(t/365,SM+EM+IM)
legend('SM','EM','IM','$N_M$');
title('mosquito population size by stages')
grid on
axis([0 tfinal/365 0 5])

%% final age distributions
SH = FileRead('Output/','SH','Final');
EH = FileRead('Output/','EH','Final');
AH = FileRead('Output/','AH','Final');
DH = FileRead('Output/','DH','Final');
RH = FileRead('Output/','RH','Final');
figure_setups;
plot(a/365,SH,'b-'); hold on;
plot(a/365,EH,'-','Color',colour_r1);
plot(a/365,AH,'-','Color',colour_r2);
plot(a/365,DH,'r-');
plot(a/365,RH,'g-');
legend('SH-age','EH-age','AH-age', 'DH-age','RH-age');
title('Age distributions by class (final time)')
grid on

%% Write output for convergence check
% Cs_t = FileRead('Output/','Cs','Rave');
% Cs_end = FileRead('Output/','Cs','Final');
% save(['Results/solution_',num2str(dt),'.mat'],'P','t','Cs_t','Cs_end')

%% Compare initial and final age distribution
SH = FileRead('Output/','SH','Final');
EH = FileRead('Output/','EH','Final');
AH = FileRead('Output/','AH','Final');
DH = FileRead('Output/','DH','Final');
RH = FileRead('Output/','RH','Final');
NH = SH+EH+AH+DH+RH;
figure_setups;
plot(a/365, NH,'g-'); hold on;
plot(a/365, P.n_tilde,'-','Color',colour_r2);
plot(a/365, P.n_tilde*exp(P.p_hat*tfinal) ,'-.');
legend('Final Age dist. (sim.)','Stable Age Dist. (initial)','Final Age dist. (theory)');
title('Population Age Distributions')

%% Immunity related figures
Cs_end = FileRead('Output/','Cs','Final');
figure_setups; 
subplot(2,2,1), plot(a,Cs_end);
xlabel('age (years)')
title(['~~~Final Immun dist, dt=', num2str(dt)])
grid on

Cs_ave = FileRead('Output/','Cs','Rave');
subplot(2,2,2), plot(t,Cs_ave);
title(['~~~Pop. Immun vs time, dt=',num2str(dt)]);
xlabel('time (days)')
grid on

Cs = FileRead('Output/','Cs','Entire'); 
subplot(2,2,3), imagesc(t/365,a/365,Cs);
set(gca,'YDir','normal');
colorbar;
ylabel('age');
xlabel('time');
title('Total Immunity');

%% Impact of immunity on the sigmoids (rho, psi, phi)
Cs = FileRead('Output/','Cs','Entire'); 
figure_setups;
subplot(2,2,1), imagesc(t/365,a/365,sigmoid_prob(Cs, 'phi'));
set(gca,'YDir','normal');
colorbar;
ylabel('age');
xlabel('time');
title('$\phi$');

subplot(2,2,2), imagesc(t/365,a/365,sigmoid_prob(Cs, 'rho'));
set(gca,'YDir','normal');
colorbar;
ylabel('age');
xlabel('time');
title('$\rho$');

subplot(2,2,3), imagesc(t/365,a/365,sigmoid_prob(Cs, 'psi'));
set(gca,'YDir','normal');
colorbar;
ylabel('age');
xlabel('time');
title('$\psi$');

figure_setups;
plot(a,P.rho)
axis_years(gca,tfinal)
title('$\rho$(age) at tfinal')
grid on

