clear all
% close all
% clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

tic

%% numerical config
tfinal = 250*365; % final time in days
age_max = 70*365; % max ages in days
P.age_max = age_max;
dt = 50; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);

P.a = a;
P.na = na;
P.nt = nt;
P.dt = dt;
P.da = da;
P.t = t;

% model parameters 
Malaria_parameters_baseline;

%% time evolution
[SH, EH, DH, AH, SM, EM, IM, Cm, Cac, Ctot] = age_structured_Malaria(na,da,nt);
PH = SH+EH+DH+AH;
PH_final = PH(:,end); % total human at age a, t = n
NH = trapz(PH,1)*da;
NM = SM+EM+IM;
%% EIR
[bh,bm] = biting_rate(NH,NM);
EIR = bh.*IM./NM*365;
EIR_EE = EIR(end)
R0 = R0_cal()
% figure_setups;
% plot(t,EIR,'b-'); hold on;

%% Population size versus time
% figure_setups;
% plot(t,trapz(SH,1)*da,'-','Color',colour_mat1); hold on;
% plot(t,trapz(EH,1)*da,'--','Color',colour_mat3);
% plot(t,trapz(AH,1)*da,'-.','Color',colour_mat2);
% plot(t,trapz(DH,1)*da,'-','Color',colour_mat7);
% plot(t,NH,'-.k')
% legend('SH-age','EH-age','AH-age', 'DH-age','$N_H$','Location','NorthWest');
% title('Population size vs time');
% axis_years(gca,tfinal); % change to x-axis to years if needed
% grid on
% axis([0 tfinal 0 max(NH)+0.1]);
%% Age profiles at tfinal
figure_setups;
plot(a/365,SH(:,end),'-','Color',colour_mat1); hold on;
plot(a/365,EH(:,end),'--','Color',colour_mat3);
plot(a/365,AH(:,end),'-','Color',colour_mat2);
plot(a/365,DH(:,end),'-','Color',colour_mat7);
plot(a/365,PH_final,'-.k');
legend('SH','EH','AH', 'DH','PH');
title('Final Age Dist.');
xlabel('age');
grid on
axis([0 age_max/365 0 max(PH_final)]);
%% Age proportions at tfinal prop
% figure_setups;
% plot(a,SH(:,end)./PH_final,'-','Color',colour_mat1); hold on;
% plot(a,EH(:,end)./PH_final,'-','Color',colour_mat3);
% plot(a,AH(:,end)./PH_final,'-','Color',colour_mat2);
% plot(a,DH(:,end)./PH_final,'-','Color',colour_mat7);
% plot(a,PH_final./PH_final,'-k');
% legend('SH','EH','AH', 'DH','PH');
% title('Final Age Dist. Proportions');
% xlabel('age');
% axis_years(gca,age_max); % change to x-axis to years if needed
% grid on
% axis([0 age_max 0 1]);
%% Population proportions versus time
% figure_setups;
% plot(t,trapz(SH,1)*da./NH,'-','Color',colour_mat1); hold on;
% plot(t,trapz(EH,1)*da./NH,'--','Color',colour_mat3);
% plot(t,trapz(AH,1)*da./NH,'-.','Color',colour_mat2);
% plot(t,trapz(DH,1)*da./NH,'-','Color',colour_mat7);
% plot(t,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1))*da./NH,'-.k');
% legend('SH-age','EH-age','AH-age', 'DH-age','$N_H$');
% title('Population proportions vs time');
% axis_years(gca,tfinal); % change to x-axis to years if needed
% xlabel('time');
% grid on
% axis([0 tfinal 0 1.1]);
%% Immunity related figures
% figure_setups;
% subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/4))./PH(:,floor(nt/4)));
% hold on;
% subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/2))./PH(:,floor(nt/2)));
% subplot(2,2,1), plot(a/365,Ctot(:,floor(3*nt/4))./PH(:,floor(3*nt/4)));
% subplot(2,2,1), plot(a/365,Ctot(:,end)./PH(:,end),'-.');
% xlabel('age')
% legend(['t = ',num2str(tfinal/(4*365))],['t = ',num2str(tfinal/(2*365))],...
%     ['t = ',num2str(3*tfinal/(4*365))],['t = ',num2str(tfinal/365)],'Location','NorthWest');
% title('$C_{total}(t)/P_H(t)$');
% grid on
% subplot(2,2,2), plot(t/365,(trapz(Ctot,1)*da)./NH);
% title('$\int C_{total}(\alpha,t)d\alpha / N_H(t)$');
% xlabel('time');
% grid on
% 
% subplot(2,2,3), imagesc(t/365,a/365,Ctot./PH);
% set(gca,'YDir','normal');
% colorbar;
% ylabel('age');
% xlabel('time');
% title('$C_{total}(\alpha,t)/P_H(\alpha,t)$');
% 
% subplot(2,2,4), plot(a/365,Ctot(:,floor(nt/4)));
% hold on;
% subplot(2,2,4), plot(a/365,Ctot(:,floor(nt/2)));
% subplot(2,2,4), plot(a/365,Ctot(:,floor(3*nt/4)));
% subplot(2,2,4), plot(a/365,Ctot(:,end),'-.');
% xlabel('age')
% legend(['t = ',num2str(tfinal/(4*365))],['t = ',num2str(tfinal/(2*365))],...
%     ['t = ',num2str(3*tfinal/(4*365))],['t = ',num2str(tfinal/365)]);
% title('$C_{total}(t)$');
%% Immunity breakdown
figure_setups;
plot(a/365,Cac(:,end)./PH_final,'-.r');
hold on;
plot(a/365,Cm(:,end)./PH_final,'-.b');
plot(a/365,Ctot(:,end)./PH_final,'-.k');
xlabel('age (years)')
legend('Acquired','Maternal','Total','Location','SouthEast');
title('Immunity per person');
axis([0 age_max/365 0 max(Ctot(:,end)./PH_final)*1.1]);
grid on
%%
% subplot(1,2,2)
% plot(a/365,Cac(:,end)./PH_final,'-r');
% hold on;
% plot(a/365,Cm(:,end)./PH_final,'-b');
% plot(a/365,Ctot(:,end)./PH_final,'-k');
% xlabel('age (years)')
% legend('$C_{ac}(\alpha,t_{max})/P_H(\alpha,t_{max})$',...
%     '$C_{m}(\alpha,t_{max})/P_H(\alpha,t_{max})$','$C_{total}(\alpha,t_{max})/P_H(\alpha,t_{max})$',...
%     'Location','SouthEast');
% title('Immun dist.');
% axis([0 age_max/365 0 max(Ctot(:,end)./PH_final)*1.1]);
% grid on
%% Mosquito infection dynamics
% figure_setups;
% plot(t,SM,'b-'); hold on;
% plot(t,EM,'-','Color',colour_r1);
% plot(t,IM,'r-.');
% plot(t,SM+EM+IM,'-.')
% legend('SM','EM','IM','$N_M$');
% title('mosquito population size by stages')
% axis_years(gca,tfinal); % change to x-axis to years if needed
% grid on
% axis([0 tfinal 0 5])
toc
