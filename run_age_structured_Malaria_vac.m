clear all
close all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

tic

%% numerical config
age_max = 80*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.dt = dt; 
P.a = a;
P.na = na;
P.da = da;

% model parameters
Malaria_parameters_baseline;
% P.betaM = 0.25; % high EIR
% P.betaM = 0.01; % low EIR
% Malaria_parameters_transform;

immunity_feedback = 1;
if immunity_feedback == 0
    P.phi_f_0 = 0.570320665853183; % value at zero
    P.phi_f_1 = 0.570320665853183; % value at L (function saturates to this value)
    
    P.rho_f_0 = 0.088575583518581; % value at zero
    P.rho_f_1 = 0.088575583518581; % value at L (function saturates to this value)  
    
    P.psi_f_0 = 0.409302219871934; % value at zero
    P.psi_f_1 = 0.409302219871934; % value at L (function saturates to this value)   
end

% [SH, EH, DH, AH, VH, SM, EM, IM, Cm, Cac, Cv, Ctot] = age_structured_Malaria_IC_vac('EE');
% PH = SH+EH+DH+AH+VH;
% PH_final = PH(:,end); % total human at age a, t = n
% NH = trapz(PH,1)*da;
%% initial condition 'init' 'EE'
[SH0, EH0, DH0, AH0, VH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0] = age_structured_Malaria_IC_vac('EE');
%% time evolution - initial run
tfinal = 1*365; t = (0:dt:tfinal)'; nt = length(t);
P.nt = nt;  P.t = t;
[SH, EH, DH, AH, VH, SM, EM, IM, Cm, Cac, Cv, Ctot] = age_structured_Malaria_vac(da,na,tfinal,...
    SH0, EH0, DH0, AH0, VH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
PH = SH+EH+DH+AH+VH;
PH_final = PH(:,end); % total human at age a, t = n
NH = trapz(PH,1)*da;
vacc = trapz(P.vp.*SH,1)*P.da*365*P.NN/1000;
%% time evolution - continuous run
% non-vac control
% cont_level = 0.5;
% P.betaA = (1-cont_level)*P.betaA;
% P.betaD = (1-cont_level)*P.betaD;
% P.betaM = (1-cont_level)*P.betaM;
% Malaria_parameters_transform;
% vac control
% Malaria_parameters_transform_vac; 
tfinal_conti = 200*365; t2 = (tfinal:dt:tfinal+tfinal_conti)'; nt = length(t2);
P.nt = nt;  P.t = t2;
SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end);
SM0 = SM(end); EM0 = EM(end); IM0 = IM(end); 
Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
[SH2, EH2, DH2, AH2, VH2, SM2, EM2, IM2, Cm2, Cac2, Cv2, Ctot2] = age_structured_Malaria_vac(da,na,tfinal_conti,...
    SH0, EH0, DH0, AH0, VH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
PH2 = SH2+EH2+DH2+AH2+VH2;
NH2 = trapz(PH2,1)*da;
vacc2 = trapz(P.vp.*SH2,1)*P.da*365*P.NN/1000;
t = [t;t2];
SH = [SH,SH2]; EH = [EH,EH2]; DH = [DH, DH2]; AH = [AH, AH2]; VH = [VH, VH2];
SM = [SM, SM2]; EM = [EM, EM2]; IM = [IM, IM2];
Cm = [Cm, Cm2]; Cac = [Cac, Cac2]; Cv = [Cv, Cv2]; Ctot = [Ctot,Ctot2];
PH = SH+EH+DH+AH+VH;
PH_final = PH(:,end); % total human at age a, t = n
NH = [NH, NH2];
vacc = [vacc, vacc2];
%% output data to .mat file for analysis
SH_EE = SH(:,end); EH_EE = EH(:,end); AH_EE = AH(:,end); DH_EE = DH(:,end); VH_EE = VH(:,end); PH_EE = PH(:,end); vp = P.vp;
Cm_EE = Cm(:,end); Cac_EE = Cac(:,end); Ctot_EE = Ctot(:,end);
save(['Results/Vaccine/vp0_',num2str(P.vp0*100),'.mat'],'t','a','vp','vacc','SH_EE','EH_EE','AH_EE','DH_EE','VH_EE','PH_EE',...
    'Cm_EE','Cac_EE','Ctot_EE');
%% EIR
% [bh,bm] = biting_rate(NH,NM);
% EIR = bh.*IM./NM*365;
% EIR_EE = EIR(end)
% tic
R0 = R0_cal()
% keyboard
% toc

% figure_setups;
% plot(t,EIR,'b-'); hold on;
%% vaccine #
figure_setups; hold on
grid on
plot(t/365,vacc)
xlim([0 max(t)/365]);
%% Population size versus time
figure_setups;
plot(t/365,trapz(SH,1)*da,'-','Color',colour_mat1); hold on;
plot(t/365,trapz(EH,1)*da,'--','Color',colour_mat3);
plot(t/365,trapz(AH,1)*da,'-.','Color',colour_mat2);
plot(t/365,trapz(DH,1)*da,'-','Color',colour_mat7);
plot(t/365,trapz(VH,1)*da,'-','Color',colour_mat6);
plot(t/365,NH,'-.k')
legend('SH-age','EH-age','AH-age', 'DH-age','VH-age','$N_H$','Location','e');
title(['Population size vs time', '~~feedback = ',num2str(immunity_feedback)]); 
grid on; grid minor
axis([0 max(t)/365 0 max(NH)+0.1]);
%% Age profiles at tfinal
figure_setups;
plot(a/365,SH(:,end),'-','Color',colour_mat1); hold on;
plot(a/365,EH(:,end),'--','Color',colour_mat3);
plot(a/365,DH(:,end),'-.','Color',colour_mat2);
plot(a/365,AH(:,end),':','Color',colour_mat7);
plot(a/365,VH(:,end),':','Color',colour_mat6);
plot(a/365,PH_final,'-k');
legend('SH','EH','DH', 'AH','VH','PH');
% title(['Final Age Dist.,~~ feedback =',num2str(immunity_feedback)]);
title(['Final Age Distribution']);
xlabel('age (years)');
grid on
axis([0 age_max/365 0 max(PH_final)]);
%% Age proportions at tfinal prop
figure_setups;
plot(a/365,SH(:,end)./PH_final,'-','Color',colour_mat1); hold on;
plot(a/365,EH(:,end)./PH_final,'--','Color',colour_mat6);
plot(a/365,DH(:,end)./PH_final,'-.','Color',colour_mat2);
plot(a/365,AH(:,end)./PH_final,':','Color',colour_mat3);
plot(a/365,VH(:,end)./PH_final,':','Color',colour_mat5);
plot(a/365,(AH(:,end)+DH(:,end))./PH_final,'r-.');
plot(a/365,PH_final./PH_final,'-k');
legend('SH/PH','EH/PH','DH/PH', 'AH/PH', 'VH/PH','(AH+DH)/PH');
title(['Final Age Dist. Proportions']); 
% title(['Final Age Dist. Proportions ~~ feedback =',num2str(immunity_feedback)]); 
xlabel('age (years)');
grid on
axis([0 P.age_max/365 0 1.1]);
xlim([0 30])
%% Population proportions versus time
% figure_setups;
% plot(t,trapz(SH,1)*da./NH,'-','Color',colour_mat1); hold on;
% plot(t,trapz(EH,1)*da./NH,'--','Color',colour_mat3);
% plot(t,trapz(AH,1)*da./NH,'-.','Color',colour_mat2);
% plot(t,trapz(DH,1)*da./NH,'-','Color',colour_mat7);
% plot(t,trapz(VH,1)*da./NH,'-','Color',colour_mat6);
% plot(t,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1)+trapz(VH,1))*da./NH,'-.k');
% legend('SH-age','EH-age','AH-age', 'DH-age', 'VH-age','$N_H$');
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
% figure_setups;
% plot(a/365,Cac(:,end),'-.r');
% hold on;
% plot(a/365,Cm(:,end),'-.b');
% plot(a/365,Cv(:,end),'-.c');
% plot(a/365,Ctot(:,end),'-.k');
% xlabel('age (years)')
% legend('Acquired','Maternal','Vaccine-derived','Total','Location','SouthEast');
% title(['Immun dist.~~ feedback =',num2str(immunity_feedback)]);
% axis([0 age_max/365 0 max(Ctot(:,end))*1.1]);
% grid on
%%
% figure_setups;
% plot(a/365,Cac(:,end)./PH_final,'-.');
% hold on;
% plot(a/365,Cm(:,end)./PH_final,'--');
% plot(a/365,Cv(:,end)./PH_final,'--');
% plot(a/365,Ctot(:,end)./PH_final,'-');
% xlabel('age (years)')
% ylabel('immunity level')
% legend('Acquired (pp)','Maternal (pp)','Vaccine-derived (pp)','Total (pp)','Location','SouthEast');
% % title(['Per-person Immun dist.~~ feedback =',num2str(immunity_feedback)]);
% title('Per-person Immun distribution');
% axis([0 age_max/365 0 max(Ctot(:,end)./PH_final)*1.1]);
% xlim([0 10])
% ylim([0 7])
% grid on
%% plot sigmoids
% figure_setups; hold on;
% plot(a/365,sigmoid_prob(Ctot(:,end)./PH_final, 'rho'),'-');
% % plot(a/365,sigmoid_prob(Ctot(:,end)./PH_final, 'psi'));
% % plot(a/365,sigmoid_prob(Ctot(:,end)./PH_final, 'phi'));
% grid on
% legend('rho (suscept.)');
% axis([0 age_max/365 0 1]);
% title(['EIR = ',num2str(EIR_EE)])
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
% toc