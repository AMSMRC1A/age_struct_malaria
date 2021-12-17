%% plot Fig3 b & c - Immunity curves at low & high setting -- results from "run_parameter_fit.m"

close all
clear all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2
%% numerical config
tfinal = 300*365; % final time in days
age_max = 80*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

Malaria_parameters_baseline;
P.betaM = 0.008; %low aEIR
Malaria_parameters_transform;
%% initial condition 'init' 'EE'
[SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
%% time evolution
[SH, EH, DH, AH, SM, EM, IM, Cm, Cac, Ctot] = age_structured_Malaria(da,na,tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
PH = SH+EH+DH+AH;
PH_final = PH(:,end); % total human at age a, t = n
NH = trapz(PH,1)*da;
NM = SM+EM+IM;

%% EIR
[bh,bm] = biting_rate(NH,NM);
EIR = bh.*IM./NM*365;
EIR_EE = EIR(end)
tic
R0 = R0_cal()
toc

%% Immunity breakdown - per person
figure_setups_3;
plot(a/365,Cac(:,end)./PH_final,'-.');
hold on;
plot(a/365,Cm(:,end)./PH_final,'--');
plot(a/365,Ctot(:,end)./PH_final,'-');
xlabel('Age (years)')
ylabel('Immunity level')
legend('$\widetilde{C}_e$ (Exposure-acquired immunity)','$\widetilde{C}_{m}$ (Maternal immunity)','$\widetilde{C}_{H}$ (Total immunity)','Location','SouthEast');
title('Per-person immunity dist.');
axis([0 age_max/365 0 max(Ctot(:,end)./PH_final)*1.1]);
xlim([0 10])
ylim([0 7])
grid on
set(gcf, 'Renderer', 'Painters');
if P.betaM == 0.008
    legend('$\widetilde{C}_e$ (Exposure-acquired immunity)','$\widetilde{C}_{m}$ (Maternal immunity)','$\widetilde{C}_{H}$ (Total immunity)','Location','SouthEast');
    print('Figures/result_immunity_low.eps', '-depsc','-r300')
else
    print('Figures/result_immunity_high.eps', '-depsc','-r300')
end
