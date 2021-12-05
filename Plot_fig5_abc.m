%% plot Fig5 ab & c - Population density & proportion at low & high setting

close all
clear all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2
%% numerical config
tfinal = 100*365; % final time in days
age_max = 80*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

% model parameters
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

%% Age profiles at tfinal
figure_setups_3;
plot(a/365,SH(:,end),'-','Color',colour_mat1); hold on;
plot(a/365,EH(:,end),':','Color',colour_mat4);
plot(a/365,AH(:,end),'--','Color',colour_mat3);
plot(a/365,DH(:,end),'-.','Color',colour_mat2);
plot(a/365,PH_final,'-k');
legend('$S_H$','$E_H$','$A_H$', '$D_H$','$P_H$','location','e');
title('~~~~~~Final age distribution'); 
xlabel('Age (years)');
ylabel('Population density')
grid on
axis([0 age_max/365 0 max(PH_final)]);
set(gcf, 'Renderer', 'Painters');
if P.betaM == 0.008
    print('Figures/solu_final_age_num_low.eps', '-depsc','-r300')
else
    print('Figures/solu_final_age_num.eps', '-depsc','-r300')
end

%% Age proportions at tfinal prop
figure_setups_3;
plot(a/365,SH(:,end)./PH_final,'-','Color',colour_mat1); hold on;
plot(a/365,EH(:,end)./PH_final,':','Color',colour_mat4);
plot(a/365,AH(:,end)./PH_final,'--','Color',colour_mat3);
plot(a/365,DH(:,end)./PH_final,'-.','Color',colour_mat2);
plot(a/365,PH_final./PH_final,'-k');
legend('$\widetilde{S}_H$','$\widetilde{E}_H$','$\widetilde{A}_H$', '$\widetilde{D}_H$','$\widetilde{P}_H$','location','e');
title('Final age dist. proportion'); 
xlabel('Age (years)');
ylabel('Fraction of population')
grid on
axis([0 P.age_max/365 0 1.1]);
xlim([0 30])

set(gcf, 'Renderer', 'Painters');
if P.betaM == 0.008
    print('Figures/solu_final_age_prop_low.eps', '-depsc','-r300')
else
    print('Figures/solu_final_age_prop.eps', '-depsc','-r300')
end
