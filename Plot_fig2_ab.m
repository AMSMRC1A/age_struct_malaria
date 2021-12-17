%% plot figure 2 a & b - results from "Malaria_parameters_transform.m"
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
Malaria_parameters_baseline;

%% figure 2a: plot calibrated fertility & death
figure_setups_3; hold on
grid on
plot(a/365, P.muH,'-','Color',colour_mat1)
plot(a/365, P.gH,'-.','Color',colour_mat2)
axis([0 age_max/365 0 1.5*10^-3]);
xlabel('Age (years)')
ylabel('Daily rate')
legend('Mortality rate $\mu_H(\alpha)$','Fertility rate $g_H(\alpha)$','Location','n')
%title(['Balanced demographics']);
set(gcf, 'Renderer', 'Painters');
print('Figures/demo_fer_mor.eps', '-depsc','-r300')
%% figure 2b: plot stable age distribution PH
figure_setups_3;
plot(a/365,P.PH_stable*365,'-k');
legend('$P_H(\alpha)$')
xlabel('Age (years)');
ylabel('Population density')
grid on
axis([0 age_max/365 0 max(P.PH_stable*365)]);
%title(['Stable age distribution']);
set(gcf, 'Renderer', 'Painters');
print('Figures/demo_PH.eps', '-depsc','-r300')