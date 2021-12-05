%% plot Fig2c fitted sigmoid curves -- results from "run_parameter_fit.m"

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
figure_setups_3; hold on
[~,~,~,~,~,~,Ctot] = steady_state('EE','fsolve');
Ctot_pp = Ctot./P.PH_stable;
cc = linspace(0,max(Ctot_pp),100);
phi_curve = sigmoid_prob(cc, 'phi');
rho_curve = sigmoid_prob(cc, 'rho');
psi_curve = sigmoid_prob(cc, 'psi');
plot(cc,phi_curve,'-')
plot(cc,rho_curve,'-.')
legend('$\phi(\widetilde{C}_{H})$','$\rho(\widetilde{C}_{H})=\psi(\widetilde{C}_{H})$','Location','e')
axis([0 max(Ctot_pp) 0 1])
xlabel('$\widetilde{C}_{H}$')
ylabel('Probability')
%title('Calibrated linking functions')
set(gcf, 'Renderer', 'Painters');
print('Figures/sigmoids_ms.eps', '-depsc','-r300')
