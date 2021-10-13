clear all;
% close all;
% clc;
format long;
global P lP
global F

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

P.a = a;
P.na = na;
P.nt = nt;
P.dt = dt;
P.da = da;
P.t = t;

%% load data
s =  load('Data/F.mat','F');
F = s.F;

%% SA setting
Malaria_parameters_baseline;
lP = 'none';
options = optimset('Display','iter','TolX',10^-3,'MaxIter',30);
x0 = [P.c_rho, P.c_psi, P.c_phi,P.k_rho, P.k_psi, P.k_phi];
% fit result - cmin cmax format 1.426205877278959  10.274569435455705  13.194094606343317   1.273174804387240   4.015779931093367   3.411148478637459
% residual ~ 8
lb = [1, 1, 1, 1, 1, 1];
ub = [500, 500, 500, 10, 10, 10];
[x,fval] = fmincon(@fun,x0,[],[], [], [], lb, ub, [], options);

keyboard

[SH,EH,DH,AH,~,~,Ctot] = steady_state('EE');
NH = trapz(SH+EH+DH+AH)*P.da;
NM = P.gM/P.muM;
[bH,bM] = biting_rate(NH,NM);
Lambda_M = bM*trapz(P.betaD*DH + P.betaA*AH)*P.da;
IM_frac_EE = P.sigma/(P.sigma+P.muM)*(Lambda_M/(Lambda_M + P.muM));
EIR = bH*IM_frac_EE*365
%% plotting rho sigmoid result
figure_setups;
z1 = sigmoid_prob(Ctot./P.PH_stable, 'rho');
plot(P.a/365,z1)
% figure_setups;
% plot(P.a/365,Ctot)
keyboard
%% plotting heatmap
% numerical config
tfinal = 100*365; % final time in days
age_max = 80*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
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

Malaria_parameters_baseline;
% P.c_rho = x(1);
% P.c_psi = x(2);
% P.c_phi = x(3);
% P.k_rho = x(4);
% P.k_psi = x(5);
% P.k_phi = x(6);
% Malaria_parameters_transform;
EIR_var = 'betaM'; % use this parameter to adjust the EIR
var_list = [0.006, 0.007, 0.01, 0.02, 0.03, 0.08, 0.1, 0.25];
final_immunity = zeros(na,length(var_list));
final_EIR = zeros(1,length(var_list));
for jj = 1:length(var_list)
    P.(EIR_var) = var_list(jj);
    Malaria_parameters_transform;
    [SH,EH,DH,AH,~,~,Ctot] = steady_state('EE');
    NH = trapz(SH+EH+DH+AH)*P.da;
    NM = P.gM/P.muM;
    [bH,bM] = biting_rate(NH,NM);
    Lambda_M = bM*trapz(P.betaD*DH + P.betaA*AH)*P.da;
    IM_frac_EE = P.sigma/(P.sigma+P.muM)*(Lambda_M/(Lambda_M + P.muM));   
    final_EIR(1,jj) = bH*IM_frac_EE*365; % annual EIR
    final_immunity(:,jj) = Ctot./P.PH_stable;
end
figure_setups;
imagesc(a/365,final_EIR,final_immunity');
xlim([0 10])
xlabel('age')
ylabel('EIR')
title('Immunity levels');
set(gca,'YDir','normal');
grid on
colormap jet
colorbar;
%% plot sigmoids
Malaria_parameters_baseline;
% P.phi_t_2 = x(1);
% P.phi_s_2 = x(2);
% P.rho_t_2 = x(3);
% P.rho_s_2 = x(4); 
% P.psi_t_2 = x(5);
% P.psi_s_2 = x(6);
% P.L = x(7);
figure_setups; hold on
[~,~,~,~,~,~,Ctot] = steady_state('EE');
cc_temp = Ctot./P.PH_stable;
cc = linspace(min(cc_temp),max(cc_temp),100);
rho_curve = sigmoid_prob(cc, 'rho');
phi_curve = sigmoid_prob(cc, 'phi');
psi_curve = sigmoid_prob(cc, 'psi');
plot(cc,rho_curve)
plot(cc,phi_curve)
plot(cc,psi_curve)
legend('rho','phi','psi')
axis([min(cc_temp) max(cc_temp) 0 1])

