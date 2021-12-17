clear all;
close all;
clc;
format long;
global P lP
global F

% numerical config
tfinal = 100*365; age_max = 80*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

%% load data - interpolated function
s =  load('Filipe_paper/F_Filipe.mat','F');
F = s.F;

%% optimization - fit to data 
Malaria_parameters_baseline;
options = optimset('Display','iter','TolX',10^-5,'MaxIter',20);
% [phi_s phi_r rho_s=psi_s rho_r=psi_r]
x0 = [2.432473210664639   1.277554702429119   3.186642715992263   1.030298116795388]; 
%====L=20=============
% res = 4.332481e-01; 31 iter, x = [0.123194800083748   0.069009892471593   0.157791442612232   0.051328791141945]; IC = [0.2 0.05 0.1 0.05];
%====L=10=============
% res = 4.376491e-01; 33 iter, x = [0.243400495658224   0.128546913213757   0.318633025120743   0.102982948018487]; IC = [0.3, 0.05, 0.3, 0.05];
%====L=5=============
% res = 4.376486e-01; 36 iter, x = [0.486500373585591   0.255614405471651   0.637327792072830   0.206052605587558]; IC = [0.5 0.1 0.5 0.1];
%====two parameter sigmoid t*L = t2 s*L=s2
% res = 4.376486e-01; x = [2.432431947045749   1.278072983365070   3.186658383357816   1.030263636242633]; IC - fitted values with L
lb = [0, 0.1, 0, 0.1];
ub = [8, 10, 8, 10];

[x,fval] = fmincon(@fun_Filipe,x0,[],[], [], [], lb, ub, [], options);

keyboard
% ----> update Malaria_parameters_baseline.m with the fitted results <-----

%% plot sigmoids with the populational average (in legend)
tfinal = 100*365; age_max = 80*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

% x = [2.432431947045749   1.278072983365070   3.186658383357816   1.030263636242633];
Malaria_parameters_baseline;
P.betaM = 0.008;
% P.phi_s_2 = x(1);
% P.phi_r_2 = x(2); 
% P.rho_s_2 = x(3);
% P.rho_r_2 = x(4); 
% P.psi_s_2 = x(3);
% P.psi_r_2 = x(4);
Malaria_parameters_transform;

figure_setups; hold on
[~,~,~,~,~,~,Ctot] = steady_state('EE','fsolve');
Ctot_pp = Ctot./P.PH_stable;
cc = linspace(0,max(Ctot_pp),100);
phi_curve = sigmoid_prob(cc, 'phi');
rho_curve = sigmoid_prob(cc, 'rho');
psi_curve = sigmoid_prob(cc, 'psi');
plot(cc,phi_curve,'-')
plot(cc,rho_curve,'--')
% plot(cc,psi_curve,'-.')
% population average sigmoids
phi_ave = trapz(sigmoid_prob(Ctot_pp, 'phi').*P.PH_stable)*P.da;
rho_ave = trapz(sigmoid_prob(Ctot_pp, 'rho').*P.PH_stable)*P.da;
psi_ave = trapz(sigmoid_prob(Ctot_pp, 'psi').*P.PH_stable)*P.da;
legend('$\phi(\tilde{C}_{H})$','$\rho(\tilde{C}_{H})=\psi(\tilde{C}_{H})$','Location','e')
axis([0 max(Ctot_pp) 0 1])
xlabel('$\tilde{C}_{H}$')
ylabel('Probability')
%title('Calibrated linking functions')
[phi_ave rho_ave psi_ave]
keyboard

%% plotting heatmap (age, EIR, immunity level)
tfinal = 100*365; age_max = 80*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

Malaria_parameters_baseline;

EIR_var = 'betaM'; % use this parameter to adjust the EIR
immunity_feedback = 1;
if immunity_feedback == 0 
    % population average sigmoids: f0 = f1 = average
    P.phi_f_0 = 0.915792480087329; % value at zero
    P.phi_f_1 = 0.915792480087329; % value at L (function saturates to this value)
    
    P.rho_f_0 = 0.114825053290306; % value at zero
    P.rho_f_1 = 0.114825053290306; % value at L (function saturates to this value)  
    
    P.psi_f_0 = 0.114825053290306; % value at zero
    P.psi_f_1 = 0.114825053290306; % value at L (function saturates to this value)   
    var_list = [0.01:0.05:1.0].^2; 
else
    var_list = [0.01:0.05:1].^2;
end
final_immunity = zeros(na,length(var_list));
final_EIR = zeros(1,length(var_list));
final_rho = zeros(na,length(var_list));
for jj = 1:length(var_list)
    P.(EIR_var) = var_list(jj);
    Malaria_parameters_transform;
    [SH,EH,DH,AH,~,~,Ctot] = steady_state('EE','fsolve');
    NH = trapz(SH+EH+DH+AH)*P.da;
    NM = P.gM/P.muM;
    [bH,bM] = biting_rate(NH,NM);
    Lambda_M = bM*trapz(P.betaD*DH + P.betaA*AH)*P.da;
    IM_frac_EE = P.sigma/(P.sigma+P.muM)*(Lambda_M/(Lambda_M + P.muM));   
    final_EIR(1,jj) = bH*IM_frac_EE*365; % annual EIR
    final_immunity(:,jj) = Ctot./P.PH_stable;
    final_rho(:,jj) = sigmoid_prob(Ctot./P.PH_stable, 'rho');
end
% save(['Results/Heatmap.mat'],'var_list','final_EIR','final_immunity')
% plot heatmap of final immunity(age, EIR)
figure_setups; 
imagesc(a/365,final_EIR,final_immunity');%
xlim([0 20])
xlabel('Age (years)')
ylabel('aEIR')
% title(['Immunity levels, feedback = ',num2str(immunity_feedback)]);
title('Immunity level per person');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:2:10);

%% plot heatmap of rho (age, EIR) = prob E->D
figure_setups; 
imagesc(a/365,final_EIR,final_rho');% ,[0 8.2] % specify the range of cvalues
xlim([0 10])
xlabel('age')
ylabel('EIR')
title(['rho (prob of E-D), feedback = ',num2str(immunity_feedback)]);
set(gca,'YDir','normal');
grid on
colormap jet
colorbar('Ticks',0:0.1:1);