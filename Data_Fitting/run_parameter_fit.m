clear all;
close all;
clc;
format long;
global P lP
global F

% numerical config
tfinal = 100*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

%% load data - interpolated function
s =  load('Filipe_paper/F_Filipe.mat','F');
F = s.F;

%% optimization - fit to data 
Malaria_parameters_baseline;
options = optimset('Display','iter','TolX',10^-5,'MaxIter',10);
          
x0 = [0.469638663283222, 0.543314482330510,  ... % phi
      0.177378797967747, 0.220419188961217,  ... % rho = psi
      8.370304463522121 % L
   ]; 
% age_max = 100;
% res = 1.197; L = 10; x = [0.423736367029315   0.547583380663965 0.173640457234346   0.175176012516704]; not too sensitive to Initial guess
% res = 1.199; L = 5;  x = [0.388255983233911   0.629539286068975 0.358933326495389   0.368027194195680]; not too sensitive to Initial guess
% res = 1.208; L = 20; x = [0.170419947143750   0.172846316120512   0.084973685532207 0.090371375233933];
% res = 1.26; L vary, x = [0.469638663283222   0.543314482330510 0.177378797967747   0.220419188961217   8.370304463522121];
% res = 1.183; L vary, x = [0.200455981024660   0.741327847456465 0.225836736227440   0.227377462176943   7.974031085468034];
lb = [0, 0.01, 0, 0.01, 1];
ub = [1, 1, 1, 1, 30];

[x,fval] = fmincon(@fun_Filipe,x0,[],[], [], [], lb, ub, [], options);

keyboard
% ----> update Malaria_parameters_baseline.m with the fitted results <-----

%% plot sigmoids with the populational average (in legend)
tfinal = 100*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

x = [0.200455981024660   0.741327847456465   0.225836736227440   0.227377462176943   7.974031085468034];
Malaria_parameters_baseline;
P.phi_t_2 = x(1);
P.phi_s_2 = x(2); 
P.rho_t_2 = x(3);
P.rho_s_2 = x(4); 
P.psi_t_2 = x(3);
P.psi_s_2 = x(4);
P.L = x(5);
Malaria_parameters_transform;

figure_setups; hold on
[~,~,~,~,~,~,Ctot] = steady_state('EE');
Ctot_pp = Ctot./P.PH_stable;
cc = linspace(0,max(Ctot_pp),100);
phi_curve = sigmoid_prob(cc, 'phi');
rho_curve = sigmoid_prob(cc, 'rho');
psi_curve = sigmoid_prob(cc, 'psi');
plot(cc,phi_curve,'-')
plot(cc,rho_curve,'--')
plot(cc,psi_curve,'-.')
% population average sigmoids
phi_ave = trapz(sigmoid_prob(Ctot_pp, 'phi').*P.PH_stable)*P.da;
rho_ave = trapz(sigmoid_prob(Ctot_pp, 'rho').*P.PH_stable)*P.da;
psi_ave = trapz(sigmoid_prob(Ctot_pp, 'psi').*P.PH_stable)*P.da;
legend(['$\phi$, pop ave = ',num2str(phi_ave,3)],['$\rho$, pop ave = ',num2str(rho_ave,3)],['$\psi$, pop ave = ',num2str(psi_ave,3)],'Location','e')
axis([0 max(Ctot_pp) 0 1])
xlabel('$\tilde{C}_{H}$')
ylabel('Probability')
[phi_ave rho_ave psi_ave]
keyboard

%% plotting heatmap (age, EIR, immunity level)
tfinal = 100*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

x = [0.357191327860647   0.706962687481590   0.194435605429801   0.165465456034519   0.094284616922569   0.434722350678851   8.577416380343156];
Malaria_parameters_baseline;
P.phi_t_2 = x(1);
P.phi_s_2 = x(2); 
P.rho_t_2 = x(3);
P.rho_s_2 = x(4); 
P.psi_t_2 = x(5);
P.psi_s_2 = x(6);
P.L = x(7);
Malaria_parameters_transform;

EIR_var = 'betaM'; % use this parameter to adjust the EIR
immunity_feedback = 1;
if immunity_feedback == 0
    % population average sigmoids: f0 = f1 = average
    P.phi_f_0 = 0.570320665853183; % value at zero
    P.phi_f_1 = 0.570320665853183; % value at L (function saturates to this value)
    
    P.rho_f_0 = 0.088575583518581; % value at zero
    P.rho_f_1 = 0.088575583518581; % value at L (function saturates to this value)  
    
    P.psi_f_0 = 0.409302219871934; % value at zero
    P.psi_f_1 = 0.409302219871934; % value at L (function saturates to this value)   
    var_list = [0.01:0.05:1.0].^2; 
else
    var_list = [0.01:0.025:0.55].^2;
end
final_immunity = zeros(na,length(var_list));
final_EIR = zeros(1,length(var_list));
final_rho = zeros(na,length(var_list));
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
    final_rho(:,jj) = sigmoid_prob(Ctot./P.PH_stable, 'rho');
end

% plot heatmap of final immunity(age, EIR)
figure_setups; 
imagesc(a/365,final_EIR,final_immunity');% ,[0 8.2] % specify the range of cvalues
xlim([0 10])
xlabel('age (years)')
ylabel('EIR')
% title(['Immunity levels, feedback = ',num2str(immunity_feedback)]);
title('Immunity level per person');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:1:8);

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