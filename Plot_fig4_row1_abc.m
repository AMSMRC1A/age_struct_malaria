%% plot Fig4 bifurcation plots using betaM -- results from "run_bifurcation_calcs.m"
tic 
clear all
close all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

%% numerical config
tfinal = 100*365; age_max = 80*365; P.age_max = age_max;
dt = 100; % time/age step size in days, default = 50; could go dt = 200 (still robust)
da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

Malaria_parameters_baseline; % model parameters - rates are in 1/day
lP = 'betaM';  % bifurcating parameters
immunity_feedback = 0;
param =[0.01:0.01:0.65].^2;
if immunity_feedback == -1 % betaM = 0.008
    P.phi_f_0 = 0.230971153687268; % value at zero
    P.phi_f_1 = 0.230971153687268; % value at L (function saturates to this value)
    P.rho_f_0 = 0.907631179204690; % value at zero
    P.rho_f_1 = 0.907631179204690 ; % value at L (function saturates to this value)
    P.psi_f_0 = 0.907631179204690 ; % value at zero
    P.psi_f_1 = 0.907631179204690 ; % value at L (function saturates to this value)
    beta0 = 0.005517283457290; % gives R0=1;
    beta_cut = 0.088276535316648; % gives R0=4;
elseif immunity_feedback == 0 % betaM = 0.25
    % average populational sigmoids f0 = f1 = average
    P.phi_f_0 = 0.915792480087329; % value at zero
    P.phi_f_1 = 0.915792480087329; % value at L (function saturates to this value)
    P.rho_f_0 = 0.114825053290306; % value at zero
    P.rho_f_1 = 0.114825053290306; % value at L (function saturates to this value
    P.psi_f_0 = 0.114825053290306; % value at zero
    P.psi_f_1 = 0.114825053290306; % value at L (function saturates to this value)
    beta0 = 0.021554911167740; % gives R0=1;
    beta_cut = 0.344878578683840; % gives R0=4;
elseif immunity_feedback == 1
    beta0 = 0.005206041495315; % gives R0=1;
    beta_cut = 0.083296663925047; % gives R0=4;
end
%% options
re_max_DFE = NaN(1,length(param));
re_max_EE = NaN(1,length(param));

A_frac_EE = NaN(1,length(param));
D_frac_EE = NaN(1,length(param));

I_frac_DFE = NaN(1,length(param));
I_frac_EE = NaN(1,length(param));


for i = 1:length(param)
    disp(['progress = ',num2str((i-1)/(length(param))*100),'%']);
    P.(lP) = param(i);
    Malaria_parameters_transform;
    %% solve for EE
    if P.(lP)>beta0 % R0>1
        FileName = ['Results/Bifur/',num2str(immunity_feedback),'/EE_',num2str(param(i),'%2.4f'),'.mat'];
        S = load(FileName,'x_EE','ee');
        x_EE = S.x_EE;
        ee = S.ee;
        I_frac_EE(i) = 1-da*trapz((x_EE(:,1)+x_EE(:,2)).*P.PH_stable);
        A_frac_EE(i) = da*trapz((x_EE(:,4)).*P.PH_stable);
        D_frac_EE(i) = da*trapz((x_EE(:,3)).*P.PH_stable);
        re_max_EE(i) = max(real(ee));
    end
    %% Solve for DFE
    FileName = ['Results/Bifur/',num2str(immunity_feedback),'/DFE_',num2str(param(i),'%2.4f'),'.mat'];
    S = load(FileName,'x_DFE','ee');
    x_DFE = S.x_DFE;
    ee = S.ee;
    I_frac_DFE(i) = 1 - da*trapz((x_DFE(:,1)+x_DFE(:,2)).*P.PH_stable);
    re_max_DFE(i) = max(real(ee));
end
%% Plot the results
% extend data to include (R0=1,EE=0) point -> beta0 value
[beta_list,ind] = sort([beta0,param]);
re_max_DFE = [0,re_max_DFE]; re_max_DFE = re_max_DFE(ind);
re_max_EE = [0,re_max_EE]; re_max_EE = re_max_EE(ind);
I_frac_DFE = [0,I_frac_DFE]; I_frac_DFE = I_frac_DFE(ind);
I_frac_EE = [0,I_frac_EE]; I_frac_EE = I_frac_EE(ind);
D_frac_EE = [0,D_frac_EE]; D_frac_EE = D_frac_EE(ind);
A_frac_EE = [0,A_frac_EE]; A_frac_EE = A_frac_EE(ind);

%% Plot with betaM on the x-axis
figure_setups_3;
hold on;
beta_DFE_EE = [beta_list,beta_list];
ind_stable = find([re_max_DFE,re_max_EE]<=0);
ind_unstable = find([re_max_DFE,re_max_EE]>=0);
I_frac = [I_frac_DFE,I_frac_EE];
h1 = plot(beta_DFE_EE(ind_stable), I_frac(ind_stable),'-','Marker','none','MarkerSize',30,'MarkerIndices',1:2:length(ind_stable));
h4 = plot(beta_list, D_frac_EE,'-','Marker','o','MarkerSize',10,'MarkerIndices',1:2:length(beta_list));
h5 = plot(beta_list, A_frac_EE,'-','Marker','^','MarkerSize',10,'MarkerIndices',1:2:length(beta_list));
h2 = plot(beta_DFE_EE(ind_unstable), I_frac(ind_unstable),'r:','Marker','^','MarkerSize',5,'MarkerIndices',1:2:length(ind_unstable));
grid on; grid minor
xlabel('$\beta_M$');
ylabel('Fraction of population');
axis([0 max(beta_list) 0 1])
% plot baseline
P.betaM = 0.25;
h3 = plot([P.betaM,P.betaM],[0,1],'m-');

set(gcf, 'Renderer', 'Painters');
if immunity_feedback == 1    
    h6 = plot([beta_cut,beta_cut],[0,1],'k-.'); [0.008, 0.03, 0.25]
    h7 = plot([0.008,0.008],[0,1],'g:');
    plot([0.03,0.03],[0,1],'g:');
    plot([0.25,0.25],[0,1],'g:');
    title('Dynamic immunity');
    legend([h6 h7], {'Row 2 cut','Row 3 cut'},'Location','e','Position', [0.615757275725849,0.638442852279118,0.283963226326078,0.162771438298906])
    print('Figures/result_bifur_1.eps', '-depsc','-r300')
elseif immunity_feedback == 0
    h6 = plot([beta_cut,beta_cut],[0,1],'k-.');
    title('Fixed high immunity');    
    print('Figures/result_bifur_0.eps', '-depsc','-r300')
elseif immunity_feedback == -1
    h6 = plot([beta_cut,beta_cut],[0,1],'k-.');
    title('Fixed low immunity');
    legend([h1 h4 h5 h2 h3], {'$D_H+A_H$','$D_H$','$A_H$','Unstable','Baseline'},'Location','e')
    print('Figures/result_bifur_-1.eps', '-depsc','-r300')
end
toc