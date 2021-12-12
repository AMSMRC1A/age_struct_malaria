%% plot Fig4 a & b - bifurcation plots -- results from "run_bifurcation_calcs.m"
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
if immunity_feedback == 0
    % average populational sigmoids f0 = f1 = average
    P.phi_f_0 = 0.915792480087329; % value at zero
    P.phi_f_1 = 0.915792480087329; % value at L (function saturates to this value)

    P.rho_f_0 = 0.114825053290306; % value at zero
    P.rho_f_1 = 0.114825053290306; % value at L (function saturates to this value)

    P.psi_f_0 = 0.114825053290306; % value at zero
    P.psi_f_1 = 0.114825053290306; % value at L (function saturates to this value)
    param = [0.01, 0.06, 0.11, 0.13, 0.16:0.05:1.0].^2; % max R0 < 7
else
    param = [0.01:0.025:0.525].^2; % max R0 = 7.3
end
%% options
plot_equilibrium = 0; % can set to zero if working with the DFE

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
    R0_cal()
    %% solve for EE
    if R0_cal()>1
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
R0_list = NaN(size(param));
for i = 1:length(param)
    P.(lP) = param(i);
    R0_list(i) = R0_cal();
end
%% extend data to include (R0=1,EE=0) point
[R0_list,ind] = sort([1,R0_list]);
re_max_DFE = [0,re_max_DFE]; re_max_DFE = re_max_DFE(ind);
re_max_EE = [0,re_max_EE]; re_max_EE = re_max_EE(ind);
I_frac_DFE = [0,I_frac_DFE]; I_frac_DFE = I_frac_DFE(ind);
I_frac_EE = [0,I_frac_EE]; I_frac_EE = I_frac_EE(ind);
D_frac_EE = [0,D_frac_EE]; D_frac_EE = D_frac_EE(ind);
A_frac_EE = [0,A_frac_EE]; A_frac_EE = A_frac_EE(ind);

%% Plot R0 as the bifurcation parameter
figure_setups_2;
hold on;
R0_DFE_EE = [R0_list,R0_list];
ind_stable = find([re_max_DFE,re_max_EE]<=0);
ind_unstable = find([re_max_DFE,re_max_EE]>=0);
I_frac = [I_frac_DFE,I_frac_EE];
h1 = plot(R0_DFE_EE(ind_stable), I_frac(ind_stable),'-','Marker','.','MarkerSize',30);
h4 = plot(R0_list, D_frac_EE,'-','Marker','o','MarkerSize',10);
h5 = plot(R0_list, A_frac_EE,'-','Marker','^','MarkerSize',10);
h2 = plot(R0_DFE_EE(ind_unstable), I_frac(ind_unstable),'r.','Marker','^','MarkerSize',5);
grid on; grid minor
xlabel('$\mathcal{R}_0(\beta_M)$');
ylabel('Fraction of population');
if immunity_feedback == 1
    title('Dynamic Immune Feedback');
else
    title('Constant High Immunity'); % could be constant low immunity case
end
P.betaM = 0.25; % "baseline value" of betaM
[~,~,D,A,~,~,~] = steady_state('EE','fsolve');
R0_baseline = R0_cal();
h3 = plot([R0_baseline,R0_baseline],[0,1],'m-');
axis([0 7.2 0 1])
xlabel('$\mathcal{R}_0 (\beta_M)$');
set(gcf, 'Renderer', 'Painters');
if immunity_feedback == 1
    legend([h1 h4 h5 h2 h3], {'$D_H+A_H$','$D_H$','$A_H$','Unstable','Baseline'},'Location','nw')
    print('Figures/result_bifur.eps', '-depsc','-r300')
else
    legend([h1 h4 h5 h2 h3], {'$D_H+A_H$','$D_H$','$A_H$','Unstable','Baseline'},'Location','e')
    print('Figures/result_bifur_0.eps', '-depsc','-r300')
end
toc