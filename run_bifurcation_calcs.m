clear all
% close all
% clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

%% numerical config
tfinal = 100*365; age_max = 80*365; P.age_max = age_max;
dt = 100; % time/age step size in days, default = 50; could go dt = 200 (still robust)
da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);

P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

% load_data is an option to use previously saved solution data if available
% setting this to 1 will speed up runtime
load_data = 1;
Malaria_parameters_baseline; % model parameters - rates are in 1/day
lP = 'betaM';  % bifurcating parameters
immunity_feedback = 0; % -1 = low immunity fixed, 1 = dynamic, 0 = high fixed
comparison = 0; % plot multiple diagrams at the end if the data is there
param =[0.01:0.01:0.65].^2;
if immunity_feedback == -1 % betaM = 0.008
    P.phi_f_0 = 0.230971153687268; % value at zero
    P.phi_f_1 = 0.230971153687268; % value at L (function saturates to this value)
    P.rho_f_0 = 0.907631179204690; % value at zero
    P.rho_f_1 = 0.907631179204690 ; % value at L (function saturates to this value)
    P.psi_f_0 = 0.907631179204690 ; % value at zero
    P.psi_f_1 = 0.907631179204690 ; % value at L (function saturates to this value)
elseif immunity_feedback == 0 % betaM = 0.25
    % average populational sigmoids f0 = f1 = average
    P.phi_f_0 = 0.915792480087329; % value at zero
    P.phi_f_1 = 0.915792480087329; % value at L (function saturates to this value)
    P.rho_f_0 = 0.114825053290306; % value at zero
    P.rho_f_1 = 0.114825053290306; % value at L (function saturates to this value
    P.psi_f_0 = 0.114825053290306; % value at zero
    P.psi_f_1 = 0.114825053290306; % value at L (function saturates to this value)
elseif immunity_feedback == 1
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
    F_prop = @(x) human_model_der_prop(x);
    options = optimoptions('fsolve','Display','none','MaxIterations',50);
    tic
    R0_cal()
    %% solve for EE
    if R0_cal()>1
        FileName = ['Results/Bifur/',num2str(immunity_feedback),'/EE_',num2str(param(i),'%2.4f'),'.mat'];
        if exist(FileName,'file') && load_data == 1 % Q_val already calculated before
            S = load(FileName,'x_EE','ee');
            x_EE = S.x_EE;
            ee = S.ee;
        else
            [S,E,D,A,Cac,~,~] = steady_state('EE','fsolve');
            x0 = [S./P.PH_stable;E./P.PH_stable;D./P.PH_stable;A./P.PH_stable;Cac./P.PH_stable];
            [xsol,err,~,~,jacobian] = fsolve(F_prop,x0,options);
            x_EE = reshape(xsol,[P.na,5]);
            jacobian([1,P.na+1,2*P.na+1,3*P.na+1, 4*P.na+1],:)=0; % zero out the rows
            jacobian(:,[1,P.na+1,2*P.na+1,3*P.na+1, 4*P.na+1])=0; % zero out the columns
            ee = eig(jacobian);
            ee(ee==0)=[];
            save(FileName,'x_EE','ee')
        end
        if norm(F_prop(x_EE),Inf) > 10^-4
            disp('EE not achieved')
            keyboard
        end
        if plot_equilibrium == 1
            % plot proportion
            figure_setups; hold on;
            plot(a,x_EE(:,1),'-','Color',colour_mat1);
            plot(a,x_EE(:,2),'-','Color',colour_mat3);
            plot(a,x_EE(:,4),'-','Color',colour_mat2);
            plot(a,x_EE(:,3),'-','Color',colour_mat7);
            legend('SH (solver)','EH (solver)','AH (solver)', 'DH (solver)');
            title('Final Age Dist. prop');
            xlabel('age');
            axis_years(gca,age_max); % change to x-axis to years if needed
            grid on
            %axis([0 age_max 0 max(sum(x_EE,2))]);
            %keyboard
        end
        I_frac_EE(i) = 1-da*trapz((x_EE(:,1)+x_EE(:,2)).*P.PH_stable);
        A_frac_EE(i) = da*trapz((x_EE(:,4)).*P.PH_stable);
        D_frac_EE(i) = da*trapz((x_EE(:,3)).*P.PH_stable);
        re_max_EE(i) = max(real(ee));
        %     if re_max_EE(i) < 0
        %        disp(['max real part of eigenvalues = ',num2str(re_max_EE(i),'%10.6f'), '; EE is stable']);
        %     else
        %        disp(['max real part of eigenvalues = ',num2str(re_max_EE(i),'%10.6f'), '; EE is unstable']);
        %     end
    end
    toc
    tic
    %% Solve for DFE
    FileName = ['Results/Bifur/',num2str(immunity_feedback),'/DFE_',num2str(param(i),'%2.4f'),'.mat'];
    if exist(FileName,'file') && load_data == 1 % Q_val already calculated before
        S = load(FileName,'x_DFE','ee');
        x_DFE = S.x_DFE;
        ee = S.ee;
    else
        [~,~,~,~,Cac,~,~] = steady_state('DFE','numerical');
        x0 = [ones(length(a),1); 0*ones(length(a),1); 0*ones(length(a),1); 0*ones(length(a),1); Cac./P.PH_stable]; % initial guess for the DFE
        [xsol_DFE,err,~,~,jacobian] = fsolve(F_prop,x0,options);
        if max(max(abs(err)))>10^-5
            disp('not converged')
            keyboard
        end
        jacobian([1, P.na+1, 2*P.na+1, 3*P.na+1, 4*P.na+1],:) = 0; % zero out the rows
        jacobian(:,[1, P.na+1, 2*P.na+1, 3*P.na+1, 4*P.na+1]) = 0; % zero out the columns
        ee = eig(jacobian);
        ee(ee==0)=[];
        x_DFE = reshape(xsol_DFE,[P.na,5]);
        save(FileName,'x_DFE','ee')
    end
    if norm(F_prop(x_DFE),Inf) > 10^-5
        disp('DFE not achieved')
        keyboard
    end
    I_frac_DFE(i) = 1 - da*trapz((x_DFE(:,1)+x_DFE(:,2)).*P.PH_stable);
    re_max_DFE(i) = max(real(ee));
    %     if re_max_DFE(i) < 0
    %        disp(['max real part of eigenvalues = ',num2str(re_max_DFE(i),'%10.6f'), '; DFE is stable']);
    %     else
    %        disp(['max real part of eigenvalues = ',num2str(re_max_DFE(i),'%10.6f'), '; DFE is unstable']);
    %     end
    toc
end
%% Plot the results
R0_list = NaN(size(param));
for i = 1:length(param)
    P.(lP) = param(i);
    R0_list(i) = R0_cal();
end
%% extend data to include (R0=1,EE=0) point - could solve for this betaM numerically?
% [R0_list,ind] = sort([1,R0_list]);
% re_max_DFE = [0,re_max_DFE]; re_max_DFE = re_max_DFE(ind);
% re_max_EE = [0,re_max_EE]; re_max_EE = re_max_EE(ind);
% I_frac_DFE = [0,I_frac_DFE]; I_frac_DFE = I_frac_DFE(ind);
% I_frac_EE = [0,I_frac_EE]; I_frac_EE = I_frac_EE(ind);
% D_frac_EE = [0,D_frac_EE]; D_frac_EE = D_frac_EE(ind);
% A_frac_EE = [0,A_frac_EE]; A_frac_EE = A_frac_EE(ind);

%% Plot with R0 on the x-axis
figure_setups;
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
title(['Immunity feedback = ',num2str(immunity_feedback)]);
% plot baseline
P.betaM = 0.25;
[~,~,D,A,~,~,~] = steady_state('EE','fsolve');
R0_baseline = R0_cal();
h3 = plot([R0_baseline,R0_baseline],[0,1],'m-');
legend([h1 h4 h5 h2 h3], {'$D_H+A_H$','$D_H$','$A_H$','Unstable','Baseline'},'Location','nw')
% legend([h1 h4 h5 h2 h3], {'$D_H+A_H$','$D_H$','$A_H$','Unstable','Baseline'},'Location','e')
axis([0 max(R0_list) 0 1])
% title(['Bifurcation']);
xlabel('$\mathcal{R}_0 (\beta_M)$');
%keyboard
% save data for comparison
save(['Results/Bifur/bifur_immune_',num2str(immunity_feedback),'.mat'],'I_frac','ind_stable','ind_unstable','R0_DFE_EE','R0_baseline','R0_list','D_frac_EE','A_frac_EE')
%keyboard
%% Plot with betaM on the x-axis
figure_setups;
beta_list = param; % for betaM on the x-axis
hold on;
beta_DFE_EE = [beta_list,beta_list];
ind_stable = find([re_max_DFE,re_max_EE]<=0);
ind_unstable = find([re_max_DFE,re_max_EE]>=0);
I_frac = [I_frac_DFE,I_frac_EE];
h1 = plot(beta_DFE_EE(ind_stable), I_frac(ind_stable),'-','Marker','.','MarkerSize',30);
h4 = plot(beta_list, D_frac_EE,'-','Marker','o','MarkerSize',10);
h5 = plot(beta_list, A_frac_EE,'-','Marker','^','MarkerSize',10);
h2 = plot(beta_DFE_EE(ind_unstable), I_frac(ind_unstable),'r.','Marker','^','MarkerSize',5);
grid on; grid minor
xlabel('$\mathcal{R}_0(\beta_M)$');
ylabel('Fraction of population');
title(['Immunity feedback = ',num2str(immunity_feedback)]);
% plot baseline
P.betaM = 0.25;
[~,~,D,A,~,~,~] = steady_state('EE','fsolve');
h3 = plot([P.betaM,P.betaM],[0,1],'m-');
% legend([h1 h4 h5 h2 h3], {'$D_H+A_H$','$D_H$','$A_H$','Unstable','Baseline'},'Location','nw')
legend([h1 h4 h5 h2 h3], {'$D_H+A_H$','$D_H$','$A_H$','Unstable','Baseline'},'Location','e')
axis([0 max(beta_list) 0 1])
% title(['Bifurcation']);
xlabel('$\beta_M$');
%% compare immunity feedback impact
if comparison == 1
    immunity_feedback = 1;
    D0 = load(['Results/Bifur/bifur_immune_',num2str(immunity_feedback),'.mat'],'I_frac','ind_stable','ind_unstable','R0_DFE_EE','R0_baseline');
    immunity_feedback = 1;
    D1 = load(['Results/Bifur/bifur_immune_',num2str(immunity_feedback),'.mat'],'I_frac','ind_stable','ind_unstable','R0_DFE_EE','R0_baseline');
    figure_setups; hold on;
    h1 = plot(D0.R0_DFE_EE(D0.ind_stable), D0.I_frac(D0.ind_stable),'--','Marker','.','MarkerSize',30);
    plot(D0.R0_DFE_EE(D0.ind_unstable), D0.I_frac(D0.ind_unstable),':','Marker','^','MarkerSize',5);
    h2 = plot(D1.R0_DFE_EE(D1.ind_stable), D1.I_frac(D1.ind_stable),'-','Marker','.','MarkerSize',30);
    plot(D1.R0_DFE_EE(D1.ind_unstable), D1.I_frac(D1.ind_unstable),'.','Marker','^','MarkerSize',5);
    grid on; grid minor
    xlabel('$\mathcal{R}_0(\beta_M)$');
    ylabel('D+A');
    title('immunity feedback comparison');
    % plot baseline
    h3 = plot([D0.R0_baseline,D0.R0_baseline],[0,1],'m--');
    plot([D1.R0_baseline,D1.R0_baseline],[0,1],'m-');
    legend([h1 h2 h3], {'immunity off','immunity on','baseline'},'Location','e')
    axis([0 7 0 1])
end