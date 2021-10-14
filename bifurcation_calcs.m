clear all
% close all
% clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

%% numerical config
tfinal = 100*365; age_max = 60*365; P.age_max = age_max;
dt = 100; % time/age step size in days, default = 50; could go dt = 200 (still robust)
da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);

P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

Malaria_parameters_baseline;
lP = 'betaM'; % bifurcating parameters
immunity_feedback = 1;
if immunity_feedback == 0
    % average populational sigmoids f0 = f1 = average
    
    P.rho_f_0 = 0.109; % value at zero
    P.rho_f_1 = 0.109; % value at L (function saturates to this value)
    
    P.phi_f_0 = 0.386; % value at zero
    P.phi_f_1 = 0.386; % value at L (function saturates to this value)
    
    P.psi_f_0 = 0.579; % value at zero
    P.psi_f_1 = 0.579; % value at L (function saturates to this value)
    param = [0.01:0.05:1.0].^2;
else
    param = [0.01:0.025:0.55].^2;
end
%% options
plot_equilibrium = 0; % can set to zero if working with the DFE

%% Check if stable numerically
re_max_EE = NaN(1,length(param));
I_frac_EE = NaN(1,length(param));
D_EE = NaN(1,length(param));
A_EE = NaN(1,length(param));
re_max_DFE = NaN(1,length(param));
I_frac_DFE = NaN(1,length(param));
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
        FileName = ['Results/Bifur/EE_',num2str(param(i),'%2.4f'),'.mat'];
        if exist(FileName,'file') % Q_val already calculated before
            S = load(FileName,'x_EE','ee');
            x_EE = S.x_EE;
            ee = S.ee;
        else
            [S,E,D,A,Cac,~,~] = steady_state('EE');
            x0 = [S./P.PH_stable;E./P.PH_stable;D./P.PH_stable;A./P.PH_stable;Cac./P.PH_stable];
            [xsol,err,~,~,jacobian] = fsolve(F_prop,x0,options);
            x_EE = reshape(xsol,[P.na,5]);
            jacobian([1,P.na+1,2*P.na+1,3*P.na+1, 4*P.na+1],:)=0; % zero out the rows
            jacobian(:,[1,P.na+1,2*P.na+1,3*P.na+1, 4*P.na+1])=0; % zero out the columns
            ee = eig(jacobian);
            ee(ee==0)=[];
            %         save(FileName,'x_EE','ee')
        end
        if norm(F_prop(x_EE),Inf) > 10^-5
            disp('EE not achieved')
            keyboard
        end
        if plot_equilibrium == 1
            % plot proportion
            figure_setups; hold on;
            plot(a,x_EE(:,1),'-','Color',colour_mat1);
            plot(a,x_EE(:,2),'-','Color',colour_mat3);
            plot(a,x_EE(:,3),'-','Color',colour_mat2);
            plot(a,x_EE(:,4),'-','Color',colour_mat7);
            plot(a,sum(x_EE,2),'-k');
            legend('SH (solver)','EH (solver)','DH (solver)', 'AH (solver)', 'PH (solver)');
            title('Final Age Dist. prop');
            xlabel('age');
            axis_years(gca,age_max); % change to x-axis to years if needed
            grid on
            axis([0 age_max 0 max(sum(x_EE,2))]);
            keyboard
        end
        I_frac_EE(i) = 1-da*trapz((x_EE(:,1)+x_EE(:,2)).*P.PH_stable);
        D_EE(i) = da*trapz(x_EE(:,3).*P.PH_stable);
        A_EE(i) = da*trapz(x_EE(:,4).*P.PH_stable);
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
    FileName = ['Results/Bifur/DFE_',num2str(param(i),'%2.4f'),'.mat'];
    if exist(FileName,'file') % Q_val already calculated before
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
        %  save(FileName,'x_DFE','ee')
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
figure_setups;
hold on;
R0_DFE_EE = [R0_list,R0_list];
ind_stable = find([re_max_DFE,re_max_EE]<0);
ind_unstable = find([re_max_DFE,re_max_EE]>0);
I_frac = [I_frac_DFE,I_frac_EE];
h1 = plot(R0_DFE_EE(ind_stable), I_frac(ind_stable),'-','Marker','.','MarkerSize',30);
h4 = plot(R0_list, D_EE,'-','Marker','.','MarkerSize',30);
h5 = plot(R0_list, A_EE,'-','Marker','.','MarkerSize',30);
h2 = plot(R0_DFE_EE(ind_unstable), I_frac(ind_unstable),'r.','Marker','^','MarkerSize',5);
grid on; grid minor
xlabel('$\mathcal{R}_0(\beta_M)$');
ylabel('Fraction of population');
title(['immunity feedback = ',num2str(immunity_feedback)]);
% plot baseline
P.betaM = 0.25;
[~,~,D,A,~,~,~] = steady_state('EE');
R0_baseline = R0_cal();
h3 = plot([R0_baseline,R0_baseline],[0,1],'m-');
legend([h1 h4 h5 h2 h3], {'D+A','D','A','unstable','baseline'},'Location','e')
axis([0 7 0 1])
keyboard
% save data for comparison
save(['Results/Bifur/bifur_immune_',num2str(immunity_feedback),'.mat'],'I_frac','ind_stable','ind_unstable','R0_DFE_EE','R0_baseline','R0_list','D_EE','A_EE')
keyboard
%% compare immunity feedback impact
immunity_feedback = 0;
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