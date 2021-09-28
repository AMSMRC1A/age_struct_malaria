clear all
close all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

%% numerical config
tfinal = 100*365; % final time in days
age_max = 60*365; % max ages in days
P.age_max = age_max;
dt = 50; % time/age step size in days, default = 50; could go dt = 200 (still robust)
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

% model parameters - rates are in 1/day
Malaria_parameters_baseline;

%% options 
plot_equilibrium = 0; % can set to zero if working with the DFE
param = [0.01:0.01:0.5].^2; % full range
% param = [0.1:0.01:0.5].^2; % R0>1
% param = [0.01:0.01:0.09].^2; % R0<1
%% Check if stable numerically
re_max_EE = NaN(1,length(param));
I_frac_EE = NaN(1,length(param));
re_max_DFE = NaN(1,length(param));
I_frac_DFE = NaN(1,length(param));
for i = 1:length(param)
    disp(['progress = ',num2str((i-1)/(length(param))*100),'%']);
    P.betaM = param(i);
    F_prop = @(x) human_model_der_prop(x);
    %% solve for EE
    tic
    FileName = ['Results/Bifur/EE_',num2str(param(i),'%2.4f'),'.mat'];
    if exist(FileName,'file') % Q_val already calculated before
        S = load(FileName,'x_EE','ee');
        x_EE = S.x_EE;
        ee = S.ee;
    else
        tfinal = 10*365; % run for 10 years to get closer to EE
        t = (0:dt:tfinal)';
        nt = length(t);
        P.nt = nt;
        P.t = t;
        [SH, EH, DH, AH, ~, ~, ~, ~, ~, ~] = age_structured_Malaria();
%         figure_setups;
%         Nh = (trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1))*da;
%         plot(t,trapz(SH,1)*da./Nh,'-','Color',colour_mat1); hold on;
%         plot(t,trapz(EH,1)*da./Nh,'--','Color',colour_mat3);
%         plot(t,trapz(AH,1)*da./Nh,'-.','Color',colour_mat2);
%         plot(t,trapz(DH,1)*da./Nh,'-','Color',colour_mat7);
%         plot(t,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1))*da./Nh,'-.k');
%         legend('SH-age','EH-age','AH-age', 'DH-age','$N_H$');
%         title('Population proportions vs time');
%         axis_years(gca,tfinal); % change to x-axis to years if needed
%         xlabel('time');
%         grid on
%         axis([0 tfinal 0 1.1]);
        x0 = [SH(:,end)./P.PH_stable;EH(:,end)./P.PH_stable;DH(:,end)./P.PH_stable;AH(:,end)./P.PH_stable];
        error_tolerance = 1e-25;
        options = optimoptions('fsolve','Display','none','OptimalityTolerance', error_tolerance);
        [xsol,err,~,~,jacobian] = fsolve(F_prop,x0,options);
        if max(max(abs(err)))>10^-5
            disp('not converged')
            keyboard
        end
        x_EE = reshape(xsol,[P.na,4]);
        jacobian([1,P.na+1,2*P.na+1,3*P.na+1],:)=0; % zero out the rows
        jacobian(:,[1,P.na+1,2*P.na+1,3*P.na+1])=0; % zero out the columns
        ee = eig(jacobian);
        ee(ee==0)=[];
        save(FileName,'x_EE','ee')
    end
    if norm(F_prop(x_EE),Inf) > 10^-5
        disp('EE not achieved')
        keyboard
    end
    if plot_equilibrium == 1
%         figure_setups; hold on;
%         plot(a,x_EE(:,1).*P.PH_stable,'-','Color',colour_mat1); 
%         plot(a,x_EE(:,2).*P.PH_stable,'-','Color',colour_mat3);
%         plot(a,x_EE(:,3).*P.PH_stable,'-','Color',colour_mat2);
%         plot(a,x_EE(:,4).*P.PH_stable,'-','Color',colour_mat7);
%         plot(a,sum(x_EE,2).*P.PH_stable,'-k');
%         legend('SH (solver)','EH (solver)','DH (solver)', 'AH (solver)', 'PH (solver)');
%         title('Final Age Dist.');
%         xlabel('age');
%         axis_years(gca,age_max); % change to x-axis to years if needed
%         grid on
%         axis([0 age_max 0 max(sum(x_EE,2).*P.PH_stable)]);
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
    I_frac_EE(i) = 1-da*trapz(x_EE(:,1).*P.PH_stable);
    re_max_EE(i) = max(real(ee));    
%     if re_max_EE(i) < 0
%        disp(['max real part of eigenvalues = ',num2str(re_max_EE(i),'%10.6f'), '; EE is stable']);
%     else
%        disp(['max real part of eigenvalues = ',num2str(re_max_EE(i),'%10.6f'), '; EE is unstable']);
%     end
    toc
    %% Solve for DFE
    tic
    FileName = ['Results/Bifur/DFE_',num2str(param(i),'%2.4f'),'.mat'];
    if exist(FileName,'file') % Q_val already calculated before
        S = load(FileName,'x_DFE','ee');
        x_DFE = S.x_DFE;
        ee = S.ee;
    else
        x0 = [ones(length(a),1); 0*ones(length(a),1); 0*ones(length(a),1); 0*ones(length(a),1)]; % initial guess for the DFE       
        options = optimoptions('fsolve','Display','none','MaxIterations',0);
        [xsol_DFE,err,~,~,jacobian] = fsolve(F_prop,x0,options);
        if max(max(abs(err)))>10^-5
            disp('not converged')
            keyboard
        end
        jacobian([1, P.na+1, 2*P.na+1, 3*P.na+1],:) = 0; % zero out the rows
        jacobian(:,[1, P.na+1, 2*P.na+1, 3*P.na+1]) = 0; % zero out the columns        
        ee = eig(jacobian);
        ee(ee==0)=[];
        x_DFE = reshape(xsol_DFE,[P.na,4]);
        save(FileName,'x_DFE','ee')
    end
    if norm(F_prop(x_DFE),Inf) > 10^-5
        disp('DFE not achieved')
        keyboard
    end
    I_frac_DFE(i) = 1 - da*trapz(x_DFE(:,1).*P.PH_stable);
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
    P.betaM = param(i);
    R0_list(i) = R0_cal();
end
figure_setups;
hold on;
plot(R0_list(re_max_EE < 0), I_frac_EE(re_max_EE < 0),'b-','Marker','.','MarkerSize',30);
plot(R0_list(re_max_EE > 0), I_frac_EE(re_max_EE > 0),'r.','Marker','^','MarkerSize',5);
h1 = plot(R0_list(re_max_DFE < 0), I_frac_DFE(re_max_DFE < 0),'b-','Marker','.','MarkerSize',30);
h2 = plot(R0_list(re_max_DFE > 0), I_frac_DFE(re_max_DFE > 0),'r.','Marker','^','MarkerSize',5);
grid on; grid minor
xlabel('$\mathcal{R}_0$');
ylabel('Fraction of infected');
axis([0 max(R0_list) 0 1])
legend([h1 h2], {'stable','unstable'},'Location','e')