clear all
% close all
% clc
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

% model parameters - rates are in 1/day
Malaria_parameters_baseline;

%% options 
error_tolerance = 1e-20;
options = optimoptions('fsolve','Display','none','OptimalityTolerance', error_tolerance);
plot_equilibrium = 1; % can set to zero if working with the DFE
param = [0.01];
% param = [0.001, 0.003, 0.005, 0.01, 0.05, 0.1, 0.15, 0.25];
for i = 1:length(param)
P.betaM = param(i);
R0_cal();
end

%% Check if stable numerically
re_max_EE = zeros(1,length(param));
sol_norm_EE = zeros(1,length(param));
re_max_DFE = zeros(1,length(param));
I_frac_DFE = zeros(1,length(param));
for i = 1:length(param)
    disp(['progress = ',num2str((i-1)/(length(param))*100),'%']);
    tic
    P.betaM = param(i);
    F = @(x) human_model_der_prop(x);
    %% solve for EE
    FileName = ['Results/Bifur/EE_',num2str(param(i),5),'.mat'];
    if exist(FileName,'file') % Q_val already calculated before
        S = load(FileName,'x_EE','ee');
        x_EE = S.x_EE;
        ee = S.ee;
    else
        x0 = [1; 0.4*ones(length(a)-1,1); 0; 0.2*ones(length(a)-1,1); 0; 0.2*ones(length(a)-1,1); 0; 0.2*ones(length(a)-1,1)]; % initial guess for the EE
        [xsol,~,~,~,jacobian] = fsolve(F,x0,options);
        x_EE = reshape(xsol,[P.na,4]);
        jacobian([1,P.na+1,2*P.na+1,3*P.na+1],:)=0; % zero out the rows
        jacobian(:,[1,P.na+1,2*P.na+1,3*P.na+1])=0; % zero out the columns
        ee = eig(jacobian);
        keyboard
        ee(ee==0)=[];
        save(FileName,'x_EE','ee')
    end
    if plot_equilibrium == 1
        % plot the solution of the nonlinear solver
        figure_setups; hold on;
        plot(a,x_EE(:,1).*P.PH_stable,'-','Color',colour_mat1); 
        plot(a,x_EE(:,2).*P.PH_stable,'-','Color',colour_mat3);
        plot(a,x_EE(:,3).*P.PH_stable,'-','Color',colour_mat2);
        plot(a,x_EE(:,4).*P.PH_stable,'-','Color',colour_mat7);
        plot(a,sum(x_EE,2).*P.PH_stable,'-k');
        legend('SH (solver)','EH (solver)','DH (solver)', 'AH (solver)', 'PH (solver)');
        title('Final Age Dist.');
        xlabel('age');
        axis_years(gca,age_max); % change to x-axis to years if needed
        grid on
        axis([0 age_max 0 max(sum(x_EE,2).*P.PH_stable)]);
    end
    sol_norm_EE(i) = da*trapz(x_EE(:,1).*P.PH_stable);
    re_max_EE(i) = max(real(ee));
    if re_max_EE(i) < 0
       disp(['max real part of eigenvalues = ',num2str(re_max_EE(i),'%10.6f'), '; EE is stable']);
    else
       disp(['max real part of eigenvalues = ',num2str(re_max_EE(i),'%10.6f'), '; EE is unstable']);
    end
    toc
%     %% Solve for DFE
%     FileName = ['Results/Bifur/DFE_',num2str(param(i),5),'.mat'];
%     if exist(FileName,'file') % Q_val already calculated before
%         S = load(FileName,'x_DFE','ee');
%         x_DFE = S.x_DFE;
%         ee = S.ee;
%         keyboard
%     else
%         x0 = [ones(length(a),1); 0*ones(length(a),1); 0*ones(length(a),1); 0*ones(length(a),1)]; % initial guess for the DFE
%         tic
%         options = optimoptions('fsolve','Display','none');
%         [xsol,err,~,~,jacobian] = fsolve(F,x0,options);
%         max(max(abs(err)))
%         ee = eig(full(jacobian));
%         sort(ee)
%         keyboard
%         jacobian([1, P.na+1, 2*P.na+1, 3*P.na+1],:)=0; % zero out the rows
%         jacobian(:,[1, P.na+1, 2*P.na+1, 3*P.na+1])=0; % zero out the columns        
%         ee(ee==0)=[];
%         toc
%         keyboard


%         F_size = @(x) human_model_der_size(x);
%         x0 = [ones(length(a),1).*P.PH_stable; 0*ones(length(a),1).*P.PH_stable; 0*ones(length(a),1).*P.PH_stable; 0*ones(length(a),1).*P.PH_stable]; % initial guess for the DFE
%         tic
%         options = optimoptions('fsolve','Display','none','MaxIterations',0);
%         [xsol,err,~,~,jacobian] = fsolve(F_size,x0,options);
%         if max(max(abs(err)))>10^-6
%             disp('not achieve steady state')
%             keyboard
%         end
%         ee = eig(full(jacobian));
%         sort(ee)
%         keyboard
%         jacobian([1, P.na+1, 2*P.na+1, 3*P.na+1],:)=0; % zero out the rows
%         jacobian(:,[1, P.na+1, 2*P.na+1, 3*P.na+1])=0; % zero out the columns
%         ee(ee==0)=[];
%         toc
%         save(FileName,'xsol','ee');    
%     end
    re_max_DFE(i) = max(real(ee));
    x_DFE = reshape(xsol,[P.na,4]);
    I_frac_DFE(i) = 1 - da*trapz(x_DFE(:,1)); % convert to years       
    if re_max_DFE(i) < 0
       disp(['max real part of eigenvalues = ',num2str(re_max_DFE(i),'%10.6f'), '; DFE is stable']);
    else
       disp(['max real part of eigenvalues = ',num2str(re_max_DFE(i),'%10.6f'), '; DFE is unstable']);
    end
end
%% Plot the results
R0_list = NaN(size(param));
for i = 1:length(param)
    P.betaM = param(i);
    R0_list(i) = R0_cal();
end
figure_setups;
hold on;
% endemic equilibrium plotting
% for j = 1:length(param)
%     if re_max_EE(j) < 0 % stable
%         scatter(R0_list(j),1-sol_norm_EE(j),'b','filled'); 
%     elseif re_max_EE(j) > 0 % unstable
%         scatter(R0_list(j),1-sol_norm_EE(j),'b');
%     end
% end
% DFE plotting
for j = 1:length(param)
    if re_max_DFE(j) < 0 % stable
        scatter(R0_list(j),I_frac_DFE(j),'b','filled'); % stable
    elseif re_max_DFE(j) > 0  % unstable
        scatter(R0_list(j),I_frac_DFE(j),'b');
    end
end
grid on;
xlabel('$R_0$');
ylabel('Fraction of infected');
%%
toc