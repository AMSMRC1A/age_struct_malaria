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
P.t = t;

% model parameters - rates are in 1/day
Malaria_parameters_baseline;

%% options
plot_equilibrium = 0; % can set to zero if working with the DFE
% bifurcating parameters
lP = 'betaM'; 
param = [0.01:0.01:0.5].^2;
% lP = 'bh';
% param = [sqrt(0.05):sqrt(0.1)/10:sqrt(5)].^2;
%% Check if stable numerically
re_max_EE = NaN(1,length(param));
I_frac_EE = NaN(1,length(param));
re_max_DFE = NaN(1,length(param));
I_frac_DFE = NaN(1,length(param));
for i = 1:length(param)
    disp(['progress = ',num2str((i-1)/(length(param))*100),'%']);
    P.(lP) = param(i);
    Malaria_parameters_transform;
    F_prop = @(x) human_model_der_prop(x);
    options = optimoptions('fsolve','Display','none','MaxIterations',50);
    tic
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
            jacobian([1,P.na+1,2*P.na+1,3*P.na+1],:)=0; % zero out the rows
            jacobian(:,[1,P.na+1,2*P.na+1,3*P.na+1])=0; % zero out the columns
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
        jacobian([1, P.na+1, 2*P.na+1, 3*P.na+1],:) = 0; % zero out the rows
        jacobian(:,[1, P.na+1, 2*P.na+1, 3*P.na+1]) = 0; % zero out the columns
        ee = eig(jacobian);
        ee(ee==0)=[];
        x_DFE = reshape(xsol_DFE,[P.na,5]);
        %         save(FileName,'x_DFE','ee')
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
h1 = plot(R0_DFE_EE(ind_stable), I_frac(ind_stable),'b-','Marker','.','MarkerSize',30);
h2 = plot(R0_DFE_EE(ind_unstable), I_frac(ind_unstable),'r.','Marker','^','MarkerSize',5);
% h1 = scatter(R0_DFE_EE(ind_stable), I_frac(ind_stable),100,'bo','filled');
% h2 = scatter(R0_DFE_EE(ind_unstable), I_frac(ind_unstable),100,'r^','filled');
grid on; grid minor
xlabel('$\mathcal{R}_0(\beta_M)$');
ylabel('D+A');
axis([0 max(R0_DFE_EE) 0 1])
legend([h1 h2], {'stable','unstable'},'Location','e')