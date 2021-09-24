clear all
% close all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

%% numerical config
tfinal = 100*365; % final time in days
age_max = 60*365; % max ages in days
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

% model parameters - rates are in 1/day
Malaria_parameters_baseline;

%% options 
error_tolerance = 1e-20;
options = optimoptions('fsolve','Display','none','OptimalityTolerance', error_tolerance);
plot_equilibrium = 1; % can set to zero if working with the DFE
param = 0.125;%linspace(0.1125,0.25,3);
%% Check if DFE is stable numerically
re_max_EE = zeros(1,length(param));
sol_norm_endemic = zeros(1,length(param));
error_EE = zeros(1,length(param));
re_max_DFE = zeros(1,length(param));
sol_norm_DFE = zeros(1,length(param));
error_DFE = zeros(1,length(param));
for i = 1:length(param)
    disp(['progress = ',num2str((i-1)/(length(param))*100),'%']);
    
    P.betaM = param(i);
    
    IC = @(x) [-x(1)+1;x(P.na+1);x(P.na*2+1);x(P.na*3+1)];
    RHS = @(x) human_model_der_fun(x);
    F = @(x) [IC(x);RHS(x)];
    %% solve for EE
    x0 = [0.4*ones(length(a),1); 0.2*ones(length(a),1); 0.2*ones(length(a),1); 0.2*ones(length(a),1)]; % initial guess for the EE
    [xsol,fval,~,~,jacobian] = fsolve(F,x0,options); 
    x_EE = reshape(xsol,[P.na,4]);
    sol_norm_endemic(i) = da*trapz(x_EE(:,1).*P.PH_stable);
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
    error_EE(i) = max(fval);
    re_max_EE(i) = max(real(eig(jacobian)));
    if re_max_EE(i) < 0
       disp(['max real part of eigenvalues = ',num2str(re_max_EE(i),'%10.6f'), '; EE is stable']);
    else
       disp(['max real part of eigenvalues = ',num2str(re_max_EE(i),'%10.6f'), '; EE is unstable']);
    end
    %% Solve for DFE
    x0 = [ones(length(a),1); 0*ones(length(a),1); 0*ones(length(a),1); 0*ones(length(a),1)]; % initial guess for the DFE
    [xsol,fval,~,~,jacobian] = fsolve(F,x0,options);
    x_DFE = reshape(xsol,[P.na,4]);
    sol_norm_DFE(i) = da*trapz(x_DFE(:,1)); % convert to years
    error_DFE(i) = max(fval);
    re_max_DFE(i) = max(real(eig(jacobian)));
    if re_max_DFE(i) < 0
       disp(['max real part of eigenvalues = ',num2str(re_max_DFE(i),'%10.6f'), '; DFE is stable']);
    else
       disp(['max real part of eigenvalues = ',num2str(re_max_DFE(i),'%10.6f'), '; DFE is unstable']);
    end
end
keyboard
%figure_setups;
% scatter(real(eig(jacobian)),imag(eig(jacobian)));
% grid on;
% xline(0,'r','LineWidth',2);
% title('Linearized spectrum at DFE');
% xlim([min(temp_eig) re_max+0.1]);
%% Plot the results
figure_setups;
hold on;
% endemic equilibrium plotting
for j = 1:length(param)
    if re_max_EE(j) < 0 % stable
        scatter(param(j),sol_norm_endemic(j),'b','filled'); 
    elseif re_max_EE(j) > 0 % unstable
        scatter(param(j),sol_norm_endemic(j),'b');
    end
end
% DFE plotting
for j = 1:length(param)
    if re_max_DFE(j) < 0 % stable
        scatter(param(j),sol_norm_DFE(j),'b','filled'); % stable
    elseif re_max_DFE(j) > 0  % unstable
        scatter(param(j),sol_norm_DFE(j),'b');
    end
end
grid on;
xlabel('$\beta_M$');
ylabel('$|| \tilde{S}_H ||$');
set(get(gca,'ylabel'),'rotation',0)
%%
toc