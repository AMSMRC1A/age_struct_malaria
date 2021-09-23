clear all
% clc
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
options = optimoptions('fsolve','Display','none','OptimalityTolerance',...
    error_tolerance);
plot_equilibrium = 1; % can set to zero if working with the DFE
param = P.betaM; %0.1125:0.025:0.25;
%% Check if DFE is stable numerically
re_max_endemic = zeros(1,length(param));
sol_norm_endemic = zeros(1,length(param));
error_endemic = zeros(1,length(param));
re_max_DFE = zeros(1,length(param));
sol_norm_DFE = zeros(1,length(param));
error_DFE = zeros(1,length(param));
for i = 1:length(param)
    disp(['progress = ',num2str((i-1)/(length(param))*100),'%']);
    IC = @(x) [-x(1,1)+1;x(1,2);x(1,3);x(1,4)];
    RHS = @(x) human_model_der_fun(x);
    F = @(x) [IC(x);RHS(x)];
    % solve for the endemic equilibrium 
    x0 = [0.4*P.PH_stable, 0.2*P.PH_stable, 0.2*P.PH_stable, 0.2*P.PH_stable]; % initial guess for the DFE
    [x_EE,fval,~,~,jacobian] = fsolve(F,x0,options);
    sol_norm_endemic(i) = da*trapz(x_EE(:,1))/365; % convert to years
    if plot_equilibrium == 1
        % plot the solution of the nonlinear solver
        figure_setups; hold on;
        plot(a,x_EE(:,1),'-.','Color',colour_mat1); 
        plot(a,x_EE(:,2),'-.','Color',colour_mat3);
        plot(a,x_EE(:,4),'-.','Color',colour_mat2);
        plot(a,x_EE(:,3),'-.','Color',colour_mat7);
%         plot(a,sum(x_EE,2),'-.k');
        legend('SH (solver)','EH (solver)','AH (solver)', 'DH (solver)');
        title('Final Age Dist. Proportions');
        xlabel('age');
        axis_years(gca,age_max); % change to x-axis to years if needed
        grid on
        axis([0 age_max 0 max(sum(x_EE,2))*1.1]);
    end
    keyboard
    error_endemic(i) = max(max(fval));
    re_max_endemic(i) = max(real(eig(jacobian)));
    if re_max_endemic(i) < 0
       disp(['max real part of eigenvalues = ',num2str(re_max_endemic(i),'%10.6f'), '; endemic equilibrium is stable']);
    else
       disp(['max real part of eigenvalues = ',num2str(re_max_endemic(i),'%10.6f'), '; endemic equilibrium is unstable']);
    end
    keyboard
    % analyze the DFE
    P2 = [ones(length(a),1), 0*ones(length(a),1), 0*ones(length(a),1), 0*ones(length(a),1)]; % initial guess for the DFE
    [x_EE,fval,~,~,jacobian] = fsolve(F_PSI,P2,options);
    sol_norm_DFE(i) = da*trapz(x_EE(:,1))/365; % convert to years
    error_DFE(i) = max(max(fval));
    temp_eig = sort(real(eig(jacobian)),'descend');
    re_max_DFE(i) = max(temp_eig);
end
%figure_setups;
% scatter(real(eig(jacobian)),imag(eig(jacobian)));
% grid on;
% xline(0,'r','LineWidth',2);
% title('Linearized spectrum at DFE');
% xlim([min(temp_eig) re_max+0.1]);
%% Plot the results
%figure_setups;
figure(9);
hold on;
error_tolerance = 1e-6;
% endemic equilibrium plotting
for j = 1:length(param)
    if re_max_endemic(j) < 0 && error_endemic(j) < error_tolerance
        scatter(param(j),sol_norm_endemic(j),'b','filled');
    elseif re_max_endemic(j) < 0 && error_endemic(j) > error_tolerance
        scatter(param(j),sol_norm_endemic(j),'r','filled');
    elseif re_max_endemic(j) > 0 && error_endemic(j) > error_tolerance
        scatter(param(j),sol_norm_endemic(j),'r');
    else
        scatter(param(j),sol_norm_endemic(j),'b');
    end
end
% DFE plotting
for j = 1:length(param)
    if re_max_DFE(j) < 0 && error_DFE(j) < error_tolerance
        scatter(param(j),sol_norm_endemic(j),'b','filled');
    elseif re_max_DFE(j) < 0 && error_DFE(j) > error_tolerance
        scatter(param(j),sol_norm_DFE(j),'r','filled');
    elseif re_max_DFE(j) > 0 && error_DFE(j) > error_tolerance
        scatter(param(j),sol_norm_DFE(j),'r');
    else
        scatter(param(j),sol_norm_DFE(j),'b');
    end
end
grid on;
xlabel('$\beta_M$');
ylabel('$|| \tilde{S}_H ||$');
set(get(gca,'ylabel'),'rotation',0)
%%
toc