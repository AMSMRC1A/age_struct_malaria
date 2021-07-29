clear all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2
tic
%% setup parameters, stable age dist, etc.
P.balance_fertility = 1; % 0 to keep original fertility, 1 for balanced birth rate so that pop. growth is zero
age_max = 60*365; % max ages in days
P.age_max = age_max;
da = 20; % age step in days
a = (0:da:age_max)';
na = length(a);

P.a = a;
P.na = na;
P.da = da;

% model parameters - rates are in 1/day
baseline_Malaria_parameters;
find_stable_age;
param = 0.2:0.1:1; % changing P.bm the mosquito desired bites
%% Solve for and check the stability of the endemic equilibrium

%% Check if DFE is stable numerically
% first row is different to account for boundary conditions

error_tolerance = 1e-20;
options = optimoptions('fsolve','Display','none','OptimalityTolerance',...
    error_tolerance);
re_max = zeros(1,length(param));
sol_norm = zeros(1,length(param));
error = zeros(1,length(param));
for i = 1:length(param)
    P.bm = param(i);
    P.phi = 1/2;
    P.psi = 1/2;
    P.rho = 1/2;
    [bH,bM] = biting_rate(1,P.gM/P.muM);
    Lambda_M = @(x) bM*P.Lambda*da*trapz(exp(-P.muH_int).*(P.betaD*x(:,3)+P.betaA*x(:,4)));
    Lambda_H = @(x) bH*P.betaM*(P.sigma/(P.sigma+P.muM))*((Lambda_M(x))/(Lambda_M(x) + P.muM));
    F_PSI = @(x) [[-x(1,1)+1; -Lambda_H(x)*x(2:end,1) + P.phi*P.rD*x(2:end,3) + P.rA*x(2:end,4) - diff(x(:,1))/da];...
    [-x(1,2); Lambda_H(x)*x(2:end,1) - P.h*x(2:end,2) - diff(x(:,2))/da];...
    [-x(1,3); P.rho*P.h*x(2:end,2) + P.psi*Lambda_H(x)*x(2:end,4) - P.rD*x(2:end,3) - diff(x(:,3))/da];...
    [-x(1,4); (1-P.rho)*P.h*x(2:end,2) - P.psi*Lambda_H(x)*x(2:end,4) + (1-P.phi)*P.rD*x(2:end,3) - P.rA*x(2:end,4) - diff(x(:,4))/da]];
    P1 = [ones(length(a),1), 0*ones(length(a),1), 0*ones(length(a),1), 0*ones(length(a),1)]; % initial guess for the DFE
    [eq_age,fval,exitflag,output,jacobian] = fsolve(F_PSI,P1,options);
    sol_norm(i) = da*trapz(eq_age(:,1));
    error(i) = max(max(fval));
    %     figure_setups;
    %     plot(a,eq_age(:,1),'-','Color',colour_mat1); hold on;
    %     plot(a,eq_age(:,2),'-','Color',colour_mat3);
    %     plot(a,eq_age(:,4),'-','Color',colour_mat2);
    %     plot(a,eq_age(:,3),'-','Color',colour_mat7);
    %     plot(a,sum(eq_age,2),'-.k');
    %     legend('SH (DFE)','EH (DFE)','AH (DFE)', 'DH (DFE)','PH (DFE)');
    %     title('Final Age Dist. Proportions (DFE)');
    % get the eigenvalues of the Jacobian and plot the linearized spectrum
    temp_eig = sort(real(eig(jacobian)),'descend');
    re_max(i) = max(temp_eig);
    if re_max < 0
        disp(['max real part of eigenvalues at DFE = ',num2str(re_max(i),'%10.6f'), '; DFE is stable']);
    else
        disp(['max real part of eigenvalues = ',num2str(re_max(i),'%10.6f'), '; DFE is unstable']);
    end
end
%figure_setups;
% scatter(real(eig(jacobian)),imag(eig(jacobian)));
% grid on;
% xline(0,'r','LineWidth',2);
% title('Linearized spectrum at DFE');
% xlim([min(temp_eig) re_max+0.1]);
%% Plot the results
figure_setups;
hold on;
error_tolerance = 1e-6;
for j = 1:length(param)
    if re_max(j) < 0 && error(j) < error_tolerance
        scatter(param(j),sol_norm(j),'b','filled');
    elseif re_max(j) < 0 && error(j) > error_tolerance
        scatter(param(j),sol_norm(j),'r','filled');
    elseif re_max(j) > 0 && error(j) > error_tolerance
        scatter(param(j),sol_norm(j),'r');
    else
        scatter(param(j),sol_norm(j),'b');
    end
end
grid on;
xlabel('$b_m$');
ylabel('$|| \tilde{S}_H ||$');
set(get(gca,'ylabel'),'rotation',0)
%%
toc