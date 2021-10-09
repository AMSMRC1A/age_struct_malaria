clear all
clc
close all
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

tic
%% numerical config
P.balance_fertility = 1; % 0 to keep original fertility, 1 for balanced birth rate so that pop. growth is zero
P.solution_plots = 0; % toggle solution plots
P.sigmoid_surfs = 0; % toggle surfaces from the sigmoids
P.age_aEIR_surf = 1;

tfinal = 200*365; % final time in days
age_max = 60*365; % max ages in days
P.age_max = age_max;
dt = 10; % time/age step size in days
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

L = [25]; % 25 seems most interesting, gives a threshold around age 5 or so
EIR_var = 'betaM'; % use this parameter to adjust the EIR
var_list = [0.05:0.01:0.07, 0.08, 0.1];
final_immunity = zeros(na,length(var_list));
final_EIR = zeros(1,length(var_list));
for ii = 1:length(L)
    % model parameters - rates are in 1/day
    Malaria_parameters_baseline;
    P.L = L(ii);
    for jj = 1:length(var_list)
        P.(EIR_var) = var_list(jj);
        Malaria_parameters_transform;
        %% This section is one simulation
        [SH, EH, DH, AH, SM, EM, IM, Cm, Cac, Ctot] = age_structured_Malaria(na,da,nt);
        PH = SH+EH+DH+AH;
        PH_final = PH(:,end); % total human at age a, t = n
        NH = trapz(PH,1)*da;
        NM = SM+EM+IM;
        %% Record data for the age vs aEIR surface to mirror the R-B paper
        if P.age_aEIR_surf == 1
            final_immunity(:,jj) = Ctot(:,end)./PH_final;
            [bH,~] = biting_rate(NH(end),NM);
            lamH = FOI_H(bH,IM(1,end),NM);
            final_EIR(1,jj) = lamH/P.betaM*365;
            disp(['EIR = ',num2str(final_EIR(1,jj),'%10.6f')]);
            R0 = R0_cal()
        end
        %% Plot the sigmoids
        if P.sigmoid_surfs == 1
            figure_setups;
            PH = SH+EH+DH+AH;
            subplot(2,2,1), imagesc(t/365,a/365,sigmoid_prob(Ctot./PH, 'phi'),[0 1]); colorbar; title('$\phi(\tilde{C}_{tot})$'); xlabel('time'); ylabel('age');
            set(gca,'YDir','normal');
            grid on
            subplot(2,2,2), imagesc(t/365,a/365,sigmoid_prob(Ctot./PH, 'rho'),[0 1]); colorbar; title('$\rho(\tilde{C}_{tot})$'); xlabel('time'); ylabel('age');
            set(gca,'YDir','normal');
            grid on
            subplot(2,2,3), imagesc(t/365,a/365,sigmoid_prob(Ctot./PH, 'psi'),[0 1]); colorbar; title('$\psi(\tilde{C}_{tot})$'); xlabel('time'); ylabel('age');
            set(gca,'YDir','normal');
            grid on
            % Immunity distribution
            PH = SH+EH+DH+AH;
            subplot(2,2,4), plot(a/365,Ctot(:,end)./PH(:,end),'-.');
            xlabel('age')
            title('$\tilde{C}_{total}(t)$');
            axis_years(gca,age_max); % change to x-axis to years if needed
            grid on
        end
        %% Plot the results
        if P.solution_plots == 1
            PH_final = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end); % total human at age a, t = n
            NH(end) = trapz(PH_final)*da;
            % Population size versus time
            figure_setups;
            Nh = (trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1))*da;
            subplot(2,2,1), plot(t,trapz(SH,1)*da,'-','Color',colour_mat1); hold on;
            subplot(2,2,1), plot(t,trapz(EH,1)*da,'--','Color',colour_mat3);
            subplot(2,2,1), plot(t,trapz(AH,1)*da,'-.','Color',colour_mat2);
            subplot(2,2,1), plot(t,trapz(DH,1)*da,'-','Color',colour_mat7);
            subplot(2,2,1), plot(t,Nh,'-.k')
            legend('$S_H$','$E_H$','$A_H$', '$D_H$','$N_H$');
            title('Population size');
            xlabel('time');
            axis_years(gca,tfinal); % change to x-axis to years if needed
            grid on
            axis([0 tfinal 0 max(Nh)+0.1]);
            % Age profiles at tfinal
            subplot(2,2,2), plot(a,SH(:,end),'-','Color',colour_mat1); hold on;
            subplot(2,2,2), plot(a,EH(:,end),'--','Color',colour_mat3);
            subplot(2,2,2), plot(a,AH(:,end),'-.','Color',colour_mat2);
            subplot(2,2,2), plot(a,DH(:,end),'-','Color',colour_mat7);
            subplot(2,2,2), plot(a,(SH(:,end)+AH(:,end)+EH(:,end)+DH(:,end)),'-.k');
            xlabel('age');
            title('Final Age Dist.');
            axis_years(gca,age_max); % change to x-axis to years if needed
            grid on
            % Age props at tfinal
            subplot(2,2,3), plot(a,SH(:,end)./PH_final,'-','Color',colour_mat1); hold on;
            subplot(2,2,3), plot(a,EH(:,end)./PH_final,'--','Color',colour_mat3);
            subplot(2,2,3), plot(a,AH(:,end)./PH_final,'-.','Color',colour_mat2);
            subplot(2,2,3), plot(a,DH(:,end)./PH_final,'-','Color',colour_mat7);
            subplot(2,2,3), plot(a,(SH(:,end)+AH(:,end)+EH(:,end)+DH(:,end))./PH_final,'-.k');
            title('Final Age Props.');
            xlabel('age');
            axis_years(gca,age_max); % change to x-axis to years if needed
            grid on
        end
    end
end
%% Plot the age vs aEIR surface to mirror the R-B paper
if P.age_aEIR_surf == 1
    figure_setups;
    imagesc(a/365,final_EIR,final_immunity'); 
    xlim([0 10])
    xlabel('age')
    ylabel('EIR')
    title('Immunity levels');
    set(gca,'YDir','normal');
    grid on
    colormap jet
    colorbar;
end
%%
toc