%% plot Fig4 row 2 - Population proportion at three scenarios

close all
clear all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2
%% numerical config
tfinal = 100*365; % final time in days
age_max = 80*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

Malaria_parameters_baseline; % model parameters - rates are in 1/day
lP = 'betaM';  % bifurcating parameters
immunity_feedback = 0; % -1 = low immunity fixed, 1 = dynamic, 0 = high fixed
if immunity_feedback == -1 
    P.phi_f_0 = 0.230971153687268; % value at zero
    P.phi_f_1 = 0.230971153687268; % value at L (function saturates to this value)
    P.rho_f_0 = 0.907631179204690; % value at zero
    P.rho_f_1 = 0.907631179204690 ; % value at L (function saturates to this value)
    P.psi_f_0 = 0.907631179204690 ; % value at zero
    P.psi_f_1 = 0.907631179204690 ; % value at L (function saturates to this value)
    P.(lP) = 0.088276535316648; % gives R0=4;
elseif immunity_feedback == 0 % betaM = 0.25
    % average populational sigmoids f0 = f1 = average
    P.phi_f_0 = 0.915792480087329; % value at zero
    P.phi_f_1 = 0.915792480087329; % value at L (function saturates to this value)
    P.rho_f_0 = 0.114825053290306; % value at zero
    P.rho_f_1 = 0.114825053290306; % value at L (function saturates to this value
    P.psi_f_0 = 0.114825053290306; % value at zero
    P.psi_f_1 = 0.114825053290306; % value at L (function saturates to this value)
    P.(lP) = 0.344878578683840;  % gives R0=4;
elseif immunity_feedback == 1
    P.(lP) = 0.083296663925047; % gives R0=4;
end

Malaria_parameters_transform;
R0_cal()
%% initial condition 'init' 'EE'
[SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
%% time evolution
[SH, EH, DH, AH, SM, EM, IM, Cm, Cac, Ctot] = age_structured_Malaria(da,na,tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
PH = SH+EH+DH+AH;
PH_final = PH(:,end); % total human at age a, t = n
NH = trapz(PH,1)*da;
NM = SM+EM+IM;

%% EIR
[bh,bm] = biting_rate(NH,NM);
EIR = bh.*IM./NM*365;
EIR_EE = EIR(end)
tic
R0 = R0_cal()
toc
%% Age proportions at tfinal prop
figure_setups_3;
plot(a/365,SH(:,end)./PH_final,'-','Color',colour_mat1); hold on;
plot(a/365,EH(:,end)./PH_final,':','Color',colour_mat4);
plot(a/365,AH(:,end)./PH_final,'--','Color',colour_mat3);
plot(a/365,DH(:,end)./PH_final,'-.','Color',colour_mat2);
plot(a/365,PH_final./PH_final,'-k');
xlabel('Age (years)');
ylabel('Fraction of population')
grid on
axis([0 P.age_max/365 0 1.1]);
xlim([0 30])
text(10.247376311844128,0.915342082307607,'$\mathcal{R}_0=4$','EdgeColor','k','FontSize',33,'LineWidth',2)
set(gcf, 'Renderer', 'Painters');
if immunity_feedback == -1
    title('Final age dist. (fixed low)'); 
    legend('$\widetilde{S}_H$','$\widetilde{E}_H$','$\widetilde{A}_H$', '$\widetilde{D}_H$','$\widetilde{P}_H$','location','e');
    print('Figures/solu_final_age_prop_-1.eps', '-depsc','-r300')
elseif immunity_feedback == 0
    title('Final age dist. (fixed high)');
    print('Figures/solu_final_age_prop_0.eps', '-depsc','-r300')
elseif immunity_feedback == 1
    title('Final age dist. (dynamic)'); 
    print('Figures/solu_final_age_prop_1.eps', '-depsc','-r300')
end
