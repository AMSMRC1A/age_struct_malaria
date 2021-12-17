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
beta_list = [0.008, 0.03, 0.25];
for ibeta = 1:3
    P.betaM = beta_list(ibeta);
    Malaria_parameters_transform;
    %% initial condition 'init' 'EE'
    [SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
    %% time evolution
    [SH, EH, DH, AH, SM, EM, IM, Cm, Cac, Ctot] = age_structured_Malaria(da,na,tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
    PH = SH+EH+DH+AH;
    PH_final = PH(:,end); % total human at age a, t = n
    NH = trapz(PH,1)*da;
    NM = SM+EM+IM;
    
    %% Age proportions at tfinal prop
    figure_setups_3;
    plot(a/365,SH(:,end)./PH_final,'-','Color',colour_mat1); hold on;
    plot(a/365,EH(:,end)./PH_final,':','Color',colour_mat4);
    plot(a/365,AH(:,end)./PH_final,'--','Color',colour_mat3);
    plot(a/365,DH(:,end)./PH_final,'-.','Color',colour_mat2);
    plot(a/365,PH_final./PH_final,'-k');
    title('Final age dist. (dynamic)'); 
    xlabel('Age (years)');
    ylabel('Fraction of population')
    grid on
    axis([0 P.age_max/365 0 1.1]);
    xlim([0 30])
    text(10.247376311844128,0.915342082307607,['$\beta_M=$',num2str(beta_list(ibeta),3)],'EdgeColor','k','FontSize',33,'LineWidth',2)
    set(gcf, 'Renderer', 'Painters');
    if ibeta == 1
        legend('$\widetilde{S}_H$','$\widetilde{E}_H$','$\widetilde{A}_H$', '$\widetilde{D}_H$','$\widetilde{P}_H$','location','e');   
    end
    print(['Figures/solu_final_age_prop_beta_',num2str(beta_list(ibeta),4),'.eps'], '-depsc','-r300')
end