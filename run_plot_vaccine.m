close all
clear all
% clc
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

DD = load(['Results/Vaccine/vp0_0.mat'],'t','a','vacc','SH_EE','EH_EE','DH_EE','AH_EE','VH_EE','PH_EE','Cm_EE','Cac_EE','Ctot_EE');
t = DD.t; a = DD.a;
vp0_list = [0,0.8];
vacc_mat = NaN(length(t),length(vp0_list)); vacc_mat(:,1) = DD.vacc';
SH_EE_mat = NaN(length(a),length(vp0_list)); EH_EE_mat = SH_EE_mat; DH_EE_mat = SH_EE_mat; AH_EE_mat = SH_EE_mat; VH_EE_mat = SH_EE_mat; PH_EE_mat = SH_EE_mat;
SH_EE_mat(:,1) = DD.SH_EE; EH_EE_mat(:,1) = DD.EH_EE; DH_EE_mat(:,1) = DD.DH_EE; AH_EE_mat(:,1) = DD.AH_EE; VH_EE_mat(:,1) = DD.VH_EE; PH_EE_mat(:,1) = DD.PH_EE;
Ctot_EE_mat = NaN(length(a),length(vp0_list)); Cm_EE_mat = Ctot_EE_mat; Cac_EE_mat = Ctot_EE_mat;
Ctot_EE_mat(:,1) = DD.Ctot_EE; Cm_EE_mat(:,1) = DD.Cm_EE; Cac_EE_mat(:,1) = DD.Cac_EE; 

for iv = 2:length(vp0_list)
    vp0 = vp0_list(iv);
    DD = load(['Results/Vaccine/vp0_',num2str(vp0*100),'.mat'],'t','a','vp','vacc','SH_EE','EH_EE','DH_EE','AH_EE','VH_EE','PH_EE',...
        'Cm_EE','Cac_EE','Ctot_EE');
    vacc_mat(:,iv) = DD.vacc'; SH_EE_mat(:,iv) = DD.SH_EE; EH_EE_mat(:,iv) = DD.EH_EE; DH_EE_mat(:,iv) = DD.DH_EE; AH_EE_mat(:,iv) = DD.AH_EE; VH_EE_mat(:,iv) = DD.VH_EE; PH_EE_mat(:,iv) = DD.PH_EE;  
    Ctot_EE_mat(:,iv) = DD.Ctot_EE; Cm_EE_mat(:,iv) = DD.Cm_EE; Cac_EE_mat(:,iv) = DD.Cac_EE;
    vp = DD.vp;
end

mrk_style = {'none','o','s','*','v','+','^'}; 
%% vaccine #
figure_setups; hold on
grid on
for iv = 1:length(vp0_list)
    plot(t/365,vacc_mat(:,iv),'Marker',mrk_style{iv},'MarkerIndices',1:100:length(t))
end
legendStrings = "vp0 = " + string(vp0_list);
legend(legendStrings)
xlim([0 max(t)/365]);
ind1 = find(vp,1,'first');
ind2 = find(vp,1,'last');
%% Age profiles at tfinal
figure_setups; hold on
for iv = 1:length(vp0_list)
    plot(a/365,DH_EE_mat(:,iv),'-','Color',colour_mat2,'Marker',mrk_style{iv},'MarkerIndices',[1:2:10,20:20:length(a)],'DisplayName',['$D_H$, $\nu_p^0 = ', num2str(vp0_list(iv)),'$']);
    plot(a/365,AH_EE_mat(:,iv),'--','Color',colour_mat3,'Marker',mrk_style{iv},'MarkerIndices',[1:2:10,20:20:length(a)],'DisplayName',['$A_H$, $\nu_p^0 = ', num2str(vp0_list(iv)),'$']);
%     plot(a/365,VH_EE_mat(:,iv),'-','Color',colour_mat5,'Marker',mrk_style{iv},'MarkerIndices',1:20:length(a),'DisplayName',['VH, vp0 = ', num2str(vp0_list(iv))]);
%     plot(a/365,PH_EE_mat(:,iv),'-k','DisplayName',['PH, vp0 = ', num2str(vp0_list(iv))]);
end
legend('AutoUpdate','off','NumColumns',2,'Location','se');
fill([a(ind1)/365 a(ind1)/365 a(ind2)/365 a(ind2)/365],[-1,1,1,-1],colour_mat5,'FaceAlpha',0.6, 'EdgeColor',colour_mat5,'EdgeAlpha', 0.6,...
    'DisplayName','Vacc. age range');
fill([8.4 8.4 9.4 9.4],[2.6*10^-5 3.2*10^-5 3.2*10^-5 2.6*10^-5],'r','FaceColor','none','EdgeColor','k','LineWidth',3);
% title(['Final age distribution']);
ylabel('Population density')
xlabel('Age (years)');
grid on
axis([0 15 0 8*10^-5]);
%% zoom in plot
figure_setups; hold on
for iv = 1:length(vp0_list)
    plot(a/365,DH_EE_mat(:,iv),'-','Color',colour_mat2,'Marker',mrk_style{iv},'MarkerIndices',1:5:length(a),...
        'LineWidth',8,'MarkerSize',20);
    plot(a/365,AH_EE_mat(:,iv),'--','Color',colour_mat3,'Marker',mrk_style{iv},'MarkerIndices',1:5:length(a),...
        'LineWidth',8,'MarkerSize',20);
end
axis([8.4 9.4 2.6*10^-5 3.2*10^-5]) 
xticks('')
yticks('')
grid off
%% Difference
figure_setups; hold on
% fill([a(ind1)/365 a(ind1)/365 a(ind2)/365 a(ind2)/365],[-1,1,1,-1],colour_mat5,'FaceAlpha',0.6, 'EdgeColor',colour_mat5,'EdgeAlpha', 0.6,...
%     'DisplayName','Vacc. age range');
for iv = 2:length(vp0_list)
    plot(a/365,DH_EE_mat(:,1)-DH_EE_mat(:,iv),'-','Color',colour_mat2,...
        'DisplayName',['Reduction in $D_H$']);
    plot(a/365,AH_EE_mat(:,1)-AH_EE_mat(:,iv),'--','Color',colour_mat3,...
        'DisplayName',['Reduction in $A_H$']);
end
legend;
% title(['Difference in distribution']);
ylabel('Population density')
xlabel('Age (years)');
grid on
axis([0 15 -2*10^-6 9*10^-6]);

%% Immunity profiles at tfinal
figure_setups; hold on
for iv = 1:length(vp0_list)
    plot(a/365,Cm_EE_mat(:,iv),'-','Color',colour_mat2,'Marker',mrk_style{iv},'MarkerIndices',[1:2:10,15:20:length(a)],...
        'DisplayName',['$C_m$, $\nu_p^0 = ', num2str(vp0_list(iv)),'$']);
    plot(a/365,Cac_EE_mat(:,iv),'--','Color',colour_mat1,'Marker',mrk_style{iv},'MarkerIndices',[1,3,10:20:length(a)],...
        'DisplayName',['$C_e$, $\nu_p^0 = ', num2str(vp0_list(iv)),'$']);
end
legend('AutoUpdate','off','NumColumns',2,'Location','se');
fill([a(ind1)/365 a(ind1)/365 a(ind2)/365 a(ind2)/365],[-1,1,1,-1],colour_mat5,'FaceAlpha',0.6, 'EdgeColor',colour_mat5,'EdgeAlpha', 0.6,...
    'DisplayName','Vacc. age range');
fill([0.4984 0.4984 1.3591 1.3591],[7.1556e-05 1.4222e-04 1.4222e-04 7.1556e-05],'r','FaceColor','none','EdgeColor','k','LineWidth',3);
% title(['Final Immunity Distribution']);
ylabel('Immunity level')
xlabel('Age (years)');
grid on
axis([0 15 0 8*10^-4]);
% axis([0 2 0.5*10^-4 2*10^-4]);
%% zoom in plot
figure_setups; hold on
for iv = 1:length(vp0_list)
    plot(a/365,Cm_EE_mat(:,iv),'-','Color',colour_mat2,'Marker',mrk_style{iv},'MarkerIndices',[1:15,15:5:length(a)],...
        'LineWidth',10,'MarkerSize',20);
    plot(a/365,Cac_EE_mat(:,iv),'--','Color',colour_mat1,'Marker',mrk_style{iv},'MarkerIndices',[1,3,10:5:length(a)],...
        'LineWidth',10,'MarkerSize',20);
end
fill([a(ind1)/365 a(ind1)/365 a(ind2)/365 a(ind2)/365],[-1,1,1,-1],colour_mat5,'FaceAlpha',0.6, 'EdgeColor',colour_mat5,'EdgeAlpha', 0.6,...
    'DisplayName','Vacc. age range');
xticks('')
yticks('')
grid off
axis([0.4984 1.3591 7.1556e-05 1.4222e-04]);
%% calculate # of DH children under 5
% baseline - no vaccine
da = a(2)-a(1);
DH_reduce_prop = NaN(length(a),1);
for ind_age = 1:length(a)
    NN = 1131950+1155574+1116436+993183+893681+1670570+590013+1867579;
    DH_reduce_prop(ind_age,1) = trapz(DH_EE_mat(1:ind_age,1)-DH_EE_mat(1:ind_age,2))*da;
end
DH_reduce = DH_reduce_prop*NN;
figure_setups; hold on
plot(a/365,DH_reduce)
plot([a(ind1)/365 a(ind1)/365],[0 2.5*10^4],'m--','DisplayName','Vaccine age - start')
plot([a(ind2)/365 a(ind2)/365],[0 2.5*10^4],'m--','DisplayName','Vaccine age - end')
grid on
xlabel('Age (years)')
ylabel('\# of $D_H$ children reduced')
xlim([0 20]);
% # of reduced for age < 3 years old
[~,ind_age1] = min(abs(a-9*30));
[~,ind_age2] = min(abs(a-3*365));
num_red = DH_reduce(ind_age2)-DH_reduce(ind_age1-1)
num_red/(NN*trapz(DH_EE_mat(ind_age1-1:ind_age2,1))*da)