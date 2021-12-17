%% plot Fig3 a - heatmap -- results from "run_parameter_fit.m"

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

Malaria_parameters_baseline;
load(['Results/Heatmap.mat'],'var_list','final_EIR','final_immunity')
% plot heatmap of final immunity(age, EIR)
figure_setups_3; hold on
imagesc(a/365,final_EIR,final_immunity');%
% indicates two cuts at 84.614242650691182 and 44.658135926297568
plot([0 20],[84.614242650691182 84.614242650691182],'k--')
plot([0 20],[44.658135926297568 44.658135926297568],'k--')
xlim([0 20])
xlabel('Age (years)')
ylabel('aEIR')
% title(['Immunity levels, feedback = ',num2str(immunity_feedback)]);
title('Immunity level per person');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:2:10,'TickLabelInterpreter','latex');
set(gcf, 'Renderer', 'Painters');
print('Figures/result_heatmap.eps', '-depsc','-r300')
