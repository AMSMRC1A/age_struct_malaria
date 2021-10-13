% return an interpolated surface function as data for fitting
% F(age(years), EIR) = susceptibility (~ rho) from Filipe's paper
close all
clear all
clc
x1 = load('Data/EIR_1.mat'); 
x2 = load('Data/EIR_10.mat');
x3 = load('Data/EIR_20.mat');
x4 = load('Data/EIR_50.mat');
x5 = load('Data/EIR_100.mat');
Data_1 = [x1.EIR_1(:,1), 1*ones(size(x1.EIR_1(:,1))), x1.EIR_1(:,2)];
Data_2 = [x2.EIR_10(:,1), 10*ones(size(x2.EIR_10(:,1))), x2.EIR_10(:,2)];
Data_3 = [x3.EIR_20(:,1), 20*ones(size(x3.EIR_20(:,1))), x3.EIR_20(:,2)];
Data_4 = [x4.EIR_50(:,1), 50*ones(size(x4.EIR_50(:,1))), x4.EIR_50(:,2)];
Data_5 = [x5.EIR_100(:,1), 100*ones(size(x5.EIR_100(:,1))), x5.EIR_100(:,2)];
Data_mat = [Data_1;Data_2;Data_3;Data_4;Data_5];
x = Data_mat(:,1); y = Data_mat(:,2); v = Data_mat(:,3);
F = scatteredInterpolant(x,y,v,'linear','nearest');
figure_setups; hold on;
plot3(x,y,v,'.')
[xq,yq] = meshgrid(0:1:50,1:10:150);
vq = F(xq,yq);
mesh(xq,yq,vq)
save('Data/F.mat','F'); 