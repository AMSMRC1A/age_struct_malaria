% return an interpolated surface function as data for fitting
% F(age(years), EIR) = immunity (~ Ctot) from RB's paper, Fig 5

close all
clear all
clc
XY = readtable('age_aEIR.csv','NumHeaderLines',1);  
D = readtable('anti_disease.csv');  ZZ = table2array(D);
x = table2array(XY(:,1)); y = table2array(XY(:,2));

[X,Y] = ndgrid(x,y);
F = griddedInterpolant(X,Y,ZZ);

figure_setups; hold on;
[xq,yq] = ndgrid(min(x):0.1:max(x),min(y):0.1:max(y));
Vq = F(xq,yq);
mesh(xq,yq,Vq)
xlabel('age')
ylabel('aEIR')
zlabel('Ctot')
axis([min(x) max(x) min(y) max(y)])
set(gca,'YDir','normal','Yscale','log');
yticks([2 5 10 25 50 100 200])
grid on
colormap jet
colorbar('Ticks',2.5:0.5:6)
% save('F_RB.mat','F'); 
