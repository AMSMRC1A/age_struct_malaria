%% Plot Sigmoidal Probability Function
% clear all
% close all
% clc

cmin = 0.1; % minimum value
cmax = 1; % maximum value
c = 5;
k = 4.015779931093367;

x = 0:.01:20; %input

y = (cmax-cmin)*(c^k./(c^k+x.^k))+cmin;

figure_setups;
plot(x,y,'linewidth',2)
hold on
plot([x(1) x(end)],[(cmin+cmax)/2 (cmin+cmax)/2],'k:')
plot([c c],[0 1],'k:')
plot([x(1) x(end)],[cmin cmin],'k:')
plot([x(1) x(end)],[cmax cmax],'k:')
hold off
text(c*1.1,0.99,'c')
text(0.9*x(end),cmin+.2*cmin,'cmin')
text(0.9*x(end),cmax+.2*cmin,'cmax')
set(gca,'fontsize',14)
xlabel('input')
ylabel('output')
axis([x(1) x(end) 0 1])
