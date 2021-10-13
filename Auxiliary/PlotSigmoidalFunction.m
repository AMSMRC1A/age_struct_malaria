%% Plot Sigmoidal Probability Function
% clear all
% close all
% clc

% 1.426205877278959  10.274569435455705  13.194094606343317   1.273174804387240   4.015779931093367   3.411148478637459

cmin = 0.1; % minimum value
cmax = 1; % maximum value
% c = 1.426205877278959; % constant at half max + min
% k = 1.273174804387240; % Hill coefficient (steepness)
c = 10.274569435455705;
k = 4.015779931093367;

L = 25; % effective range is [0,L]
f_0 = 0.1; % value at zero
f_1 = 0.9; % value at L (function saturates to this value)
t_2 = 0.2; % threshold value (as a fraction of L)
s_2 = 0.15; % sigmoid steepness, smaller is steeper

x = 0:.01:8; %input

% y =  (f_0 + (f_1-f_0)./(- exp(-(x./L-t_2)./s_2))); 
y = (cmax-cmin)*(c^k./(c^k+x.^k))+cmin;

figure_setups;
plot(x,y,'linewidth',2)
hold on
plot([x(1) x(end)],[(cmin+cmax)/2 (cmin+cmax)/2],'k:')
plot([c c],[0 1],'k:')
% plot([x(1) x(end)],[f_0 f_0],'k:')
% plot([x(1) x(end)],[f_1 f_1],'k:')

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
