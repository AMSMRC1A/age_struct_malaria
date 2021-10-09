%% Plot Sigmoidal Probability Function
clear all
% close all
clc

cmin = 0.1; % minimum value
cmax = 0.6; % maximum value
c = 1; % constant at half max + min
k = 2; % Hill coefficient (steepness)

L = 25; % effective range is [0,L]
f_0 = 0.1; % value at zero
f_1 = 0.9; % value at L (function saturates to this value)
t_2 = 0.2; % threshold value (as a fraction of L)
s_2 = 0.15; % sigmoid steepness, smaller is steeper

x = 0:.01:10; %input

% y =  (f_0 + (f_1-f_0)./(- exp(-(x./L-t_2)./s_2))); 
y = (cmax-cmin)*(1-c^k./(c^k+x.^k))+cmin;

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
