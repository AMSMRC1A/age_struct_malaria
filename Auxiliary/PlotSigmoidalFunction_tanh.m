%% Plot Sigmoidal Probability Function
% clear all
close all
% clc

L = 6; % effective range is [0,L]
f_0 = 0.01; % value at zero
f_1 = 1; % value at L (function saturates to this value)
t_2 = 0.993686837551587; % threshold value (as a fraction of L)
s_2 = 0.030339261357744; % sigmoid steepness, smaller is steeper
x = 1.5:0.01:L; %input

y_inc = f_0 + (f_1-f_0)./(1 + exp(-(x./L-t_2)./s_2));
y_dec = f_1 + (f_0-f_1)./(1 + exp(-(x./L-t_2)./s_2)); % switch f0 f1

figure_setups; hold on
plot(x,y_inc,'linewidth',2)
plot(x,y_dec,'linewidth',2)
plot(x,f_0*ones(size(x)),'k--')
plot(x,f_1*ones(size(x)),'k--')
axis([min(x) max(x) 0, 1])

