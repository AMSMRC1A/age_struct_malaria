% close all 
clear all
% clc
%% Kenya per-capita annual birth rate
rawdata = [0 0 0 0 96 206 183 148 100 38 9 0 0 0 0 0 0 0 0]/1000/365/2;
%% years old
% Age Group: <1 year, 1-4 years, 5-9 years, 10-14 years,..., 80-84 years and 85+ years
alpha = [0.5 2.5 7 12 17 22 27 32 37 42 47 52 57 62 67 72 77 82 87.5].*365;
gH = rawdata;
%% estimate parameter values for the mortality rate function (2.*cc.*normpdf((age-zz)./ww).*normcdf(alpha.*(age-zz)./ww)./ww);
modelfun = @(b,x) (2.*b(1).*normpdf((x./365-b(2))./b(4)).*normcdf(b(3).*(x./365-b(2))./b(4))./b(4))/365/2;

% initial guess for parameter values
b_0 = [4.6, 25, 28, 13.5];

options = statset('Display','final','TolFun',1e-16, 'MaxIter', 1e6, 'MaxFunEvals', 1e6);
md1 = fitnlm(alpha, gH, modelfun, b_0, 'Options', options)

% extract estimated parameter values
b_est = md1.Coefficients{:, 1};

% examine residuals
% figure(1)
% plotResiduals(md1, 'fitted')

% slice plot
% one can drag the vertical dashed line to see the effect of a change in alpha on mu
% plotSlice(md1)
%% 
figure_setups; hold on
scatter(alpha/365,gH)
age = 0:0.1:100*365; % years
b = b_est;
gH_fit =  modelfun(b_est,age); % natural human mortality rate
plot(age/365,gH_fit)
plot(age/365,modelfun(b_0,age));
legend('data','estimate','initial')


