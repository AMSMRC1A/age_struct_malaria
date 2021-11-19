% close all 
clear all
% clc
%% Kenya mortality data 
% downloaded from GHO life tables
% nMx - age-specific death rate between ages x and x+n for both sexes in 2019
% Age Group: <1 year, 1-4 years, 5-9 years, 10-14 years,..., 80-84 years and 85+ years
nMx = [0.033586496	0.002975525	0.001067792	0.000915904	0.001320595	0.001742293	0.002260166	0.003252213	0.004840494	0.006919824	0.00937546	0.012254916	0.015811554	0.021992071	0.031078758	0.047578104	0.072404337	0.112801155 0.195538827];
%%
mu = nMx/365;%-log(1-nqx)/365;
alpha = [0.5 2.5 7 12 17 22 27 32 37 42 47 52 57 62 67 72 77 82 88]*365;

%% estimate parameter values for the mortality rate function b0 + b1*exp(-b2*a) + b3*exp(b4*a)
modelfun = @(b,x) (b(1) + b(2)*exp(-b(3)*x./365) + b(4)*exp(b(5)*x./365))./365;

% initial guess for parameter values
b_0 = [0, 0.05, 0.505, 0.01, 0.055];

options = statset('Display','final','TolFun',1e-12, 'MaxIter', 1e4, 'MaxFunEvals', 1e4);
md1 = fitnlm(alpha, mu, modelfun, b_0, 'Options', options)

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
scatter(alpha/365,mu,'ko')
age = 0:20:100*365; % days
muH =  modelfun(b_est,age);
plot(age/365,muH)
muH =  modelfun(b_0,age);
plot(age/365,muH)
legend('data','estimate','initial')
title('Mortality fit Kenya')
xlabel('age')
ylabel('daily death rate')

