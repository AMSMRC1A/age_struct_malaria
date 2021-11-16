% close all 
clear all
clc
%% Burkina Faso mortality data 
% downloaded from https://apps.who.int/gho/data/view.main.61860?lang=en
% nMx - age-specific death rate between ages x and x+n for both sexes in 2019
% Age Group: <1 year, 1-4 years, 5-9 years, 10-14 years,..., 80-84 years,
% and 85+ years
% nMx = [0.0565423270000000;0.00905588700000000;0.00221164700000000;0.00132328300000000;0.00146924100000000;0.00197064900000000;0.00233488800000000;0.00286785500000000;0.00378169300000000;0.00538959500000000;0.00751589400000000;0.0113435140000000;0.0162343470000000;0.0245338440000000;0.0352633890000000;0.0549726570000000;0.0801675880000000;0.119584264000000;0.189767204000000]./365;
% nqx = [0.054389607	0.035453008	0.010997429	0.006594598	0.007319319 0.00980494	0.011606689	0.014237201	0.018731373	0.026589707	0.036886386  0.055153485	0.0780058	0.115580155	0.162032413	0.241652576	0.333915034  0.460307488]; 
%% Kenya
nMx = [0.033586496	0.002975525	0.001067792	0.000915904	0.001320595	0.001742293	0.002260166	0.003252213	0.004840494	0.006919824	0.00937546	0.012254916	0.015811554	0.021992071	0.031078758	0.047578104	0.072404337	0.112801155	0.195538827 1];
% nqx = [0.032814997	0.011817706	0.005324746	0.004569058	0.006581245	0.008673684	0.011237335	0.016129918	0.023913094	0.034010747	0.045803721	0.0594531	0.076051537	0.104229781	0.144190626	0.212602463	0.306535446	0.439941112];
%%
mu = nMx/365;%-log(1-nqx)/365;
alpha = [0.5 2.5 7 12 17 22 27 32 37 42 47 52 57 62 67 72 77 82]*365;
% alpha = [0 3 7.5 12.5 17.5 22.5 27.5 32.5 37.5 42.5 47.5 52.5 57.5 62.5 67.5 72.5 77.5 82.5 87.5];%[0:1:100]';
%nn = length(alpha);
%mu = zeros(nn,1);
%mu(1) = rawdata(1);
%mu(2:5) = rawdata(2)*ones(1, 4);
%nj = floor(84/5);
%for j = 1:nj
%    mu(j*5+1:j*5+5) = rawdata(j+2)*ones(1, 5);
%end
%mu(nj*5+6:nn) = rawdata(end)*ones(1,nn-nj*5-5);

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

