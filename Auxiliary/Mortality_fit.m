close all 
clear all

%% Burkina Faso mortality data 
% downloaded from https://apps.who.int/gho/data/view.main.61860?lang=en
% nMx - age-specific death rate between ages x and x+n for both sexes in 2019
% Age Group: <1 year, 1-4 years, 5-9 years, 10-14 years,..., 80-84 years,
% and 85+ years
rawdata = [0.0565423270000000;0.00905588700000000;0.00221164700000000;0.00132328300000000;0.00146924100000000;0.00197064900000000;0.00233488800000000;0.00286785500000000;0.00378169300000000;0.00538959500000000;0.00751589400000000;0.0113435140000000;0.0162343470000000;0.0245338440000000;0.0352633890000000;0.0549726570000000;0.0801675880000000;0.119584264000000;0.189767204000000];

alpha = [0.5 2.5 7 12 17 22 27 32 37 42 47 52 57 62 67 72 77 82 87.5];
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
mu = rawdata;

%% estimate parameter values for the mortality rate function b0 + b1*exp(-b2*a) + b3*exp(b4*a)
modelfun = @(b,x)(b(1) + b(2)*exp(-b(3)*x) + b(4)*exp(b(5)*x));

% initial guess for parameter values
b_0 = [0, 0.05, 0.505, 0.01, 0.055];

options = statset('Display','final','TolFun',1e-12, 'MaxIter', 1e4, 'MaxFunEvals', 1e4);
md1 = fitnlm(alpha, mu, modelfun, b_0, 'Options', options)

% extract estimated parameter values
b_est = md1.Coefficients{:, 1};

% examine residuals
figure(1)
plotResiduals(md1, 'fitted')

% slice plot
% one can drag the vertical dashed line to see the effect of a change in alpha on mu
plotSlice(md1)
%% 
figure_setups; hold on
scatter(alpha*365,rawdata/365)
age = 0:20:80*365; % days
b = b_est;
muH =  b(1) + b(2)*exp(-b(3)*age/365) + b(4)*exp(b(5)*age/365); % natural human mortality rate
muH = muH/365;
plot(age,muH)
b = b_0;
muH =  b(1) + b(2)*exp(-b(3)*age/365) + b(4)*exp(b(5)*age/365); % natural human mortality rate
muH = muH/365;
plot(age,muH)
legend('data','estimate','initial')
