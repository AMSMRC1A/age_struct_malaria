global P

a = P.a;

P.rD = 1/180; % recovery rate for DH
P.rA = 1/360; % recovery rate for AH
P.h = 1/15; % incubation rate in human

%% immunity parameters/rates
P.dac = 5*365; % half life of acquired immunity
P.dm = 1*365; % half life of maternal immunity
P.c1 = 1; % weight for acquired immunity
P.c2 = 1; % weight for maternal immunity
P.cS = 1; % SH weight % previously 2x10^5
P.cE = 1; % EH weight 
P.cA = 1; % AH weight % previously 10^3;
P.cD = 1; % DH weight
P.cV = 0; % weight for vaccination
P.m = 1; % fraction of new-born immunity relative to mother’s

%% mosquito related parameters/rates
P.bm = 3;%7;%0.67; % desired biting rate per mosquito
P.bh = 5;%18; % tolerated biting rate per human

P.betaM = 0.01;%0.25; % infectivity of mosquitoes
P.betaD = 0.01;%0.35; % infectivity of DH
P.betaA = 0.01;%0.03; % infectivity of AH

P.gM = 0.5;  % ?? recruitment rate of mosquitoes
P.muM = 1/10; % natural mortality rate of mosquitoes
P.sigma = 1/10; % incubation rate for mosquitoes

%% mortality functions (placeholder parameters)
b0 = 0;
b1 = 0.05;
b2 = 0.505;
b3 = 0.01;
b4 = 0.055;
muH =  b0 + b1*exp(-b2*a/365) + b3*exp(b4*a/365); % natural human mortality rate
muH = muH/365;
%% fertility rate (placeholder parameters)
cc = 4.6;
zz = 25;
alpha = 28;
ww = 13.5;
gH =  2*cc.*normpdf((a/365-zz)/ww).*normcdf(alpha*(a/365-zz)/ww)/ww; % 0.05*ones(size(a)); % human fertility rate
gH = gH/365;
%%
muD = 0*ones(size(a));% 0.05*ones(size(a)); % disease-induced human mortality rate
w = 1/180*ones(size(a)); % 1/50 transition rate RH -> SH, may also be a function on time t
v = 0*ones(size(a)); % 0.01 ?? vaccination rate

P.muH = muH; %
P.muD = muD; % 
P.gH = gH; %
P.w = w;
P.v = v;