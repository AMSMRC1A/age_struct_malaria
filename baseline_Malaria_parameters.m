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
P.cS = 0.75; % SH weight
P.cE = 0.1; % EH weight 
P.cA = 0.1; % AH weight
P.cD = 0.05; % DH weight
P.cV = 0; % weight for vaccination
P.m = 1; % fraction of new-born immunity relative to mother’s

%% mosquito related parameters/rates
P.bm = 2;%7;%0.67; % desired biting rate per mosquito
P.bh = 0.5;%18; % tolerated biting rate per human

P.betaM = 0.125; % infectivity of mosquitoes
P.betaD = 0.35; % infectivity of DH
P.betaA = 0.03; % infectivity of AH

P.gM = 0.5;  % ?? recruitment rate of mosquitoes
P.muM = 1/10; % natural mortality rate of mosquitoes
P.sigma = 1/10; % incubation rate for mosquitoes

%% mortality functions (placeholder parameters)
P.b0 = 0;
P.b1 = 0.05;
P.b2 = 0.505;
P.b3 = 0.01;
P.b4 = 0.055;
muH =  P.b0 + P.b1*exp(-P.b2*a/365) + P.b3*exp(P.b4*a/365); % natural human mortality rate
muH = muH/365;
muH_int_fun = @(age) (age./365).*P.b0 + (P.b1./P.b2).*(1-exp(-P.b2.*age./365)) + (P.b3./P.b4).*(-1+exp(P.b4.*age./365));        
%% fertility rate (placeholder parameters)
P.cc = 4.6;
P.zz = 25;
P.alpha = 28;
P.ww = 13.5;
gH =  2*P.cc.*normpdf((a/365-P.zz)/P.ww).*normcdf(P.alpha*(a/365-P.zz)/P.ww)/P.ww; % 0.05*ones(size(a)); % human fertility rate
gH = gH/365;
gH_fun = @(age) (2.*P.cc.*normpdf((age./365-P.zz)./P.ww).*normcdf(P.alpha.*(age./365-P.zz)./P.ww)./P.ww)./365;
%%
muD = 0*ones(size(a));% 0.05*ones(size(a)); % disease-induced human mortality rate
w = 1/180*ones(size(a)); % 1/50 transition rate RH -> SH, may also be a function on time t
v = 0*ones(size(a)); % 0.01 ?? vaccination rate

P.muH = muH; %
P.muD = muD; % 
P.gH = gH; %
P.w = w;
P.v = v;
P.gH_fun = gH_fun;
P.muH_int_fun = muH_int_fun;