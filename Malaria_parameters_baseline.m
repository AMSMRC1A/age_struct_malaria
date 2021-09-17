global P

P.balance_fertility = 1; % 0 to keep original fertility, 1 for balanced birth rate so that pop. growth is zero

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
P.cV = 0.1; % weight for vaccination
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

%% fertility rate (placeholder parameters)
P.cc = 4.6;
P.zz = 25;
P.alpha = 28;
P.ww = 13.5;

%%
P.muD = 0;  % disease-induced human mortality rate
P.v0 = 0.5; % vaccination rate parameter (constant rate in age)
P.v0_lower = 0.1;  P.v0_upper = 0.8;
Malaria_parameters_transform;