global P
P.balance_fertility = 1; % 0 to keep original fertility, 1 for balanced birth rate so that pop. growth is zero

P.rD = 1/180; % recovery rate for DH
P.rA = 1/360; % recovery rate for AH
P.h = 1/15; % incubation rate in human

%% immunity parameters/rates
P.dac = 5*365; % half life of acquired immunity
P.dm = 0.25*365; % half life of maternal immunity
P.c1 = 1; % weight for acquired immunity
P.cS = 0.75; % SH weight
P.cE = 0.1; % EH weight  ~~ AH
P.cA = 0.1; % AH weight
P.cD = 0.05; % DH weight
P.cV = 0.75; % weight for vaccination ~~ SH
P.m = 1; % fraction of new-born immunity relative to mother’s
P.uc = 10; % Duration in which immunity is not boosted
%% progression probabilities parameters, sigmoid parameters
% fitted values     
P.L = 7.974031085468034; % effective range is [0,L]
P.phi_f_0 = 0.01; % value at zero
P.phi_f_1 = 1; % value at L (function saturates to this value)
P.phi_t_2 = 0.200455981024660; % threshold value (as a fraction of L)
P.phi_s_2 = 0.741327847456465; % sigmoid steepness, smaller is steeper
% rho = psi
P.rho_f_0 = 0.01; % value at zero
P.rho_f_1 = 1; % value at L (function saturates to this value)
P.rho_t_2 = 0.225836736227440; % threshold value (as a fraction of L)
P.rho_s_2 = 0.227377462176943; % sigmoid steepness, smaller is steeper

P.psi_f_0 = 0.01; % value at zero
P.psi_f_1 = 1; % value at L (function saturates to this value)
P.psi_t_2 = 0.225836736227440; % threshold value (as a fraction of L)
P.psi_s_2 = 0.227377462176943; % sigmoid steepness, smaller is steeper

%% mosquito related parameters/rates
P.bh = 5; P.bh_lower = 2; P.bh_upper = 19;% tolerated biting rate per human
P.bm = 0.6; P.bm_lower = 0.4; P.bm_upper = 0.7; % desired biting rate per mosquito

% P.betaM = 0.08;
P.betaM = 0.25; P.betaM_lower = 0.125; P.betaM_upper = 0.3; % infectivity of mosquitoes
P.betaD = 0.35; P.betaD_lower = 0.3; P.betaD_upper = 0.4; % infectivity of DH
P.betaA = 0.03; P.betaA_lower = 0.02; P.betaA_upper = 0.04; % infectivity of AH

P.muM = 1/10; P.muM_lower = 1/15; P.muM_upper = 1/7; % natural mortality rate of mosquitoes
P.MHr = 5; P.MHr_lower = 2; P.MHr_upper = 10; % mosquito-human population ratio gm = muM*ratio
P.sigma = 1/15; P.sigma_lower = 1/30; P.sigma_upper = 1/10; % incubation rate for mosquitoes

%% mortality functions (placeholder parameters)
% P.b0 = 0;
% P.b1 = 0.05;
% P.b2 = 0.505;
% P.b3 = 0.01;
% P.b4 = 0.055;

P.b0 = 0.0011;
P.b1 = 0.0902;
P.b2 = 0.9788;
P.b3 = 0.0001;
P.b4 = 0.0824;            

%% fertility rate (placeholder parameters) -- probability not used later...
P.cc = 4.6;
P.zz = 25;
P.alpha = 28;
P.ww = 13.5;

%%
P.v0 = 0; P.v0_lower = 0.1;  P.v0_upper = 0.8; % vaccination rate parameter (constant rate in age)
Malaria_parameters_transform;