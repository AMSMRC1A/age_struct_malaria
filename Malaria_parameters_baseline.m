global P
%% 
P.verbose = 1; % turn on the warning messages. Error messages from the check routines will display regardless
P.balance_fertility = 0; % balanced fertility or not
P.balance_mortality = 1; % balanced mortality or not

%% vaccine related parameters
P.vb0 = 0; P.vb0_lower = 0.1;  P.vb0_upper = 0.8; % vaccination rate for boosting immunity (Cv) (boosting)
P.dv = 5*365; % Half-life of vaccine-boosted immunity (Cv)
%
P.vp0 = 0; P.vp0_lower = 0.1;  P.vp0_upper = 0.8; % vaccination rate for sterilizing immunity (VH) (infection protection)
P.w = 1/(0.66*365); % Waning rate for infection-protection immunity (VH) for children 
P.e = 0.73; % Vaccine efficacy for the infection-protection immunity (VH) for children 

%%
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
P.phi_f_0 = 0.01; 
P.phi_f_1 = 1;
P.phi_s_2 = 2.432473210664639; 
P.phi_r_2 = 1.277554702429119; 
% rho = psi
P.rho_f_0 = 0.01; 
P.rho_f_1 = 1; 
P.rho_s_2 = 3.186642715992263; 
P.rho_r_2 = 1.030298116795388; 
P.psi_f_0 = 0.01; 
P.psi_f_1 = 1; 
P.psi_s_2 = 3.186642715992263; 
P.psi_r_2 = 1.030298116795388; 

%% mosquito related parameters/rates
P.bh = 5; P.bh_lower = 2; P.bh_upper = 19;% tolerated biting rate per human
P.bm = 0.6; P.bm_lower = 0.4; P.bm_upper = 0.7; % desired biting rate per mosquito

P.betaM = 0.25; P.betaM_lower = 0.125; P.betaM_upper = 0.3; % infectivity of mosquitoes
P.betaD = 0.35; P.betaD_lower = 0.3; P.betaD_upper = 0.4; % infectivity of DH
P.betaA = 0.03; P.betaA_lower = 0.02; P.betaA_upper = 0.04; % infectivity of AH

P.muM = 1/10; P.muM_lower = 1/15; P.muM_upper = 1/7; % natural mortality rate of mosquitoes
P.MHr = 5; P.MHr_lower = 2; P.MHr_upper = 10; % mosquito-human population ratio gm = muM*ratio
P.sigma = 1/15; P.sigma_lower = 1/30; P.sigma_upper = 1/10; % incubation rate for mosquitoes

%% mortality rate parameters (initial estimates from raw data)
% use rawdata - Kenya - nMx
P.b0 = 0.0024214446844162;
P.b1 = 0.0887924178445357;
P.b2 = 2.09862983723212;
P.b3 = 6.87709371762464e-05;
P.b4 = 0.0901695513967616;

%% fertility rate parameters (initial estimates)
% use rawdata - Kenya
P.cc = 4.024086261410830;
P.zz = 17.963601264000353;
P.alpha = 4.083610527018673;
P.ww = 13.196127635937707;

%%
Malaria_parameters_transform;