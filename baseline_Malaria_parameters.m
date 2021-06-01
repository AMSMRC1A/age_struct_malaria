global P
global a % needed for age-dependent parameters

% from parameter table 
P.phi = 0.5; % 0.5 ??
P.rho = 0.5; %  prob. of severely infected, EH -> DH
P.rD = 1/180; % recovery rate for DH
P.rA = 1/360; % recovery rate for AH
P.h = 1/15; % incubation rate in human
P.psi = 0.5; % prob. AH -> DH

P.bm = 0.67; % desired biting rate per mosquito
P.bh = 18; % tolerated biting rate per human

P.betaM = 0.25; % infectivity of mosquitoes
P.betaD = 0.35; % infectivity of DH
P.betaA = 0.03; % infectivity of AH

P.gM = 0.1;  % ?? recruitment rate of mosquitoes
P.muM = 1/10; % natural mortality rate of mosquitoes
P.sigma = 1/10; % incubation rate for mosquitoes

% age-dependent parameters
gH = 1/(50*365)*ones(size(a)); %0.05*ones(size(a)); % natural human birth rate
muH = 0*ones(size(a));% 0.05*ones(size(a)); % natural human mortality rate
muD = 0*ones(size(a));% 0.05*ones(size(a)); % disease-induced human mortality rate

w = 1/50*ones(size(a)); % 1/50 transition rate RH -> SH, may also be a function on time t
v = 0.01*ones(size(a)); % 0.01 ?? vaccination rate

P.muH = muH; % need to parameterize the birth-death rates, so that the ratio M to H is reasonable
P.muD = muD;
P.gH = gH;
P.w = w;
P.v = v;