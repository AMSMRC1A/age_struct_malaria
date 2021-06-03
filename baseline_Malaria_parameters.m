global P
global a % needed for age-dependent parameters

P.rD = 1/180; % recovery rate for DH
P.rA = 1/360; % recovery rate for AH
P.h = 1/15; % incubation rate in human

% immunity parameters
P.ds = 5;
P.dm = 0.25;
P.c1 = 1;
P.c2 = 0.1;
P.Cm0 = 0.4; % intial maternal immunity level at age = 0
P.Cm = P.Cm0*exp(-a./P.dm); % maternal immunity level at age a

P.bm = 0.67; % desired biting rate per mosquito
P.bh = 18; % tolerated biting rate per human

P.betaM = 0.25; % infectivity of mosquitoes
P.betaD = 0.35; % infectivity of DH
P.betaA = 0.03; % infectivity of AH

P.gM = 0.1;  % ?? recruitment rate of mosquitoes
P.muM = 1/10; % natural mortality rate of mosquitoes
P.sigma = 1/10; % incubation rate for mosquitoes

%% mortality functions (placeholder parameters)
b0 = 0;
b1 = 0.05;
b2 = 0.505;
b3 = 0.01;
b4 = 0.055;
muH =  0*ones(size(a));%b0 + b1*exp(-b2*a/365) + b3*exp(b4*a/365); % natural human mortality rate
%% fertility rate (placeholder parameters)
cc = 4.6;
zz = 25;
alpha = 28;
ww = 13.5;
gH =  0*ones(size(a));%2*cc.*normpdf((a/365-zz)/ww).*normcdf(alpha*(a/365-zz)/ww)/ww; % 0.05*ones(size(a)); % human fertility rate

%%
muD = 0*ones(size(a));% 0.05*ones(size(a)); % disease-induced human mortality rate

w = 1/50*ones(size(a)); % 1/50 transition rate RH -> SH, may also be a function on time t
v = 0.01*ones(size(a)); % 0.01 ?? vaccination rate

P.muH = muH; %
P.muD = muD; % 
P.gH = gH; %
P.w = w;
P.v = v;