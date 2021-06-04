global P
global a % needed for age-dependent parameters

% All rates are parameterized in days


P.rD = 1/180/P.scale; % recovery rate for DH
P.rA = 1/360/P.scale; % recovery rate for AH
P.h = 1/15/P.scale; % incubation rate in human

% immunity parameters
P.ds = 5*365*P.scale; % mean immunity length
P.dm = 0.25*365*P.scale; % mean maternal immunity length
P.c1 = 1;
P.c2 = 1;
P.c3 = 1;

P.bm = .67/P.scale;%0.67; % desired biting rate per mosquito
P.bh = 18/P.scale; % tolerated biting rate per human

P.betaM = 0.25; % infectivity of mosquitoes
P.betaD = 0.35; % infectivity of DH
P.betaA = 0.03; % infectivity of AH

P.gM = 0.5/P.scale;  % ?? recruitment rate of mosquitoes
P.muM = 1/10/P.scale; % natural mortality rate of mosquitoes
P.sigma = 1/10/P.scale; % incubation rate for mosquitoes

%% mortality functions (placeholder parameters)
b0 = 0;
b1 = 0.05;
b2 = 0.505;
b3 = 0.01;
b4 = 0.055;
muH =  b0 + b1*exp(-b2*a/365/P.scale) + b3*exp(b4*a/365/P.scale); % natural human mortality rate
muH = muH/365/P.scale;
%% fertility rate (placeholder parameters)
cc = 4.6;
zz = 25;
alpha = 28;
ww = 13.5;
gH =  2*cc.*normpdf((a/365/P.scale-zz)/ww).*normcdf(alpha*(a/365/P.scale-zz)/ww)/ww; % 0.05*ones(size(a)); % human fertility rate
gH = gH/365/P.scale;
%%
muD = 0/P.scale*ones(size(a));% 0.05*ones(size(a)); % disease-induced human mortality rate

w = 1/180/P.scale*ones(size(a)); % 1/50 transition rate RH -> SH, may also be a function on time t
v = 0/P.scale*ones(size(a)); % 0.01 ?? vaccination rate

P.muH = muH; %
P.muD = muD; % 
P.gH = gH; %
P.w = w;
P.v = v;