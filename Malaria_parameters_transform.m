global P

a = P.a;

P.c2 = P.c1; % weight for maternal immunity
P.c3 = P.c1; % weight for vaccine-derived immunity

P.muD = 0*ones(size(a));  % disease-induced human mortality rate

P.rho = sigmoid_prob(zeros(size(a)), 'rho');
P.phi = sigmoid_prob(zeros(size(a)), 'phi');
P.psi = sigmoid_prob(zeros(size(a)), 'psi');

%% mortality functions 
muH =  P.b0 + P.b1*exp(-P.b2*a/365) + P.b3*exp(P.b4*a/365); % natural human mortality rate
muH = muH/365;
muH_int_fun = @(age) (age./365).*P.b0 + (P.b1./P.b2).*(1-exp(-P.b2.*age./365)) + (P.b3./P.b4).*(-1+exp(P.b4.*age./365));        

%% fertility rate 
gH_fun = @(age) (2.*P.cc.*normpdf((age./365-P.zz)./P.ww).*normcdf(P.alpha.*(age./365-P.zz)./P.ww)./P.ww)./365/2;
gH =  gH_fun(a); % human fertility rate

%%
P.gM = P.muM*P.MHr; % recruitment rate of mosquitoes;
%%
vb_fun = @(age) P.vb0*(age<5*365); % constant vaccination rate, **if changed, need to update Cac steadystate as well**
vb = vb_fun(a); 
vp_fun = @(age) P.vp0*(age<5*365); % constant vaccination rate, **if changed, need to update Cac steadystate as well**
vp = vp_fun(a); 

P.muH = muH; %
P.gH = gH; %
P.vb = vb;
P.vb_fun = vb_fun;
P.vp = vp;
P.vp_fun = vp_fun;
P.gH_fun = gH_fun;
P.muH_int_fun = muH_int_fun;
P.muH_int = muH_int_fun(a);

% Update the fertility and stable age dist. if balanced option is selected
if P.balance_fertility == 1
    balance_fertility;
end

if P.balance_mortality == 1
    balance_mortality;
end

find_stable_age;
%% plot death
figure_setups;
plot(a/365, P.muH, a/365, P.gH)
title(['balanced birth and death Kenya']);
axis([0 age_max/365 0 2*10^-3]);
%% plot stable age distribution PH
figure_setups;
plot(a/365,P.PH_stable,'-k');
title(['Age Distribution']);
xlabel('age (years)');
title('Kenya')
grid on
axis([0 age_max/365 0 max(P.PH_stable)]);
%% flip
figure_setups;
plot(P.PH_stable,a/365,'-k');
title(['Age Distribution']);
ylabel('age (years)');
title('Kenya')
grid on
axis([0 max(P.PH_stable) 0 age_max/365 ]);
% keyboard