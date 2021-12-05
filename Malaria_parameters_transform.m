function Malaria_parameters_transform
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
%% vaccination functions - baseline vaccine (default = 0)
vb_fun = @(age) P.vb0.*ones(size(age));
vb = vb_fun(a);
vp_fun = @(age) P.vp0.*ones(size(age));
vp = vp_fun(a);

% approximation of theta at DFE - needed for analytical purpose: DFE, R0, bifurcation
pi_fun = @(x) P.w+P.e.*vp_fun(x);
pi_int_a = intf(pi_fun,P.a);
pi_int_fun = @(x) interp1(P.a,pi_int_a,x);
exp_pi_int_a = intf(@(x) P.w.*exp(pi_int_fun(x)),P.a);
exp_pi_int = @(x) interp1(P.a,exp_pi_int_a,x);
theta_fun = @(x) exp(-pi_int_fun(x)).*(1+exp_pi_int(x));
theta0 = theta_fun(P.a);

P.muH = muH;
P.gH = gH;
P.vb = vb;
P.vb_fun = vb_fun;
P.vp = vp;
P.vp_fun = vp_fun;
P.theta_fun = theta_fun;
P.theta = theta0;
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

%% plot calibrated fertility & death
% figure_setups; hold on
% grid on
% plot(a/365, P.muH, a/365, P.gH)
% xlabel('Age (years)')
% ylabel('Daily rates')
% %title(['Balanced demographics']);
% legend('Mortality rate $\mu_H(\alpha)$','Fertility rate $g_H(\alpha)$')
% %% plot stable age distribution PH
% figure_setups;
% plot(a/365,P.PH_stable,'-k');
% %title(['Stable age distribution']);
% legend('$P_H(\alpha)$')
% xlabel('Age (years)');
% ylabel('Population density')
% grid on
% axis([0 age_max/365 0 max(P.PH_stable)]);
%
% %% flip
% figure_setups;
% plot(P.PH_stable,a/365,'-k');
% title(['Age Distribution']);
% ylabel('age (years)');
% title('Kenya')
% grid on
% axis([0 max(P.PH_stable) 0 age_max/365 ]);
% % keyboard
end

function intfx = intf(f,xs)
ns = length(xs);
fx = f(xs);
e = ones(ns,ns+1)/2;
e(2:end,1:end-1)=1;e(1,end)=0;
A = spdiags(e,-ns:0,ns,ns);
intfx= A*fx;
end

