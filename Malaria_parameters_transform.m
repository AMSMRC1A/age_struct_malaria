global P

a = P.a;

P.muD = 0*ones(size(a));  % disease-induced human mortality rate

P.rho = sigmoid_prob(zeros(size(a)), 'rho');
P.phi = sigmoid_prob(zeros(size(a)), 'phi');
P.psi = sigmoid_prob(zeros(size(a)), 'psi');

%% mortality functions (Burkina Faso parameters)
muH =  P.b0 + P.b1*exp(-P.b2*a/365) + P.b3*exp(P.b4*a/365); % natural human mortality rate
muH = muH/365;
muH_int_fun = @(age) (age./365).*P.b0 + (P.b1./P.b2).*(1-exp(-P.b2.*age./365)) + (P.b3./P.b4).*(-1+exp(P.b4.*age./365));        

%% fertility rate (Burkina Faso parameters)
gH_fun = @(age) (2.*P.cc.*normpdf((age./365-P.zz)./P.ww).*normcdf(P.alpha.*(age./365-P.zz)./P.ww)./P.ww)./365;
gH =  gH_fun(a); % human fertility rate

%%
P.gM = P.muM*P.MHr; % recruitment rate of mosquitoes;
%%
v_fun = @(age) P.v0*ones(size(age)); % constant vaccination rate, **if changed, need to update Cac steadystate as well**
v = v_fun(a); 

P.muH = muH; %
P.gH = gH; %
P.v = v;
P.v_fun = v_fun;
P.gH_fun = gH_fun;
P.muH_int_fun = muH_int_fun;
P.muH_int = muH_int_fun(a);

% Update the fertility and stable age dist. if balanced option is selected
if P.balance_fertility == 1
    balance_fertility;
end

find_stable_age;
