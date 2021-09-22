global P

a = P.a;

%% mortality functions (placeholder parameters)
muH =  P.b0 + P.b1*exp(-P.b2*a/365) + P.b3*exp(P.b4*a/365); % natural human mortality rate
muH = muH/365;
muH_int_fun = @(age) (age./365).*P.b0 + (P.b1./P.b2).*(1-exp(-P.b2.*age./365)) + (P.b3./P.b4).*(-1+exp(P.b4.*age./365));        

%% fertility rate (placeholder parameters)
gH_fun = @(age) (2.*P.cc.*normpdf((age./365-P.zz)./P.ww).*normcdf(P.alpha.*(age./365-P.zz)./P.ww)./P.ww)./365;
gH =  gH_fun(a); % human fertility rate

%%
v_fun = @(age) P.v0*ones(size(age)); % constant vaccination rate, **if changed, need to update Cac steadystate as well**
v = v_fun(a); 

P.muH = muH; %
P.gH = gH; %
P.v = v;
P.v_fun = v_fun;
P.gH_fun = gH_fun;
P.muH_int_fun = muH_int_fun;