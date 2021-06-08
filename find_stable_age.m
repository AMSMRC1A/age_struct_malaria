%% 
global P
global a

da = 1;
a_max = 50*365; % in years
a = 0:da:a_max;

% setup the birth and mortality functions
baseline_Malaria_parameters;
muH_int = (a/365)*b0 + (b1/b2)*(1-exp(-b2*a/365)) + (b3/b4)*(-1+exp(b4*a/365));
P.muH_int = muH_int;%/365;

% Look for a zero in [LP,RP] to find p_hat for the stable age distribution
F = @(p) da.*trapz(P.gH.*exp(-p*a-P.muH_int)) - 1;
options = optimset('Display','iter','TolX',1e-8); % show iterations
p0 = [0 1]; % [LP,RP]
p = fzero(F,p0,options);

%% Construct the stable age distribution using p_hat
Lambda = 1/(da.*trapz(exp(-p*a-P.muH_int)));
n_tilde = Lambda*exp(-p*a-P.muH_int); 

figure_setups;
plot(a,n_tilde);
axis_years(gca,a_max)
xlabel('age');
ylabel('pop. density')
title('Stable Age Distribution (?)')
trapz(a,n_tilde) % sanity check, should be = 1

