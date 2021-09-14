function find_stable_age()
global P
da = P.da;

b0 = 0;
b1 = 0.05;
b2 = 0.505;
b3 = 0.01;
b4 = 0.055;
muH_int = @(age) (age./365).*b0 + (b1./b2).*(1-exp(-b2.*age./365)) + (b3./b4).*(-1+exp(b4.*age./365)); % natural human mortality rate

cc = 4.6;
zz = 25;
alpha = 28;
ww = 13.5;
gH = @(age) (2.*cc.*normpdf((age./365-zz)./ww).*normcdf(alpha.*(age./365-zz)./ww)./ww)./365; 

% Look for a zero in [LP,RP] to find p_hat for the stable age distribution
% use exact function values to reduce numerical errors
F = @(p) integral(@(age) gH(age).*exp(-p.*age- muH_int(age)),0,Inf) - 1;
options = optimset('TolX',1e-18); % don't show iterations
p0 = [0 10^-3]; % [LP,RP]
p = fzero(F,p0,options);
P.p_hat = p;

%% Construct the stable age distribution using p_hat
a = P.a;
P.muH_int = muH_int(a);
P.Lambda = 1/(da.*trapz(exp(-P.p_hat*a-P.muH_int))); % needed for R0 calculation later so make global
n_tilde = P.Lambda*exp(-P.p_hat*a-P.muH_int);
P.n_tilde = n_tilde; % need this elsewhere in Malaria_IC

%% Update the fertility and stable age dist. if balanced option is selected
if P.balance_fertility == 1
    da_fine = 10;
    a_fine = (0:da_fine:P.age_max)';
    gH_fine = gH(a_fine);
    F2 = @(x) x.*(da_fine.*trapz(x.*exp(-muH_int(a_fine))) - 1); % add zero solution to artificially exclude negatives
    options = optimoptions('fsolve','Display','none','OptimalityTolerance', 1e-25);
    [balanced_births_fine,fval,exitflag,output,jacobian] = fsolve(F2,gH_fine,options); % start from current fertility
    % check if new growth factor is approximately zero
    balanced_births = interp1(a_fine,balanced_births_fine,a);
    F3 = @(p) da.*trapz(balanced_births.*exp(-p*a-P.muH_int)) - 1;
    options = optimset('TolX',1e-25);
    p0 = [-10^-3 10^-3]; % [LP,RP]
    p_new = fzero(F3,p0,options);
    P.p_hat = p_new;
    %     figure_setups;
    %     plot(a,P.gH,'b');
    %     hold on;
    %     plot(a,balanced_births,'-.r');
    %     axis_years(gca,P.age_max);
    %     legend('original fertility','balanced fertility');
    
    % update the fertility with the balanced one and update the stable age
    % dist as well
    P.gH = balanced_births;
    P.Lambda = 1/(da.*trapz(exp(-P.muH_int)));
    n_tilde = P.Lambda*exp(-P.muH_int);
    %Lambda = 1/(da.*trapz(exp(-P.p_hat*a-P.muH_int)));
    %n_tilde = Lambda*exp(-P.p_hat*a-P.muH_int);
    P.n_tilde = n_tilde;
end
end
%%
% disp(['q = ',num2str(P.p_hat)]); % we want this as close to zero as possible
% trapz(a,P.n_tilde) % sanity check, should be = 1 for proper normalization
