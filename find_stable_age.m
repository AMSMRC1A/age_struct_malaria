function find_stable_age()
% find the stable age distribution for a general birth and death rates 

global P

% natural human mortality rate
muH_int = P.muH_int_fun; 
% human fertility rate
gH = P.gH_fun;

% Look for a zero in [LP,RP] to find p_hat for the stable age distribution
% use exact function values to reduce numerical errors
F = @(p) integral(@(age) gH(age).*exp(-p.*age- muH_int(age)),0,Inf) - 1;
options = optimset('TolX',1e-25); % don't show iterations
p0 = [-10^-3 10^-3]; % [LP,RP]
p = fzero(F,p0,options);
P.p_hat = p;

if P.balance_fertility == 1
    if abs(p)>10^-3
        disp('recheck the balance fertility')
    end
end

% Construct the stable age distribution using p_hat
P.K = 1/integral(@(a) exp(-P.p_hat*a-P.muH_int_fun(a)),0,Inf);
P.PH_stable = P.K*exp(-P.p_hat*P.a-P.muH_int);
P.PH_stable_fun = @(a) P.K.*exp(-P.p_hat.*a-muH_int(a));

end
%%
% disp(['q = ',num2str(P.p_hat)]); % we want this as close to zero as possible
% trapz(a,P.n_tilde) % sanity check, should be = 1 for proper normalization
