function find_stable_age()
% find the stable age distribution for a general birth and death rates 

global P

% Look for a zero in [LP,RP] to find p_hat for the stable age distribution
% use exact function values to reduce numerical errors
F = @(p) P.da.*trapz(P.gH.*exp(-p*P.a-P.muH_int)) - 1;
options = optimset('TolX',1e-25); % don't show iterations
p0 = [-10^-3 10^-3]; % [LP,RP]
p = fzero(F,p0,options);
P.p_hat = p;
if P.balance_fertility == 1 || P.balance_mortality == 1
    if abs(p)>10^-3
        disp('recheck the balance fertility')
        keyboard
    end
end

% Construct the stable age distribution using p_hat
% P.K = 1/integral(@(a) exp(-P.p_hat*a-P.muH_int_fun(a)),0,Inf);
P.K = 1/trapz(P.a, exp(-P.p_hat.*P.a-P.muH_int));
P.PH_stable = P.K*exp(-P.p_hat*P.a-P.muH_int);
P.PH_stable_fun = @(a) P.K.*exp(-P.p_hat.*a-P.muH_int_fun(a));

end
%%
% disp(['q = ',num2str(P.p_hat)]); % we want this as close to zero as possible
% trapz(P.a,P.PH_stable) // integral(@(a) P.PH_stable_fun(a),0,Inf) % sanity check, should be = 1 for proper normalization
