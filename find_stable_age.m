%%
global P

da = P.da;

% setup the birth and mortality functions
baseline_Malaria_parameters;
muH_int = (a/365)*b0 + (b1/b2)*(1-exp(-b2*a/365)) + (b3/b4)*(-1+exp(b4*a/365));
P.muH_int = muH_int;

% Look for a zero in [LP,RP] to find p_hat for the stable age distribution
F = @(p) da.*trapz(P.gH.*exp(-p*a-P.muH_int)) - 1;
% options = optimset('Display','iter','TolX',1e-12); % show iterations
options = optimset('TolX',1e-18); % don't show iterations
p0 = [0 1]; % [LP,RP]
p = fzero(F,p0,options);
P.p_hat = p;

%% Construct the stable age distribution using p_hat
P.Lambda = 1/(da.*trapz(exp(-P.p_hat*a-P.muH_int))); % needed for R0 calculation later so make global
n_tilde = P.Lambda*exp(-P.p_hat*a-P.muH_int);
P.n_tilde = n_tilde; % need this elsewhere in Malaria_IC

%% Update the fertility and stable age dist. if balanced option is selected
if P.balance_fertility == 1
    F2 = @(x) x.*(da.*trapz(x.*exp(-P.muH_int)) - 1); % add zero solution to artificially exclude negatives
    error_tolerance = 1e-25;
    options = optimoptions('fsolve','Display','none','OptimalityTolerance',...
        error_tolerance);
    [balanced_births,fval,exitflag,output,jacobian] = fsolve(F2,P.gH,options); % start from current fertility
    %balanced_births = sqrt(balanced_births);
    
    % check if new growth factor is approximately zero
    F3 = @(p) da.*trapz(balanced_births.*exp(-p*a-P.muH_int)) - 1;
    options = optimset('TolX',1e-25);
    p0 = [0 1]; % [LP,RP]
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
%%
%figure_setups;
%plot(a,P.n_tilde);
%axis_years(gca,P.age_max)
%xlabel('age');
%ylabel('pop. density')
%title('Stable Age Distribution');
disp(['q = ',num2str(P.p_hat)]); % we want this as close to zero as possible
%trapz(a,P.n_tilde) % sanity check, should be = 1 for proper normalization
