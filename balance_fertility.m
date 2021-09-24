function balance_fertility
% calculate a balanced fertility rate, given a death rate
global P

if exist('balanced_births.mat','file') % need to re-run the mat file if the death rates are updated!!
    load('balanced_births.mat','fun_balanced_births')
else
    disp('calculating balanced fertility profile...hit enter')
    pause
    gH = P.gH_fun;
    muH_int = P.muH_int_fun;
    da_fine = 5; % fine grid for approximating new birth
    a_fine = (0:da_fine:P.age_max)';
    gH_fine = gH(a_fine);
    F2 = @(x) x.*(da_fine.*trapz(x.*exp(-muH_int(a_fine))) - 1); % add zero solution to artificially exclude negatives
    options = optimoptions('fsolve','Display','none','OptimalityTolerance', 1e-25);
    [balanced_births_fine,fval,exitflag,output,jacobian] = fsolve(F2,gH_fine,options); % start from current fertility
    fun_balanced_births = griddedInterpolant(a_fine,balanced_births_fine);
    save('balanced_births.mat','fun_balanced_births');
end

P.gH_fun = fun_balanced_births;
P.gH = fun_balanced_births(P.a);

%% check if new growth factor is approximately zero
% F3 = @(p) da.*trapz(balanced_births.*exp(-p*a-P.muH_int)) - 1;
% options = optimset('TolX',1e-25);
% p0 = [-10^-3 10^-3]; % [LP,RP]
% p_new = fzero(F3,p0,options);
% P.p_hat = p_new;

%     figure_setups;
%     plot(a,P.gH,'b');
%     hold on;
%     plot(a,balanced_births,'-.r');
%     axis_years(gca,P.age_max);
%     legend('original fertility','balanced fertility');
end