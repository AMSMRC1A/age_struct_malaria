function balance_fertility
% calculate a balanced fertility rate, given a death rate
global P

if exist('balanced_births.mat','file') % need to re-run the mat file if the death rates are updated!!
    load('balanced_births.mat','fun_balanced_births')
    if P.verbose==1; disp('Previously calculated balanced birth rate loaded.'); end
else
    if P.verbose==1; disp('calculating new balanced fertility profile.'); end
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

gH_old = P.gH;

P.gH_fun = fun_balanced_births;
P.gH = fun_balanced_births(P.a);

%% check if new growth factor is approximately zero
% F3 = @(p) P.da.*trapz(P.gH.*exp(-p*P.a-P.muH_int)) - 1;
% options = optimset('TolX',1e-25);
% p0 = [-10^-3 10^-3]; % [LP,RP]
% p_new = fzero(F3,p0,options);
% P.p_hat = p_new;
% %
%     figure_setups;
%     plot(P.a/365,gH_old,'b');
%     hold on;
%     plot(P.a/365,P.gH,'-.r');
%     legend('original fertility','balanced fertility');
% %     keyboard
end