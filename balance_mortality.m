function balance_mortality
% calculate a balanced mortality rate, given a fertility rate
global P

if exist('balanced_mortality.mat','file') % need to re-run the mat file if the death rates are updated!!
    load('balanced_mortality.mat','fun_balanced_muH_int')
    disp('Previously calculated balanced mortality rate loaded.')

else
    disp('calculating new balanced mortality profile.')    
    
    gH = P.gH_fun;
    da_fine = 5; % fine grid for approximating new birth
    a_fine = (0:da_fine:P.age_max)';
    gH_fine = gH(a_fine);
   
    F2 = @(x) (da_fine.*trapz(gH_fine.*exp(-x.*P.muH_int_fun(a_fine))) - 1); % add zero solution to artificially exclude negatives
    options = optimoptions('fsolve','Display','none','OptimalityTolerance', 1e-25);
    [balanced_coef,fval,exitflag,output,jacobian] = fsolve(F2,1,options); % start from current mortality
    
    balanced_M_fine = balanced_coef.*P.muH_int_fun(a_fine);
    fun_balanced_muH_int = griddedInterpolant(a_fine,balanced_M_fine);
    save('balanced_mortality.mat','fun_balanced_muH_int');
end

muH_old = P.muH;

P.muH_int_fun = fun_balanced_muH_int;
P.muH_int = fun_balanced_muH_int(P.a);
P.muH = (fun_balanced_muH_int(P.a+1)-fun_balanced_muH_int(P.a-1))/2;

%% check if new growth factor is approximately zero
% F3 = @(p) P.da.*trapz(P.gH.*exp(-p*P.a-P.muH_int)) - 1;
% options = optimset('TolX',1e-25);
% p0 = [-10^-3 10^-3]; % [LP,RP]
% p_new = fzero(F3,p0,options);
% P.p_hat = p_new;
    figure_setups;
    plot(P.a/365,muH_old,'b');
    hold on;
    plot(P.a/365,P.muH,'-.r');
    legend('original mortality','balanced mortality');
    keyboard
end