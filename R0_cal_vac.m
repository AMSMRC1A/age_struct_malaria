function [R0,RHM,RMH] = R0_cal_vac()
global P

% Stability of DFE when q = 0
if P.balance_fertility == 1 || P.balance_mortality == 1
    alphamax = P.age_max; %% Inf
    muH_int = P.muH_int_fun;
    [~,~,~,~,~,~,~,~, CH] = steady_state_vac('DFE','handle'); % steady state for CH - total immunity (function handle)
    CH = @(a) CH(a)./P.PH_stable_fun(a);
    sigmoid_rho = sigmoid_prob_fun('rho'); % return a function handle
    sigmoid_phi = sigmoid_prob_fun('phi'); % return a function handle
    rho = @(a) sigmoid_rho(CH(a));
    phi = @(a) sigmoid_phi(CH(a));
    %% For D part ----- integral (alpha, a, s)
    D_alpha_a_x = @(alpha,a,x) exp(-muH_int(alpha)).*exp(-P.rD.*(alpha-a)).*rho(a).*P.h.*exp(-P.h.*(a-x)).*P.theta_fun(x);
    amax = @(alpha) alpha;
    xmax = @(alpha,a) a;
    D_int = integral3(D_alpha_a_x,0,alphamax,0,amax,0,xmax);
    %%  For A part --- integral on (alpha, a, x) + integral on (alpha, a, x, s)...
    A1_alpha_a_x = @(alpha,a,x) exp(-muH_int(alpha)).*exp(-P.rA.*(alpha-a)).*...
        (1-rho(a)).*P.h.*exp(-P.h.*(a-x)).*P.theta_fun(x);
    A2_alpha_a_x_s = @(alpha,a,x,s) exp(-muH_int(alpha)).*exp(-P.rA.*(alpha-a)).*...
        P.rD.*(1-phi(a)).*exp(-P.rD.*(a-x)).*rho(x).*P.h.*exp(-P.h.*(x-s)).*P.theta_fun(s);
    amax = @(alpha) alpha;
    xmax = @(alpha,a) a;
    smax = @(alpha,a,x) x;
    A1_int = integral3(A1_alpha_a_x,0,alphamax,0,amax,0,xmax);
    A2_int = integralN(A2_alpha_a_x_s,0,alphamax,0,amax,0,xmax,0,smax);
    A_int = A1_int+A2_int;
    %% calculate R0
    [bH,bM] = biting_rate(1,P.gM/P.muM);  % assume NH=1; NH(end) for numerical simulation is > 1
    K = 1/ integral(@(a) exp(-muH_int(a)), 0, alphamax);
    RHM = bH*K*(P.betaD*D_int+P.betaA*A_int);
    RMH = bM*P.betaM*P.sigma/(P.sigma+P.muM)./P.muM;
    R0 = RHM*RMH;
    R0 = sqrt(R0);
end
end



