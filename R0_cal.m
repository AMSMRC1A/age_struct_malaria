function [R0,RHM,RMH] = R0_cal()
global P

% Stability of DFE when q = 0
if P.balance_fertility == 1 || P.balance_mortality == 1
    alphamax = P.age_max; %% Inf
    muH_int = P.muH_int_fun;
    [~,~,~,~,~,~,CH] = steady_state('DFE','handle'); % steady state for CH - total immunity (function handle)
    CH = @(a) CH(a)./P.PH_stable_fun(a);
    sigmoid_rho = sigmoid_prob_fun('rho'); % return a function handle
    sigmoid_phi = sigmoid_prob_fun('phi'); % return a function handle
    rho = @(a) sigmoid_rho(CH(a)); 
    phi = @(a) sigmoid_phi(CH(a)); 
%% For D part ----- double integral (alpha, a)       
    D_alpha_a = @(alpha,a) exp(-muH_int(alpha)).*exp(-P.rD.*(alpha-a)).*rho(a).*(1-exp(-P.h.*a));
    amax = @(alpha) alpha;
    D_int = integral2(D_alpha_a,0,alphamax,0,amax);
    
%%  For A part --- double + triple integral on (alpha, a, x)...
    A1_alpha_a = @(alpha,a,x) exp(-muH_int(alpha)).*exp(-P.rA.*(alpha-a)).*...
        (1-rho(a)).*(1-exp(-P.h.*a));
    A2_alpha_a_x = @(alpha,a,x) exp(-muH_int(alpha)).*exp(-P.rA.*(alpha-a)).*...
        P.rD.*(1-phi(a)).*exp(-P.rD.*(a-x)).*rho(x).*(1-exp(-P.h.*x));
    amax = @(alpha) alpha;
    xmax = @(alpha,a) a;
    A1_int = integral2(A1_alpha_a,0,alphamax,0,amax);
    A2_int = integral3(A2_alpha_a_x,0,alphamax,0,amax,0,xmax);
    A_int = A1_int+A2_int;
%% calculate R0
    [bH,bM] = biting_rate(1,P.gM/P.muM);  % assume NH=1; NH(end) for numerical simulation is > 1
    Lambda = 1/ integral(@(a) exp(-muH_int(a)), 0, alphamax);
    RHM = bH*Lambda*(P.betaD*D_int+P.betaA*A_int);
    RMH = bM*P.betaM*P.sigma/(P.sigma+P.muM)./P.muM;
    R0 = RHM*RMH;
    R0 = sqrt(R0);
end
end