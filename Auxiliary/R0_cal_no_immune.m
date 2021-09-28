function R0 = R0_cal_no_immune() % beta version: constant immunity parameter 
global P


% Stability of DFE when q = 0
if P.balance_fertility == 1
    alphamax = P.age_max; %% Inf
    muH_int = P.muH_int_fun;
%% For D part ----- double integral (alpha, a)
    % error = truncation on age_max & trapz integration
    % if taking age_max = inf, truncate error on age dominates; 
    % age_max = 60 method 1 worst, 10^2, method 2 ~ method 3, 10^-9
    % age_max = 100 method 1 worst, 10^-3, method 2 ~ method 3, 10^-9
    % if taking age_max = 60years, then trapz error dominates; 
    % P.da = 5 method 1 not too bad 10^-3, method 2 ~ method 3, 10^-7
    % P.da = 0.1 method 1 not too bad 10^-7, method 2 ~ method 3, 10^-7
    % use method 2
    %% method - 1 integrate on a, trapz on alpha
%     D_age = @(age) exp(-P.rD.*age).*integral(@(a) exp(P.rD.*a).*P.rho(1).*(1-exp(-P.h.*a)), 0, age);
%     D_age_exp = @(age) exp(-muH_int(age)).*D_age(age);
%     D = NaN(size(P.a));     
%     for ia = 1:length(P.a)
%         D(ia,1) = D_age_exp(P.a(ia));
%     end
%     D_int1 = trapz(D)*P.da;
    
    %% method - 2 double integral of conti function with Inf limits
    D_alpha_a = @(alpha,a) exp(-muH_int(alpha)).*exp(-P.rD.*(alpha-a)).*P.rho(1).*(1-exp(-P.h.*a));
    amax = @(alpha) alpha;
    D_int = integral2(D_alpha_a,0,alphamax,0,amax,'AbsTol',10^-10,'RelTol',10^-10);
    
    %% method - 3 reference, exact on a (n/a for non-constant immunity), then integrate alpha
%     F_P = @(a) -P.rho(1)*P.h*((P.h).*exp(-(P.rD).*a) - P.h...
%         -(P.rD).*exp(-(P.h).*a) + P.rD )./((P.h)*(P.h-P.rD)*(P.rD));
%     D_inte = integral(@(age) exp(-muH_int(age)).*F_P(age),0,Inf,'AbsTol',10^-10,'RelTol',10^-10);
% %     D_inte = P.da*trapz(exp(-P.muH_int).*F_P(0,0:P.da:P.age_max));
%     D_int1 - D_inte
%     D_int - D_inte
%     keyboard
%%  For A part --- triple integral (alpha, a, x)...
    % use Inf method 1 and 2, 10^-4
    % use age_max method 1 and 2, 10^-6
    %% method 1  double + triple integral on (alpha, a, x)
    A1_alpha_a = @(alpha,a,x) exp(-muH_int(alpha)).*exp(-P.rA.*(alpha-a)).*...
        (1-P.rho(1)).*(1-exp(-P.h.*a));
    A2_alpha_a_x = @(alpha,a,x) exp(-muH_int(alpha)).*exp(-P.rA.*(alpha-a)).*...
        P.rD.*(1-P.phi(1)).*exp(-P.rD.*(a-x)).*P.rho(1).*(1-exp(-P.h.*x));
    amax = @(alpha) alpha;
    xmax = @(alpha,a) a;
    A1_int = integral2(A1_alpha_a,0,alphamax,0,amax,'AbsTol',10^-10,'RelTol',10^-10);
    A2_int = integral3(A2_alpha_a_x,0,alphamax,0,amax,0,xmax,'AbsTol',10^-10,'RelTol',10^-10);
    A_int = A1_int+A2_int;
    
    %% method 2 - integrate on alpha, exact (a,x) - n/a for non-constant immunity
%     H_P = @(a) (((1-P.rho(1))*P.h)./(P.h)).*( ( -P.h + (P.h)*exp(-(P.rA)*a) + P.rA - (P.rA)*exp(-(P.h)*a) )./((P.rA)*(P.rA - P.h)) )...
%         + ((1-P.phi(1))*P.rD*P.rho(1)*P.h/((P.h)*(P.rD-P.h)*(P.rD))).*...
%         ( ((P.h -P.rD)/(P.rA)-(P.rD)/(P.h-P.rA)+(P.h)./(P.rD-P.rA))*exp(-(P.rA).*a) + ...
%         (P.h)*exp(-(P.rD).*a)./(P.rA-P.rD) + (P.rD-P.h)/(P.rA) + (P.rD)*exp(-(P.h).*a)/(P.h-P.rA) );
%     A_inte = integral(@(alpha) exp(-muH_int(alpha)).*H_P(alpha),0,alphamax)
%     A_int-A_inte

%% calculate R0
%     C_star = (bM*P.betaM*P.sigma./(P.sigma+P.muM)./P.muM)*(bH*P.Lambda);
%     zeta_P = @(p) C_star*da*trapz(exp(-P.muH_int).*(P.betaD.*F_P(p,0:da:age_max)' + P.betaA.*H_P(p,0:da:age_max)'));
    [bH,bM] = biting_rate(1,P.gM/P.muM);  % assume NH=1; NH(end) for numerical simulation is > 1
    RHM = bH*P.Lambda*(P.betaD*D_int+P.betaA*A_int);
    RMH = bM*P.betaM*P.sigma/(P.sigma+P.muM)./P.muM;
    R0 = RHM*RMH;
end
end