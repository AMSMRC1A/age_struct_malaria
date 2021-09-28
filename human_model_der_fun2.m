function xt = human_model_der_fun2(x)
% x = number of population: for general setting (including changing population)
% return xt for 2:end
global P
% X(alpha,group)
xmat = reshape(x,[P.na,4]);
SH = xmat(:,1);
EH = xmat(:,2);
DH = xmat(:,3);
AH = xmat(:,4);

da = P.da;

[bH,bM] = biting_rate(1,P.gM/P.muM);

P.phi = sigmoid_prob(NaN(size(P.a)), 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(NaN(size(P.a)), 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(NaN(size(P.a)), 'psi'); % prob. AH -> DH

% Lambda_M = bM*P.K*da*trapz(exp(-P.muH_int).*(P.betaD*DH + P.betaA*AH));
Lambda_M = bM*da*trapz(P.betaD*DH + P.betaA*AH);
Lambda_H = bH*P.betaM*(P.sigma/(P.sigma+P.muM))*(Lambda_M/(Lambda_M + P.muM));

SHa = -Lambda_H.*SH(2:end) + P.rD*P.phi(2:end).*DH(2:end) + P.rA*AH(2:end) - diff(SH)/da - P.muH(2:end).*SH(2:end);
EHa = Lambda_H.*SH(2:end) - P.h*EH(2:end) - diff(EH)/da - P.muH(2:end).*EH(2:end);
DHa = P.rho(2:end).*P.h.*EH(2:end) + P.psi(2:end).*Lambda_H.*AH(2:end) - P.rD*DH(2:end) - diff(DH)/da - P.muH(2:end).*DH(2:end);
AHa = (1-P.rho(2:end)).*P.h.*EH(2:end) - P.psi(2:end).*Lambda_H.*AH(2:end) - P.rA*AH(2:end) + P.rD*(1-P.phi(2:end)).*DH(2:end)  - diff(AH)/da - P.muH(2:end).*AH(2:end);

% include boundary condition
birth = P.PH_stable(1);
SHa = [SH(1)-birth; SHa];
EHa = [EH(1); EHa];
DHa = [DH(1); DHa];
AHa = [AH(1); AHa];
xt = [SHa; EHa; DHa; AHa];
end