function xt = human_model_der_fun(x)
% x = proportion of population: for general setting (including changing population)
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

phi = 1/2;
psi = 1/2;
rho = 1/2;

% Lambda_M = bM*P.K*da*trapz(exp(-P.muH_int).*(P.betaD*DH + P.betaA*AH));
Lambda_M = bM*da*trapz(P.betaD*DH + P.betaA*AH);
Lambda_H = bH*P.betaM*(P.sigma/(P.sigma+P.muM))*((Lambda_M)/(Lambda_M + P.muM));

SHa = -Lambda_H.*SH(2:end) + phi*P.rD*DH(2:end) + P.rA*AH(2:end) - diff(SH)/da;
EHa = Lambda_H.*SH(2:end) - P.h*EH(2:end) - diff(EH)/da;
DHa = rho*P.h*EH(2:end) + psi*Lambda_H.*AH(2:end) - P.rD*DH(2:end) - diff(DH)/da;
AHa = (1-rho)*P.h*EH(2:end) - psi*Lambda_H.*AH(2:end) - P.rA*AH(2:end) + (1-phi)*P.rD*DH(2:end)  - diff(AH)/da;

% include boundary condition
SHa = [SH(1)-1; SHa];
EHa = [EH(1); EHa];
DHa = [DH(1); DHa];
AHa = [AH(1); AHa];
xt = [SHa; EHa; DHa; AHa];

end