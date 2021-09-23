function xt = human_model_der_fun(x)
% x = proportion of population: for general setting (including changing population)
% return xt for 2:end
global P
% X(alpha,group)
SH = x(:,1);
EH = x(:,2);
DH = x(:,3);
AH = x(:,4);

da = P.da;

[bH,bM] = biting_rate(1,P.gM/P.muM);

phi = 1/2;
psi = 1/2;
rho = 1/2;

Lambda_M = bM*P.K*da*trapz(exp(-P.muH_int).*(P.betaD*DH + P.betaA*AH));
Lambda_H = bH*P.betaM*(P.sigma/(P.sigma+P.muM))*((Lambda_M)/(Lambda_M + P.muM));

SHt = -Lambda_H.*SH(2:end) + phi*P.rD*DH(2:end) + P.rA*AH(2:end) - diff(SH)/da;
EHt = Lambda_H.*SH(2:end) - P.h*EH(2:end) - diff(EH)/da;
DHt = rho*P.h*EH(2:end) + psi*Lambda_H.*AH(2:end) - P.rD*DH(2:end) - diff(DH)/da;
AHt = (1-rho)*P.h*EH(2:end) - psi*Lambda_H.*AH(2:end) + (1-phi)*P.rD*DH(2:end) - P.rA*AH(2:end) - diff(AH)/da;
xt = [SHt; EHt; DHt; AHt];

end