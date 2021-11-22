function xt = human_model_der_prop(x)
% x = proportion of population: for general setting (including changing population)
% return xt for 2:end
global P
% X(alpha,group)
xmat = reshape(x,[P.na,5]);
SH = xmat(:,1);
EH = xmat(:,2);
DH = xmat(:,3);
AH = xmat(:,4);
Cac = xmat(:,5);

da = P.da;

[bH,bM] = biting_rate(1,P.gM/P.muM);

Lambda_M = bM*P.K*da*trapz(exp(-P.muH_int).*(P.betaD*DH + P.betaA*AH));
Lambda_H = bH*P.betaM*(P.sigma/(P.sigma+P.muM))*(Lambda_M/(Lambda_M + P.muM));

% update progression probability based on immunity Ctot
% PH_temp = SH + EH + AH + DH  -> P.PH_stable; should be the actual PH, not PH tilde
Cm0 = P.m*trapz(P.gH.*Cac.*P.PH_stable)*P.da/P.PH_stable(1);
Cm = Cm0.*exp(-P.a./P.dm);
Ctot = P.c1*Cac + P.c2*Cm;
P.phi = sigmoid_prob(Ctot, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Ctot, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Ctot, 'psi'); % prob. AH -> DH

SHa = -Lambda_H.*SH(2:end) + P.rD*P.phi(2:end).*DH(2:end) + P.rA*AH(2:end) - ...
     [(-3*SH(1)+4*SH(2)-SH(3))/(2*da); (SH(3:end)-SH(1:end-2))/(2*da)]; % diff(SH)/da;
EHa = Lambda_H.*SH(2:end) - P.h*EH(2:end) - ...
     [(-3*EH(1)+4*EH(2)-EH(3))/(2*da); (EH(3:end)-EH(1:end-2))/(2*da)]; % diff(EH)/da;
DHa = P.rho(2:end).*P.h.*EH(2:end) + P.psi(2:end).*Lambda_H.*AH(2:end) - P.rD*DH(2:end) -...
     [(-3*DH(1)+4*DH(2)-DH(3))/(2*da); (DH(3:end)-DH(1:end-2))/(2*da)]; %  diff(DH)/da;
AHa = (1-P.rho(2:end)).*P.h.*EH(2:end) - P.psi(2:end).*Lambda_H.*AH(2:end) - P.rA*AH(2:end) + P.rD*(1-P.phi(2:end)).*DH(2:end)  - ...
     [(-3*AH(1)+4*AH(2)-AH(3))/(2*da); (AH(3:end)-AH(1:end-2))/(2*da)]; % diff(AH)/da;
Caca = f(Lambda_H)*(P.cS*SH(2:end)+P.cE*EH(2:end)+P.cA*AH(2:end)+P.cD*DH(2:end)) + P.cV*P.vb(2:end).*SH(2:end) - 1/P.dac*Cac(2:end) - ...
     [(-3*Cac(1)+4*Cac(2)-Cac(3))/(2*da); (Cac(3:end)-Cac(1:end-2))/(2*da)]; % diff(Cac)/da;

Cac0 = P.cV*P.vb(1);
% include boundary condition
SHa = [-SH(1)+1; SHa];
EHa = [-EH(1); EHa];
DHa = [-DH(1); DHa];
AHa = [-AH(1); AHa];
Caca = [-Cac(1)+Cac0; Caca];

xt = [SHa; EHa; DHa; AHa; Caca];

end