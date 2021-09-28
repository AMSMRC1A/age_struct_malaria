function xt = human_model_der_prop_no_da(x)
% x = proportion of population: for general setting (including changing population)
% return xt for 2:end
global P
% X(alpha,group)
xmat = reshape(x,[P.na-1,4]);
SH = xmat(:,1);
EH = xmat(:,2);
DH = xmat(:,3);
AH = xmat(:,4);

da = P.da;

[bH,bM] = biting_rate(1,P.gM/P.muM);

P.phi = sigmoid_prob(NaN(P.na-1,1), 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(NaN(P.na-1,1), 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(NaN(P.na-1,1), 'psi'); % prob. AH -> DH

% Lambda_M = bM*P.K*da*trapz(exp(-P.muH_int).*(P.betaD*DH + P.betaA*AH));
Lambda_M = bM*da*trapz(P.betaD*DH + P.betaA*AH);
Lambda_H = bH*P.betaM*(P.sigma/(P.sigma+P.muM))*(Lambda_M/(Lambda_M + P.muM));

SHa = -Lambda_H.*SH + P.rD*P.phi .*DH  + P.rA*AH ;
EHa = Lambda_H.*SH - P.h*EH ;
DHa = P.rho .*P.h.*EH  + P.psi .*Lambda_H.*AH  - P.rD*DH ;
AHa = (1-P.rho ).*P.h.*EH  - P.psi .*Lambda_H.*AH  - P.rA*AH  + P.rD*(1-P.phi ).*DH ;

% include boundary condition
% SHa = [SH(1)-1; SHa];
% EHa = [EH(1); EHa];
% DHa = [DH(1); DHa];
% AHa = [AH(1); AHa];

xt = [SHa; EHa; DHa; AHa];

end