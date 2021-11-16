function [Cm,Cac,Cv,CH] = Immunity_IC_vac()
global P

na = P.na;

Cm = 0*ones(na,1); 
Cac = 0*ones(na,1);
Cv = 0*ones(na,1);
CH = P.c1*Cac+P.c2*Cm+P.c3*Cv;

P.phi = sigmoid_prob(CH, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(CH, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(CH, 'psi'); % prob. AH -> DH

end