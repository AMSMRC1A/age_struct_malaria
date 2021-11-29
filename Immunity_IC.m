function [Cm,Cac,CH] = Immunity_IC()
global P

na = P.na;

Cm = 0*ones(na,1); 
Cac = 0*ones(na,1);
CH = P.c1*Cac+P.c2*Cm;
P.phi = sigmoid_prob(CH./P.PH_stable, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(CH./P.PH_stable, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(CH./P.PH_stable, 'psi'); % prob. AH -> DH

end