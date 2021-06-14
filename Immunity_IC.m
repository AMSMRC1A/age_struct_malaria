function [Cm,Cac,Cs] = Immunity_IC()
global P

na = P.na;

Cm = 0*ones(na,1); 
Cac = 0*ones(na,1);
Cs = P.c1*Cac+P.c2*Cm;
P.phi = sigmoid_prob(Cs, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Cs, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Cs, 'psi'); % prob. AH -> DH

end