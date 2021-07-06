function [Cm,Cac,Ctot] = Immunity_IC()
global P

na = P.na;

Cm = 0*ones(na,1); 
Cac = 0*ones(na,1);
Ctot = P.c1*Cac+P.c2*Cm;
P.phi = sigmoid_prob(Ctot, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Ctot, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Ctot, 'psi'); % prob. AH -> DH

end