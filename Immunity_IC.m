function [Cm,Cac,Ctot] = Immunity_IC()
% Cm(na,np)
% P.phi(na,np)
global P

na = P.na;
np = P.np;

Cm = 0*ones(na,np); 
Cac = 0*ones(na,np);
Ctot = P.c1*Cac+P.c2*Cm;
P.phi = sigmoid_prob(Ctot, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Ctot, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Ctot, 'psi'); % prob. AH -> DH

end