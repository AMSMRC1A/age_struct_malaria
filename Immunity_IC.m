function Cs = Immunity_IC()

global P
global a
na = length(a);

Cs = 0*ones(na,1); 
P.phi = sigmoid_prob(Cs(:,1), 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Cs(:,1), 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Cs(:,1), 'psi'); % prob. AH -> DH

end