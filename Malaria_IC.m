function [SH,EH,DH,AH,RH,SM,EM,IM] = Malaria_IC(NH,NM)
global a
na = length(a);

% with age structure
SH = 0.99*NH/na;
EH = 0.01*NH/na;
DH = 0*NH/na;
AH = 0*NH/na;
RH = 0*NH/na;


% for mosquitoes - assume at equilibrium
[SM,EM,IM] = mosquito_ODE(DH,AH,NH,NM);

end