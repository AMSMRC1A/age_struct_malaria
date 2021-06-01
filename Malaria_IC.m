function [SH,EH,DH,AH,RH,SM,EM,IM] = Malaria_IC(NH,NM)
global a
na = length(a);

% with age structure
SH = 0.9*NH/na;
EH = 0.1*NH/na;
DH = 0*ones(na,1);
AH = 0*ones(na,1);
RH = 0*ones(na,1);

% for mosquitoes - assume at equilibrium
[SM,EM,IM] = mosquito_ODE(DH,AH,NH,NM);

end