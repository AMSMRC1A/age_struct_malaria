function [SH,EH,DH,AH,SM,EM,IM] = Malaria_IC(NH,NM)
global P

% with age structure
SH = 0.97*P.PH_stable*NH; %0.9*NH/na/da*ones(na,1); % cell averages
EH = 0.01*P.PH_stable*NH; %0.1/na/da*ones(na,1);
DH = 0.01*P.PH_stable*NH;
AH = 0.01*P.PH_stable*NH;

% for mosquitoes - assume at equilibrium
[SM,EM,IM] = mosquito_ODE(DH,AH,NH,NM);

end

