function [SH,EH,DH,AH,SM,EM,IM] = Malaria_IC(NH,NM)
global P

find_stable_age; % calculate the stable age distribution

% with age structure
SH = 0.99*P.n_tilde*NH; %0.9*NH/na/da*ones(na,1); % cell averages
EH = 0.01*P.n_tilde*NH; %0.1/na/da*ones(na,1);
DH = 0.0*P.n_tilde*NH;
AH = 0.0*P.n_tilde*NH;

% for mosquitoes - assume at equilibrium
[SM,EM,IM] = mosquito_ODE(DH,AH,NH,NM);

end

