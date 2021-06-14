function [SH,EH,DH,AH,RH,SM,EM,IM] = Malaria_IC(NH,NM)
global P
na = P.na;
da = P.da;

find_stable_age; % calculate the stable age distribution

% with age structure
SH = 0.9*P.n_tilde; %0.9*NH/na/da*ones(na,1); % cell averages
EH = 0.1*P.n_tilde; %0.1*NH/na/da*ones(na,1);
DH = 0*NH/na/da*ones(na,1);
AH = 0*NH/na/da*ones(na,1);
RH = 0*NH/na/da*ones(na,1);

% for mosquitoes - assume at equilibrium
[SM,EM,IM] = mosquito_ODE(DH,AH,NH,NM);

end

