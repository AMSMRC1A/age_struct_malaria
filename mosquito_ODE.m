function [SM,EM,IM] = mosquito_ODE(DH,AH,NH,NM)
% evolve mosquito quantities, keep at the steady state %global P
global P
[~,bM] = biting_rate(NH,NM);
lamM = FOI_M(bM,DH,AH,NH);
SM = P.gM/(lamM+P.muM);
IM = P.sigma/(P.sigma+P.muM)*lamM/(lamM+P.muM)*P.gM/P.muM;
EM = P.muM/P.sigma*IM;

end